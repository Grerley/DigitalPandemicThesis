"""Age- and sex-structured S-A-I-R model (thesis Chapters 4, 5, 7 and 9).

Groups g = (sex, age band): sex in (girls, boys); age bands 10-12, 13-15, 16-19, 20-24 (widths 3, 3, 4, 5 years).
Young people enter at age 10 into S, half girls and half boys, move to the next band at rate 1 / width and
leave the population at 25. The stationary population share of a band is proportional to its width.

State: X[g, j] for j in (S, A, I, R), as fractions of the whole population aged 10 to 24.

For group g with population share n_g and per-capita prevalences y_h = X[h] / n_h:

    lambda_g = b * sum_h C[g, h] * (y_h[I] + eta * y_h[A])          social exposure
    inflow to A from S in g = beta_g * (lambda_g + alpha_t) * X[g, S]

beta_g = phi (girls) or 1 (boys), times psi for the youngest band (10 to 12) and 1 otherwise: a vulnerability
multiplier that scales both routes of acquisition. The other rates (sigma, gamma1, gamma2, delta, eta) are
shared by all groups. C is a row-stochastic mixing matrix:

    C[g, h] = eps_a * 1[same age] * (eps_s * 1[same sex] + (1 - eps_s) * share of h's sex in that age band)
              + (1 - eps_a) * n_h

so eps_a is the share of contacts reserved for the same age band and eps_s, within that, for the same sex.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

from .model import Params

WIDTHS = np.array([3.0, 3.0, 4.0, 5.0])
AGES = ["10-12", "13-15", "16-19", "20-24"]
SEXES = ["girls", "boys"]
G = 8                      # groups: index = sex * 4 + age
K = 4                      # states


def gid(sex: int, age: int) -> int:
    return sex * 4 + age


@dataclass(frozen=True)
class Struct:
    phi: float = 1.0       # girls' vulnerability relative to boys
    psi: float = 1.0       # vulnerability of 10-to-12-year-olds relative to older bands
    eps_a: float = 0.7     # share of contacts within the same age band
    eps_s: float = 0.5     # within the same age band, share within the same sex

    def with_(self, **kw):
        return replace(self, **kw)


def pop_shares():
    w = WIDTHS / WIDTHS.sum()
    return np.array([0.5 * w[a] for s in range(2) for a in range(4)])


N_G = pop_shares()


def mixing(st: Struct):
    C = np.zeros((G, G))
    for s in range(2):
        for a in range(4):
            g = gid(s, a)
            for s2 in range(2):
                for a2 in range(4):
                    h = gid(s2, a2)
                    same_age = st.eps_a * (a == a2) * (st.eps_s * (s == s2) + (1 - st.eps_s) * 0.5)
                    C[g, h] = same_age + (1 - st.eps_a) * N_G[h]
    return C


def vulnerability(st: Struct):
    return np.array([(st.phi if s == 0 else 1.0) * (st.psi if a == 0 else 1.0) for s in range(2) for a in range(4)])


def ageing_rates():
    return np.array([1.0 / WIDTHS[a] for s in range(2) for a in range(4)])


def rhs_factory(p: Params, st: Struct, alpha_t=None, b_t=None, access_t=None):
    """access_t(t) -> length-G vector of the share of each group with social media access (default all).
    Young people without access neither acquire problematic use nor contribute to others' social exposure."""
    C = mixing(st); beta = vulnerability(st); mu = ageing_rates()
    entry = np.zeros(G); entry[gid(0, 0)] = entry[gid(1, 0)] = 0.5 / WIDTHS.sum()   # births into the 10-12 band
    nxt = [gid(s, a + 1) if a < 3 else None for s in range(2) for a in range(4)]

    def f(t, x):
        X = x.reshape(G, K)
        al = p.alpha if alpha_t is None else alpha_t(t)
        b = p.b if b_t is None else b_t(t)
        y = X / N_G[:, None]
        m = np.ones(G) if access_t is None else access_t(t)
        lam = b * C @ (m * (y[:, 2] + p.eta * y[:, 1]))
        inc = m * beta * (lam + al) * X[:, 0]
        dX = np.empty_like(X)
        dX[:, 0] = -inc + p.gamma1 * X[:, 1]
        dX[:, 1] = inc - (p.sigma + p.gamma1) * X[:, 1] + p.delta * X[:, 3]
        dX[:, 2] = p.sigma * X[:, 1] - p.gamma2 * X[:, 2]
        dX[:, 3] = p.gamma2 * X[:, 2] - p.delta * X[:, 3]
        out = mu[:, None] * X
        dX -= out
        for g in range(G):
            if nxt[g] is not None:
                dX[nxt[g]] += out[g]
        dX[:, 0] += entry
        return dX.ravel()
    return f


def initial_state(I0=0.05):
    X = np.zeros((G, K)); X[:, 0] = N_G * (1 - I0); X[:, 2] = N_G * I0
    return X.ravel()


def simulate(p: Params, st: Struct, x0, years, t_eval=None, alpha_t=None, b_t=None, access_t=None,
             rtol=1e-8, atol=1e-11, max_step=0.25):
    if t_eval is None:
        t_eval = np.linspace(0, years, int(years) + 1)
    sol = solve_ivp(rhs_factory(p, st, alpha_t, b_t, access_t), (0, years), np.asarray(x0, float), t_eval=t_eval,
                    method="LSODA", rtol=rtol, atol=atol, max_step=max_step)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.t, sol.y


def equilibrium(p: Params, st: Struct, x0=None):
    """Integrate towards the steady state, then polish with a root-finder."""
    x0 = initial_state() if x0 is None else x0
    _, y = simulate(p, st, x0, 300, t_eval=[300], max_step=np.inf)
    x = y[:, -1]
    f = rhs_factory(p, st)
    sol = root(lambda z: f(0, z), x, method="hybr", tol=1e-13)
    return sol.x if sol.success and np.all(sol.x > -1e-12) else x


def group_prev(x):
    X = np.asarray(x).reshape(G, K)
    return X[:, 2] / N_G


def hbsc_obs(x):
    """HBSC observation operator: ages 11, 13 and 15 in equal numbers. Age 11 lies in the 10-12 band and
    ages 13 and 15 in the 13-15 band, so the HBSC figure is weighted 1/3 on the first band and 2/3 on the second;
    girls and boys are equally represented."""
    y = group_prev(x)
    band = lambda s, a: y[gid(s, a)]
    girls = band(0, 0) / 3 + 2 * band(0, 1) / 3
    boys = band(1, 0) / 3 + 2 * band(1, 1) / 3
    age11 = 0.5 * (band(0, 0) + band(1, 0))
    age13_15 = 0.5 * (band(0, 1) + band(1, 1))
    return {"total": 0.5 * (girls + boys), "girls": girls, "boys": boys, "age11": age11, "age13_15": age13_15}


def whole_prev(x):
    X = np.asarray(x).reshape(G, K)
    return X[:, 2].sum()


def R0_struct(p: Params, st: Struct):
    """Spectral radius of F V^{-1} at the problem-free state (alpha = 0). Infected variables (A, I, R) per group."""
    C = mixing(st); beta = vulnerability(st); mu = ageing_rates()
    n = 3 * G
    F = np.zeros((n, n)); V = np.zeros((n, n))
    idx = lambda g, j: 3 * g + j          # j: 0 A, 1 I, 2 R
    nxt = [gid(s, a + 1) if a < 3 else None for s in range(2) for a in range(4)]
    for g in range(G):
        for h in range(G):
            coef = p.b * beta[g] * N_G[g] * C[g, h] / N_G[h]
            F[idx(g, 0), idx(h, 0)] += coef * p.eta
            F[idx(g, 0), idx(h, 1)] += coef
        V[idx(g, 0), idx(g, 0)] = p.sigma + p.gamma1 + mu[g]
        V[idx(g, 0), idx(g, 2)] = -p.delta
        V[idx(g, 1), idx(g, 1)] = p.gamma2 + mu[g]
        V[idx(g, 1), idx(g, 0)] = -p.sigma
        V[idx(g, 2), idx(g, 2)] = p.delta + mu[g]
        V[idx(g, 2), idx(g, 1)] = -p.gamma2
        if nxt[g] is not None:
            for j in range(3):
                V[idx(nxt[g], j), idx(g, j)] -= mu[g]
    K_ = F @ np.linalg.inv(V)
    return float(max(abs(np.linalg.eigvals(K_))))
