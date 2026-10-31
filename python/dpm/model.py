"""S-A-I-R model of problematic digital use with population turnover (thesis v2, Chapter 4).

State x = (S, A, I, R) are fractions of the youth population aged 10 to 24.
All rates are per year.

    dS/dt = nu - (lam + alpha) S + gamma1 A - nu S
    dA/dt = (lam + alpha) S - (sigma + gamma1 + nu) A + delta R
    dI/dt = sigma A - (gamma2 + nu) I
    dR/dt = gamma2 I - (delta + nu) R
    lam   = b (I + eta A)

nu      turnover (entry at age 10 into S, exit at age 25 from every state)
b       effective social-exposure rate (per-contact transmission x contacts)
eta     relative influence of At-Risk compared with Addicted individuals
alpha   exogenous (non-social) acquisition rate
sigma   progression A -> I
gamma1  remission A -> S
gamma2  recovery I -> R
delta   relapse R -> A
"""
from __future__ import annotations

from dataclasses import dataclass, asdict, replace
import math
import numpy as np
from scipy.integrate import solve_ivp

STATES = ("S", "A", "I", "R")


@dataclass(frozen=True)
class Params:
    b: float
    eta: float
    alpha: float
    sigma: float
    gamma1: float
    gamma2: float
    delta: float
    nu: float

    def with_(self, **kw) -> "Params":
        return replace(self, **kw)

    def as_dict(self) -> dict:
        return asdict(self)


# ---------------------------------------------------------------- dynamics
def rhs(t, x, p: Params, alpha_t=None):
    S, A, I, R = x
    a = p.alpha if alpha_t is None else alpha_t(t)
    lam = p.b * (I + p.eta * A)
    inc = (lam + a) * S
    dS = p.nu - inc + p.gamma1 * A - p.nu * S
    dA = inc - (p.sigma + p.gamma1 + p.nu) * A + p.delta * R
    dI = p.sigma * A - (p.gamma2 + p.nu) * I
    dR = p.gamma2 * I - (p.delta + p.nu) * R
    return [dS, dA, dI, dR]


def simulate(p: Params, x0, years: float, t_eval=None, alpha_t=None, rtol=1e-10, atol=1e-12):
    """Integrate the ODE with an adaptive Runge-Kutta (Dormand-Prince) method."""
    if t_eval is None:
        t_eval = np.linspace(0, years, int(years * 12) + 1)
    sol = solve_ivp(rhs, (0, years), list(x0), t_eval=t_eval, args=(p, alpha_t),
                    method="RK45", rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.t, sol.y


# ---------------------------------------------------------------- threshold
def _blocks(p: Params):
    a = p.sigma + p.gamma1 + p.nu
    c = p.gamma2 + p.nu
    d = p.delta + p.nu
    D = a * c * d - p.sigma * p.gamma2 * p.delta     # det V
    return a, c, d, D


def R0(p: Params) -> float:
    """Closed form of rho(F V^-1) at the disease-free state (alpha = 0), infected set (A, I, R)."""
    a, c, d, D = _blocks(p)
    return p.b * d * (p.eta * c + p.sigma) / D


def R0_ngm(p: Params) -> float:
    """Spectral radius of the next-generation matrix computed numerically (check of R0)."""
    a, c, d, _ = _blocks(p)
    F = np.array([[p.b * p.eta, p.b, 0.0], [0, 0, 0], [0, 0, 0]])
    V = np.array([[a, 0, -p.delta], [-p.sigma, c, 0], [0, -p.gamma2, d]])
    return float(max(abs(np.linalg.eigvals(F @ np.linalg.inv(V)))))


def relapse_multiplier(p: Params) -> float:
    """1 / (1 - q), q = probability that an at-risk spell returns to A through I and R."""
    a, c, d, _ = _blocks(p)
    q = (p.sigma / a) * (p.gamma2 / c) * (p.delta / d)
    return 1.0 / (1.0 - q)


# ---------------------------------------------------------------- equilibria
def equilibrium(p: Params):
    """Unique equilibrium in the simplex.

    alpha = 0: disease-free (1,0,0,0) if R0 <= 1, else the endemic state with S* = 1/R0.
    alpha > 0: the unique positive root of a quadratic in A*; no disease-free state exists.
    """
    a, c, d, D = _blocks(p)
    m = 1 + p.sigma / c + p.sigma * p.gamma2 / (c * d)
    kappa = p.b * (p.sigma + p.eta * c) / c
    q = D / (c * d)
    if p.alpha == 0:
        r0 = R0(p)
        if r0 <= 1:
            return np.array([1.0, 0, 0, 0])
        A = (1 - 1 / r0) / m
    else:
        B = kappa - p.alpha * m - q
        if kappa == 0:
            A = p.alpha / (p.alpha * m + q)
        else:
            A = (B + math.sqrt(B * B + 4 * kappa * m * p.alpha)) / (2 * kappa * m)
    I = p.sigma * A / c
    R = p.sigma * p.gamma2 * A / (c * d)
    S = 1 - A - I - R
    return np.array([S, A, I, R])


def jacobian(p: Params, x):
    S, A, I, R = x
    lam = p.b * (I + p.eta * A)
    J = np.zeros((4, 4))
    J[0] = [-(lam + p.alpha) - p.nu, -p.b * p.eta * S + p.gamma1, -p.b * S, 0]
    J[1] = [lam + p.alpha, p.b * p.eta * S - (p.sigma + p.gamma1 + p.nu), p.b * S, p.delta]
    J[2] = [0, p.sigma, -(p.gamma2 + p.nu), 0]
    J[3] = [0, 0, p.gamma2, -(p.delta + p.nu)]
    return J


def stable(p: Params, x) -> bool:
    """Local stability on the invariant simplex: eigenvalues of the 3x3 reduced Jacobian."""
    J = jacobian(p, x)
    # restrict to S + A + I + R = 1 by eliminating S = 1 - A - I - R
    T = np.array([[-1, -1, -1], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
    Jr = (J @ T)[1:]
    return bool(np.all(np.linalg.eigvals(Jr).real < 0))


# ---------------------------------------------------------------- sensitivity
def elasticity(f, p: Params, name: str, h: float = 1e-6) -> float:
    """Normalised forward sensitivity (d f / d theta) * theta / f by central difference."""
    v = getattr(p, name)
    if v == 0:
        return float("nan")
    up = f(p.with_(**{name: v * (1 + h)}))
    dn = f(p.with_(**{name: v * (1 - h)}))
    return (up - dn) / (2 * h * v) * v / f(p)


def problematic(x) -> float:
    """Model counterpart of the HBSC problematic social media use indicator: A + I."""
    return float(x[1] + x[2])
