"""Extended hidden Markov estimation for the Chapter 6 recovery study (workstream recx).

Generalises dpm.recovery.baum_welch in four ways, without changing that module:

1. Emission families: Gaussian score (as in dpm.recovery), binomial 0 to 9 symptom count, and a
   fully flexible categorical distribution over 0 to M-1. In every family Susceptible (S) and
   Recovered (R) share one emission distribution.
2. An optional baseline item z observed at wave 1 only (for example a question on past problematic
   use), with state-specific probabilities P(z = 1 | state) estimated freely in the M-step.
3. Frequency weights, so that identical discrete sequences are evaluated once.
4. The xi accumulation is a single matrix product, which is several times faster than the einsum
   in dpm.recovery for large panels.

The structural zeros of the transition matrix are preserved by EM: an entry that starts at zero has
zero expected count at every iteration, so it stays at zero.
"""
from __future__ import annotations

import numpy as np
from scipy.special import gammaln

from .recovery import ALLOWED, K

NMAX = 9  # items of the Social Media Disorder Scale
_LOGC = gammaln(NMAX + 1) - gammaln(np.arange(NMAX + 1) + 1) - gammaln(NMAX - np.arange(NMAX + 1) + 1)


# ---------------------------------------------------------------- data generation
def simulate_states(P, pi0, n, waves, rng):
    s = np.zeros((n, waves), int)
    s[:, 0] = rng.choice(K, size=n, p=pi0)
    cum = np.cumsum(P, axis=1)
    for t in range(1, waves):
        u = rng.random(n)
        s[:, t] = (u[:, None] > cum[s[:, t - 1]]).sum(1)
    s = np.minimum(s, K - 1)
    return s


def emit_gauss(states, mu, sd, rng):
    return rng.normal(mu[states], sd[states])


def emit_binom(states, q, rng):
    return rng.binomial(NMAX, q[states])


def emit_betabinom(states, q, rho, rng):
    """Beta-binomial count with state mean 9 q and intra-person item correlation rho."""
    a = q * (1 - rho) / rho
    b = (1 - q) * (1 - rho) / rho
    p = rng.beta(a[states], b[states])
    return rng.binomial(NMAX, p)


def baseline_item(states0, h, rng):
    """Binary baseline item at wave 1 with P(z = 1 | state) = h[state]."""
    return (rng.random(states0.shape[0]) < h[states0]).astype(int)


def pmf_binom(q):
    k = np.arange(NMAX + 1)
    return np.exp(_LOGC + k * np.log(q[:, None]) + (NMAX - k) * np.log1p(-q[:, None]))


def pmf_betabinom(q, rho):
    k = np.arange(NMAX + 1)
    a = q * (1 - rho) / rho
    b = (1 - q) * (1 - rho) / rho
    lb = gammaln(k + a[:, None]) + gammaln(NMAX - k + b[:, None]) - gammaln(NMAX + a[:, None] + b[:, None]) \
        - (gammaln(a) + gammaln(b) - gammaln(a + b))[:, None]
    return np.exp(_LOGC + lb)


# ---------------------------------------------------------------- compression
def compress(y, z=None):
    """Unique rows of a discrete panel (with the baseline item appended) and their counts."""
    X = y if z is None else np.column_stack([y, z])
    U, inv, cnt = np.unique(X, axis=0, return_inverse=True, return_counts=True)
    if z is None:
        return U, None, cnt.astype(float)
    return U[:, :-1], U[:, -1], cnt.astype(float)


# ---------------------------------------------------------------- emission evaluation
def _emission(y, fam, th):
    if fam == "gauss":
        z = (y[..., None] - th["mu"]) / th["sd"]
        return np.exp(-0.5 * z * z) / (th["sd"] * np.sqrt(2 * np.pi)) + 1e-300
    if fam == "binom":
        return pmf_binom(th["q"]).T[y] + 1e-300          # y integer, table lookup
    if fam == "cat":
        return th["e"].T[y] + 1e-300                     # e[K, M]
    raise ValueError(fam)


def _m_emission(y, w_obs, fam, th):
    """M-step for emissions with S and R tied. w_obs[n, T, K] are posterior weights times frequency."""
    W = w_obs.reshape(-1, K)
    yy = y.reshape(-1)
    Wt = W.copy(); Wt[:, 0] += Wt[:, 3]; Wt[:, 3] = Wt[:, 0]     # pooled S and R weights
    tot = np.ones(len(yy)) @ Wt
    if fam == "gauss":
        yf = yy.astype(float)
        mu = (yf @ Wt) / tot
        var = ((yf * yf) @ Wt) / tot - mu * mu
        return dict(mu=mu, sd=np.sqrt(np.maximum(var, 1e-4)))
    if fam == "binom":
        q = (yy.astype(float) @ Wt) / (NMAX * tot)
        return dict(q=np.clip(q, 1e-6, 1 - 1e-6))
    if fam == "cat":
        M = th["e"].shape[1]
        e = np.zeros((K, M))
        for k in range(K):
            e[k] = np.bincount(yy, weights=Wt[:, k], minlength=M)
        e = np.maximum(e, 1e-10)
        return dict(e=e / e.sum(1, keepdims=True))
    raise ValueError(fam)


# ---------------------------------------------------------------- EM
def _step(y, fam, prm, z, w, h_fix=None):
    """One EM update. Returns the updated parameters and the log-likelihood at the input parameters."""
    n, T = y.shape
    P, pi, th, h = prm["P"], prm["pi"], prm["th"], prm["h"]
    B = _emission(y, fam, th)                                   # n, T, K
    if z is not None:
        B = B.copy(); B[:, 0] = B[:, 0] * np.where(z[:, None] == 1, h, 1 - h) + 1e-300
    one = np.ones(K)
    al = np.empty((n, T, K)); c = np.empty((n, T))
    a = pi * B[:, 0]; c[:, 0] = a @ one; al[:, 0] = a / c[:, 0, None]
    for t in range(1, T):                                       # scaled forward recursion
        a = (al[:, t - 1] @ P) * B[:, t]
        c[:, t] = a @ one; al[:, t] = a / c[:, t, None]
    be = np.empty((n, T, K)); be[:, -1] = 1.0
    for t in range(T - 2, -1, -1):                              # scaled backward recursion
        be[:, t] = (B[:, t + 1] * be[:, t + 1]) @ P.T / c[:, t + 1, None]
    ll = float((w[:, None] * np.log(c)).sum())
    wg = al * be * w[:, None, None]                             # posterior state weights
    L = (al[:, :-1] * w[:, None, None]).reshape(-1, K)
    Rr = (B[:, 1:] * be[:, 1:] / c[:, 1:, None]).reshape(-1, K)
    xi = P * (L.T @ Rr)                                         # expected transition counts
    out = dict(P=xi / xi.sum(1, keepdims=True), pi=wg[:, 0].sum(0) / w.sum(),
               th=_m_emission(y, wg, fam, th), h=None)
    if z is not None:
        h = (wg[:, 0] * z[:, None]).sum(0) / wg[:, 0].sum(0)
        if h_fix is not None:                                   # known operating characteristics
            for k, v in h_fix.items():
                h[k] = v
        out["h"] = np.clip(h, 1e-8, 1 - 1e-8)
    return out, ll


def _pack(prm):
    parts = [prm["P"].ravel(), prm["pi"]] + [prm["th"][k].ravel() for k in sorted(prm["th"])]
    if prm["h"] is not None:
        parts.append(prm["h"])
    return np.concatenate(parts)


def _unpack(v, like):
    out = {}; i = 0
    out["P"] = v[i:i + K * K].reshape(K, K); i += K * K
    out["pi"] = v[i:i + K]; i += K
    th = {}
    for k in sorted(like["th"]):
        m = like["th"][k].size
        th[k] = v[i:i + m].reshape(like["th"][k].shape); i += m
    out["th"] = th
    out["h"] = None if like["h"] is None else v[i:i + K]
    return out


def _valid(prm, fam):
    if prm["P"].min() < 0 or prm["pi"].min() < 0:
        return False
    if prm["h"] is not None and (prm["h"].min() <= 0 or prm["h"].max() >= 1):
        return False
    th = prm["th"]
    if fam == "gauss":
        return th["sd"].min() > 1e-2
    if fam == "binom":
        return th["q"].min() > 0 and th["q"].max() < 1
    return th["e"].min() > 0


def em(y, fam, P0, pi0, th0, z=None, h0=None, wts=None, iters=1000, tol=1e-8, rule="abs", accel=False, h_fix=None):
    """Baum-Welch for equal-length sequences y[n, T] with S/R-tied emissions and structural zeros.

    z: optional baseline item at wave 1 (n,), with initial P(z=1|state) h0 (K,).
    h_fix: optional {state index: P(z=1|state)} held fixed (known sensitivity and specificity).
    wts: frequency weights (n,).
    rule "abs": stop when the log-likelihood gain per sequence falls below tol (strict; used for the
    extended study). rule "rel": stop when the gain falls below tol times |log-likelihood|, the rule of
    dpm.recovery.baum_welch (with tol = 1e-7 there).
    accel: SQUAREM extrapolation (Varadhan and Roland, 2008, scheme S3) with a monotonicity safeguard.
    Each accelerated cycle costs three E-steps; `iters` counts E-steps in both modes.
    """
    n, T = y.shape
    w = np.ones(n) if wts is None else wts
    N = w.sum()
    prm = dict(P=P0.copy(), pi=pi0.copy(), th={k: v.copy() for k, v in th0.items()},
               h=None if z is None else h0.copy())
    if h_fix is not None:
        for k, v in h_fix.items():
            prm["h"][k] = min(max(v, 1e-8), 1 - 1e-8)
    ll_old = -np.inf; nE = 0
    while nE < iters:
        if not accel:
            new, ll = _step(y, fam, prm, z, w, h_fix); nE += 1
        else:
            p1, ll = _step(y, fam, prm, z, w, h_fix)
            p2, _ = _step(y, fam, p1, z, w, h_fix); nE += 2
            x0, x1, x2 = _pack(prm), _pack(p1), _pack(p2)
            r = x1 - x0; v = x2 - x1 - r
            nv = np.linalg.norm(v)
            cand = p2
            if nv > 0:
                alpha = -np.linalg.norm(r) / nv
                while alpha < -1:
                    xp = x0 - 2 * alpha * r + alpha ** 2 * v
                    pp = _unpack(xp, prm)
                    if _valid(pp, fam):
                        pp["P"] = pp["P"] / pp["P"].sum(1, keepdims=True); pp["pi"] = pp["pi"] / pp["pi"].sum()
                        if fam == "cat":
                            pp["th"]["e"] = pp["th"]["e"] / pp["th"]["e"].sum(1, keepdims=True)
                        cand = pp
                        break
                    alpha = (alpha - 1) / 2
            new, ll_c = _step(y, fam, cand, z, w, h_fix); nE += 1     # stabilising EM step
            if ll_c < ll:                                       # safeguard: fall back to plain EM
                new = p2
        gain = ll - ll_old
        if (rule == "abs" and gain / N < tol) or (rule == "rel" and gain < tol * abs(ll)):
            break
        ll_old = ll
        prm = new
    return dict(P=new["P"], pi=new["pi"], th=new["th"], h=new["h"], loglik=ll, iters=nE)


def generic_start(fam, rng, y=None, M=NMAX + 1):
    """Starting values that do not depend on the generating truth (as in dpm.recovery.fit_multistart)."""
    P0 = np.where(ALLOWED, rng.uniform(0.05, 0.3, (K, K)), 0.0)
    np.fill_diagonal(P0, 0); np.fill_diagonal(P0, 1 - P0.sum(1))
    pi0 = rng.dirichlet(np.ones(K) * 5)
    if fam == "gauss":
        m = np.sort(rng.uniform(0, 9, K))
        th = dict(mu=np.array([m[0], m[1], m[3], m[0]]), sd=np.full(K, y.std()))
    elif fam == "binom":
        qq = np.sort(rng.uniform(0.02, 0.9, 3))
        th = dict(q=np.array([qq[0], qq[1], qq[2], qq[0]]))
    elif fam == "cat":
        if M == NMAX + 1:
            qq = np.sort(rng.uniform(0.02, 0.9, 3))
            e = 0.8 * pmf_binom(np.array([qq[0], qq[1], qq[2], qq[0]])) + 0.2 / M
        else:  # coarse categories: increasing probability of the top category by state severity
            base = np.sort(rng.uniform(0.02, 0.9, 3))
            lev = np.array([base[0], base[1], base[2], base[0]])
            e = np.stack([np.linspace(1 - l, l, M) for l in lev]) + 0.05
        th = dict(e=e / e.sum(1, keepdims=True))
    else:
        raise ValueError(fam)
    return P0, pi0, th


def fit(y, fam, rng, starts=2, z=None, wts=None, iters=1000, tol=1e-8, M=NMAX + 1, rule="abs", accel=False,
        h_fix=None):
    """Multi-start EM from generic starting values; keeps the fit with the highest likelihood."""
    best = None
    for _ in range(starts):
        P0, pi0, th0 = generic_start(fam, rng, y, M)
        h0 = None if z is None else rng.uniform(0.2, 0.8, K)
        f = em(y, fam, P0, pi0, th0, z=z, h0=h0, wts=wts, iters=iters, tol=tol, rule=rule, accel=accel,
                 h_fix=h_fix)
        if best is None or f["loglik"] > best["loglik"]:
            best = f
    return best


# ---------------------------------------------------------------- continuous-time fit
def _ct_unpack(x):
    hz = np.exp(x[:5])                                   # incidence, gamma1, sigma, gamma2, delta
    a = np.concatenate([[0.0], x[5:8]]); pi = np.exp(a - a.max()); pi /= pi.sum()
    mu = np.array([x[8], x[9], x[10], x[8]]); sd = np.exp(np.array([x[11], x[12], x[13], x[11]]))
    return hz, pi, mu, sd


def ct_matrix(hz):
    """Exact one-year transition matrix exp(Q) of the continuous-time chain with yearly hazards hz."""
    from scipy.linalg import expm
    lam, g1, sg, g2, de = hz
    Q = np.zeros((K, K))
    Q[0, 1] = lam; Q[1, 0] = g1; Q[1, 2] = sg; Q[2, 3] = g2; Q[3, 1] = de
    np.fill_diagonal(Q, -Q.sum(1))
    return expm(Q)


def _ct_nll(x, y):
    hz, pi, mu, sd = _ct_unpack(x)
    P = ct_matrix(hz)
    B = _emission(y, "gauss", dict(mu=mu, sd=sd))
    one = np.ones(K)
    a = pi * B[:, 0]; c = a @ one; ll = np.log(c).sum(); a = a / c[:, None]
    for t in range(1, y.shape[1]):
        a = (a @ P) * B[:, t]; c = a @ one; ll += np.log(c).sum(); a = a / c[:, None]
    return -ll


def fit_ct(y, start):
    """Maximum likelihood for the continuous-time hidden Markov model with Gaussian emissions tied for S and R.
    The transition matrix is exp(Q) with the generator of Chapter 4 (no turnover within the panel), so
    multiple transitions within a year are allowed. `start` is a one-step fit (from em) whose annual
    probabilities are converted to hazards; the optimiser is L-BFGS-B with numerical gradients."""
    from scipy.optimize import minimize
    Pm = start["P"]
    lam = -np.log(1 - Pm[0, 1]); g2 = -np.log(1 - Pm[2, 3]); de = -np.log(1 - Pm[3, 1])
    leave = Pm[1, 2] + Pm[1, 0]; hA = -np.log(1 - leave)
    hz0 = np.array([lam, hA * Pm[1, 0] / leave, hA * Pm[1, 2] / leave, g2, de])
    pi0 = np.clip(start["pi"], 1e-4, None); pi0 /= pi0.sum()
    th = start["th"]
    x0 = np.concatenate([np.log(np.clip(hz0, 1e-4, None)), np.log(pi0[1:] / pi0[0]),
                         th["mu"][:3], np.log(th["sd"][:3])])
    r = minimize(_ct_nll, x0, args=(y,), method="L-BFGS-B", options=dict(maxiter=2000))
    hz, pi, mu, sd = _ct_unpack(r.x)
    return dict(hz=hz, pi=pi, mu=mu, sd=sd, loglik=-r.fun, converged=bool(r.success), nit=r.nit)
