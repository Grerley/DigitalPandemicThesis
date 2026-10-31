"""Parameter-recovery study for a hidden Markov model of S-A-I-R transitions (thesis v2, Chapter 6).

Purpose: establish what a panel study would need (sample size, waves, measurement noise) for the
transition probabilities to be recoverable. The synthetic cohort is generated from KNOWN annual
transition probabilities; nothing estimated here enters the calibrated parameter set.

Observation model per wave: a symptom score (continuous proxy for a 0 to 9 Social Media Disorder Scale
count) with state-specific Gaussian emissions. S and R share the same emission distribution, so R is
separated from S only through the transition structure (R can be entered only from I).
"""
from __future__ import annotations

import numpy as np

K = 4  # S, A, I, R
# allowed transitions (row -> column); others are structural zeros
ALLOWED = np.array([[1, 1, 0, 0],
                    [1, 1, 1, 0],
                    [0, 0, 1, 1],
                    [0, 1, 0, 1]], bool)


def annual_matrix(sigma, gamma1, gamma2, delta, incidence):
    """Annual transition matrix from yearly hazards (competing risks within A)."""
    P = np.zeros((K, K))
    pSA = 1 - np.exp(-incidence)
    P[0, 1] = pSA; P[0, 0] = 1 - pSA
    hA = sigma + gamma1
    leave = 1 - np.exp(-hA)
    P[1, 2] = leave * sigma / hA; P[1, 0] = leave * gamma1 / hA; P[1, 1] = 1 - leave
    pIR = 1 - np.exp(-gamma2); P[2, 3] = pIR; P[2, 2] = 1 - pIR
    pRA = 1 - np.exp(-delta); P[3, 1] = pRA; P[3, 3] = 1 - pRA
    return P


def generate(P, pi0, mu, sd, n, waves, rng):
    states = np.zeros((n, waves), int)
    states[:, 0] = rng.choice(K, size=n, p=pi0)
    for t in range(1, waves):
        u = rng.random(n)
        cum = np.cumsum(P[states[:, t - 1]], axis=1)
        states[:, t] = (u[:, None] > cum).sum(1)
    y = rng.normal(mu[states], sd[states])
    return states, y


def _emis(y, mu, sd):
    z = (y[..., None] - mu) / sd
    return np.exp(-0.5 * z * z) / (sd * np.sqrt(2 * np.pi)) + 1e-300


def baum_welch(y, P0, pi0, mu0, sd0, tie_SR=True, iters=500, tol=1e-7):
    """Scaled forward-backward EM for many equal-length sequences y[n, T].
    Structural zeros are preserved because EM never moves mass onto a zero entry."""
    n, T = y.shape
    P, pi, mu, sd = P0.copy(), pi0.copy(), mu0.copy(), sd0.copy()
    ll_old = -np.inf
    for it in range(iters):
        B = _emis(y, mu, sd)                                  # n, T, K
        alpha = np.zeros((n, T, K)); c = np.zeros((n, T))
        alpha[:, 0] = pi * B[:, 0]; c[:, 0] = alpha[:, 0].sum(1); alpha[:, 0] /= c[:, 0, None]
        for t in range(1, T):
            alpha[:, t] = (alpha[:, t - 1] @ P) * B[:, t]
            c[:, t] = alpha[:, t].sum(1); alpha[:, t] /= c[:, t, None]
        beta = np.ones((n, T, K))
        for t in range(T - 2, -1, -1):
            beta[:, t] = (B[:, t + 1] * beta[:, t + 1]) @ P.T / c[:, t + 1, None]
        gamma = alpha * beta
        xi = np.einsum('ntj,jk,ntk->jk', alpha[:, :-1], P, B[:, 1:] * beta[:, 1:] / c[:, 1:, None])
        ll = np.log(c).sum()
        pi = gamma[:, 0].mean(0)
        P = xi / xi.sum(1, keepdims=True)
        w = gamma.reshape(-1, K); yy = y.reshape(-1)
        mu = (w * yy[:, None]).sum(0) / w.sum(0)
        var = (w * (yy[:, None] - mu) ** 2).sum(0) / w.sum(0)
        if tie_SR:  # S and R share one emission distribution
            ws = w[:, [0, 3]].sum()
            m = (w[:, 0] * yy + w[:, 3] * yy).sum() / ws
            v = (w[:, 0] * (yy - m) ** 2 + w[:, 3] * (yy - m) ** 2).sum() / ws
            mu[[0, 3]] = m; var[[0, 3]] = v
        sd = np.sqrt(var)
        if ll - ll_old < tol * abs(ll):
            break
        ll_old = ll
    return dict(P=P, pi=pi, mu=mu, sd=sd, loglik=ll, iters=it + 1)


def fit_multistart(y, P_true, pi_true, mu, sd, rng, starts=3):
    """Multi-start EM; starts are perturbations of a generic initial guess (not of the truth)."""
    best = None
    for s in range(starts):
        P0 = np.where(ALLOWED, rng.uniform(0.05, 0.3, (K, K)), 0.0)
        np.fill_diagonal(P0, 0); P0 = np.where(ALLOWED, P0, 0); np.fill_diagonal(P0, 1 - P0.sum(1))
        pi0 = rng.dirichlet(np.ones(K) * 5)
        mu0 = np.sort(rng.uniform(0, 9, K)); mu0 = np.array([mu0[0], mu0[1], mu0[3], mu0[0]])
        sd0 = np.full(K, y.std())
        fit = baum_welch(y, P0, pi0, mu0, sd0)
        if best is None or fit['loglik'] > best['loglik']:
            best = fit
    return best
