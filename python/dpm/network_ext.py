"""Extensions of the network model for the deepened Chapter 8 (used by ext_networks.py).

simulate_ext generalises dpm.network.simulate (same monthly update, same order of draws within a step) to:
  * node-specific exogenous acquisition rates alpha_i (homophily without contagion),
  * node-specific susceptibility multipliers w_i on social exposure (targeting),
  * a complex-contagion rule: social exposure acts on a Susceptible node only when at least
    `min_affected` of its neighbours are At-Risk or Addicted,
  * explicit initial states and snapshots of the full state vector, with a person identifier that
    changes when a node ages out and is replaced (so that two "survey waves" follow the same people).
It does not change dpm.network.
"""
from __future__ import annotations

import math
import numpy as np
import networkx as nx

from .model import Params

DT = 1.0 / 12.0


def simulate_ext(G, p: Params, beta: float, months: int, x0=None, seed: int = 0, *, init_states=None,
                 alpha_vec=None, sigma_vec=None, w=None, min_affected: int | None = None, record_every: int = 12,
                 snapshots=(), Adj=None):
    """Monthly stochastic S-A-I-R run on graph G.

    Returns dict with 'traj' (rows: recorded months, columns: S, A, I, R fractions), 'months' and
    'snaps' {month: (states, person_id)}."""
    rng = np.random.default_rng(seed)
    n = G.number_of_nodes()
    if Adj is None:
        Adj = nx.to_scipy_sparse_array(G, format="csr", dtype=float)
    if init_states is not None:
        st = np.asarray(init_states, dtype=np.int64).copy()
    else:
        st = rng.choice(4, size=n, p=np.asarray(x0) / np.sum(x0))
    pid = np.arange(n)
    next_id = n
    if alpha_vec is None:
        am = 1 - math.exp(-p.alpha * DT)
    else:
        am = 1 - np.exp(-np.asarray(alpha_vec) * DT)
    if sigma_vec is None:
        pA = 1 - math.exp(-(p.sigma + p.gamma1) * DT)
        shareI = p.sigma / (p.sigma + p.gamma1)
    else:
        sv = np.asarray(sigma_vec, float)
        pA = 1 - np.exp(-(sv + p.gamma1) * DT)
        shareI = sv / (sv + p.gamma1)
    pI = 1 - math.exp(-p.gamma2 * DT)
    pR = 1 - math.exp(-p.delta * DT)
    pexit = 1 - math.exp(-p.nu * DT)
    log1mb = math.log(1 - beta)
    snaps = {}
    snapset = set(snapshots)
    out, rec_m = [], []
    for t in range(months + 1):
        if t % record_every == 0:
            out.append(np.bincount(st, minlength=4) / n); rec_m.append(t)
        if t in snapset:
            snaps[t] = (st.copy(), pid.copy())
        if t == months:
            break
        Ivec = (st == 2).astype(float); Avec = (st == 1).astype(float)
        exposure = Adj @ (Ivec + p.eta * Avec)
        if min_affected is not None:
            m_aff = Adj @ (Ivec + Avec)
            exposure = np.where(m_aff >= min_affected, exposure, 0.0)
        if w is not None:
            exposure = exposure * w
        p_inf = 1 - (1 - am) * np.exp(log1mb * exposure)
        u = rng.random(n)
        new = st.copy()
        S = st == 0; A = st == 1; I = st == 2; R = st == 3
        new[S & (u < p_inf)] = 1
        leaveA = A & (u < pA)
        new[leaveA & (u < pA * shareI)] = 2
        new[leaveA & (u >= pA * shareI)] = 0
        new[I & (u < pI)] = 3
        new[R & (u < pR)] = 1
        out_ = rng.random(n) < pexit
        new[out_] = 0
        k = int(out_.sum())
        if k:
            pid[out_] = np.arange(next_id, next_id + k); next_id += k
        st = new
    return {"traj": np.array(out), "months": np.array(rec_m), "snaps": snaps}


def stationary_b0(p: Params, alpha, sigma=None):
    """Stationary state distribution (S, A, I, R) of one isolated node (b = 0) with exogenous rate alpha and
    progression rate sigma (default p.sigma), including turnover. Vectorised. This is the b = 0 equilibrium
    of the ODE (Theorem 5.8)."""
    alpha = np.asarray(alpha, float)
    sg = p.sigma if sigma is None else np.asarray(sigma, float)
    a = sg + p.gamma1 + p.nu; c = p.gamma2 + p.nu; d = p.delta + p.nu
    m = 1 + sg / c + sg * p.gamma2 / (c * d)
    omega = (a * c * d - sg * p.gamma2 * p.delta) / (c * d)
    A = alpha / (alpha * m + omega)
    I = sg * A / c
    R = sg * p.gamma2 * A / (c * d)
    S = 1 - A - I - R
    return np.stack([S, A, I, R], axis=-1)


def homophilous_graph(z, k: float, h: float, groups: int, seed: int):
    """Stochastic block model on trait groups: nodes are split into `groups` equal groups by the rank of the
    latent trait z; a share h of each node's expected ties lie within its own group and 1 - h are spread
    uniformly over all nodes. Mean degree k."""
    n = len(z)
    rng = np.random.default_rng(seed)
    order = np.argsort(z)
    g = np.empty(n, int); g[order] = np.arange(n) * groups // n
    size = n / groups
    p_in = h * k / (size - 1) + (1 - h) * k / (n - 1)
    p_out = (1 - h) * k / (n - 1)
    edges = []
    iu = np.triu_indices(n, 1)
    # draw all pairs at once (n = 2000 gives about two million pairs)
    same = g[iu[0]] == g[iu[1]]
    pr = np.where(same, p_in, p_out)
    keep = rng.random(len(pr)) < pr
    G = nx.Graph(); G.add_nodes_from(range(n))
    G.add_edges_from(zip(iu[0][keep].tolist(), iu[1][keep].tolist()))
    return G


def acquaintance_nodes(G, m: int, rng):
    """Acquaintance selection (Cohen, Havlin and ben-Avraham, 2003): pick a node at random, then protect one of
    its neighbours chosen at random; repeat until m distinct nodes are protected. Needs no degree data."""
    nodes = np.arange(G.number_of_nodes())
    nbrs = [list(G.neighbors(i)) for i in nodes]
    chosen = set()
    while len(chosen) < m:
        i = int(rng.integers(len(nodes)))
        if not nbrs[i]:
            continue
        chosen.add(int(nbrs[i][rng.integers(len(nbrs[i]))]))
    return np.fromiter(chosen, int)


def logit_fit(X, y, iters: int = 50):
    """Logistic regression by iteratively reweighted least squares. X includes the intercept column.
    Returns coefficients and model-based standard errors."""
    X = np.asarray(X, float); y = np.asarray(y, float)
    beta = np.zeros(X.shape[1])
    for _ in range(iters):
        eta = X @ beta
        mu = 1 / (1 + np.exp(-eta))
        W = mu * (1 - mu)
        H = X.T @ (X * W[:, None])
        g = X.T @ (y - mu)
        step = np.linalg.solve(H, g)
        beta = beta + step
        if np.max(np.abs(step)) < 1e-10:
            break
    mu = 1 / (1 + np.exp(-(X @ beta)))
    H = X.T @ (X * (mu * (1 - mu))[:, None])
    se = np.sqrt(np.diag(np.linalg.inv(H)))
    return beta, se


def stationary_discrete_b0(p: Params, alpha, sigma=None):
    """Exact stationary distribution of the monthly chain that simulate_ext runs for an isolated node
    (b = 0), with node-specific alpha and sigma. Used to start runs at stationarity. Vectorised."""
    alpha = np.atleast_1d(np.asarray(alpha, float))
    sg = np.full_like(alpha, p.sigma) if sigma is None else np.broadcast_to(np.asarray(sigma, float), alpha.shape)
    n = len(alpha)
    am = 1 - np.exp(-alpha * DT)
    pA = 1 - np.exp(-(sg + p.gamma1) * DT); sh = sg / (sg + p.gamma1)
    pI = 1 - math.exp(-p.gamma2 * DT); pR = 1 - math.exp(-p.delta * DT); px = 1 - math.exp(-p.nu * DT)
    M = np.zeros((n, 4, 4))
    M[:, 0, 0] = 1 - am; M[:, 0, 1] = am
    M[:, 1, 0] = pA * (1 - sh); M[:, 1, 1] = 1 - pA; M[:, 1, 2] = pA * sh
    M[:, 2, 2] = 1 - pI; M[:, 2, 3] = pI
    M[:, 3, 1] = pR; M[:, 3, 3] = 1 - pR
    P = (1 - px) * M; P[:, :, 0] += px
    A = np.transpose(P, (0, 2, 1)) - np.eye(4)[None]
    A[:, 3, :] = 1.0
    rhs = np.zeros((n, 4, 1)); rhs[:, 3, 0] = 1.0
    return np.linalg.solve(A, rhs)[:, :, 0]
