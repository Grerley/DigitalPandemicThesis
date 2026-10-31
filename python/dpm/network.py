"""Stochastic S-A-I-R dynamics on explicit contact networks (thesis v2, Chapter 8).

Discrete monthly steps. Each node is one young person. Per month:
  * S -> A with probability 1 - (1 - alpha_m) * (1 - beta) ** (m_I + eta * m_A),
    where m_I, m_A are the numbers of Addicted and At-Risk neighbours.
    This is a probability by construction (the v1 form 1 - (1-beta)^m + mu could exceed 1).
  * A, I and R leave by competing hazards (sigma, gamma1 | gamma2 | delta), converted with
    P(leave) = 1 - exp(-h dt) and split in proportion to the hazards.
  * Turnover: with probability 1 - exp(-nu dt) a node ages out and is replaced, in the same
    network position, by a new Susceptible ten-year-old.
"""
from __future__ import annotations

import math
import numpy as np
import networkx as nx
import scipy.sparse as sp

from .model import Params

DT = 1.0 / 12.0  # one month in years


def make_graph(kind: str, n: int, k: float, seed: int, ws_p: float = 0.1):
    if kind == "ER":
        return nx.fast_gnp_random_graph(n, k / (n - 1), seed=seed)
    if kind == "WS":
        return nx.connected_watts_strogatz_graph(n, int(round(k)), ws_p, seed=seed)
    if kind == "BA":
        return nx.barabasi_albert_graph(n, int(round(k / 2)), seed=seed)
    raise ValueError(kind)


def degree_moments(G):
    d = np.array([deg for _, deg in G.degree()], float)
    return d.mean(), (d ** 2).mean()


def beta_from_b(p: Params, k: float) -> float:
    """Per-edge monthly transmission probability matching the mean-field exposure rate b = beta_rate * k."""
    rate_per_edge_per_month = p.b / k / 12.0
    return 1.0 - math.exp(-rate_per_edge_per_month)


def transmissibility(p: Params, beta: float) -> float:
    """Probability that an index case who has just entered A influences a given susceptible neighbour
    before leaving A and I (relapse and the neighbour's own state changes ignored)."""
    theta = -math.log(1.0 - beta) / DT                    # per-year hazard per edge
    a = p.sigma + p.gamma1 + p.nu
    c = p.gamma2 + p.nu
    escape_A = a / (a + p.eta * theta)                   # no transmission during the A spell
    escape_I = c / (c + theta)                            # no transmission during an I spell
    return 1.0 - escape_A * ((p.gamma1 + p.nu) / a + (p.sigma / a) * escape_I)


def R0_network(p: Params, beta: float, G) -> float:
    """Excess-degree branching number T (<k^2> - <k>) / <k> (Newman, 2002), with the relapse multiplier."""
    from .model import relapse_multiplier
    k1, k2 = degree_moments(G)
    return transmissibility(p, beta) * (k2 - k1) / k1 * relapse_multiplier(p)


def simulate(G, p: Params, beta: float, months: int, x0, seed: int, record_every: int = 1):
    rng = np.random.default_rng(seed)
    n = G.number_of_nodes()
    Adj = nx.to_scipy_sparse_array(G, format="csr", dtype=float)
    st = rng.choice(4, size=n, p=np.asarray(x0) / np.sum(x0))   # 0 S, 1 A, 2 I, 3 R
    am = 1 - math.exp(-p.alpha * DT)
    hA = p.sigma + p.gamma1 + p.nu
    pA = 1 - math.exp(-(p.sigma + p.gamma1) * DT)
    shareI = p.sigma / (p.sigma + p.gamma1)
    pI = 1 - math.exp(-p.gamma2 * DT)
    pR = 1 - math.exp(-p.delta * DT)
    pexit = 1 - math.exp(-p.nu * DT)
    log1mb = math.log(1 - beta)
    out = []
    for t in range(months + 1):
        if t % record_every == 0:
            out.append(np.bincount(st, minlength=4) / n)
        if t == months:
            break
        exposure = Adj @ ((st == 2) + p.eta * (st == 1))
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
        # turnover (independent draw)
        new[rng.random(n) < pexit] = 0
        st = new
    return np.array(out)
