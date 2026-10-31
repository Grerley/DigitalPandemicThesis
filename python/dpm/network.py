"""Stochastic S-A-I-R dynamics on explicit contact networks (thesis Chapter 8).

Discrete monthly steps. Each node is one young person. Per month:
  * S -> A with probability 1 - (1 - alpha_m) * (1 - beta) ** (m_I + eta * m_A),
    where m_I, m_A are the numbers of Addicted and At-Risk neighbours.
    This is a probability by construction.
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


def _rates(p: Params, beta: float):
    theta = -math.log(1.0 - beta) / DT                    # per-year exposure hazard along one edge
    a = p.sigma + p.gamma1 + p.nu
    c = p.gamma2 + p.nu
    d = p.delta + p.nu
    return theta, a, c, d


def transmissibility(p: Params, beta: float) -> float:
    """Single-spell transmissibility: probability that a node which has just entered A exposes a given
    susceptible neighbour successfully before its first exit from A and I (relapse ignored)."""
    theta, a, c, d = _rates(p, beta)
    escape_A = a / (a + p.eta * theta)
    escape_I = c / (c + theta)
    return 1.0 - escape_A * ((p.gamma1 + p.nu) / a + (p.sigma / a) * escape_I)


def career_transmissibility(p: Params, beta: float) -> float:
    """Probability that a node which has just entered A exposes a given susceptible neighbour successfully
    at some point in its whole career, counting returns to A through relapse (the neighbour is held
    susceptible and remission to S ends the career). Repeated spells re-expose the SAME neighbour,
    so relapse enters through the geometric sum X / (1 - Y), not through the mean-field multiplier
    1 / (1 - q), which would allow values above one.
      X = P(no transmission and the career ends without relapse)
      Y = P(no transmission and a relapse back into A)."""
    theta, a, c, d = _rates(p, beta)
    hA = a + p.eta * theta
    hI = c + theta
    X = (p.gamma1 + p.nu) / hA + (p.sigma / hA) * (p.nu / hI + (p.gamma2 / hI) * (p.nu / d))
    Y = (p.sigma / hA) * (p.gamma2 / hI) * (p.delta / d)
    return 1.0 - X / (1.0 - Y)


def R0_network(p: Params, beta: float, G) -> float:
    """Heuristic network reproduction number T_career (<k^2> - <k>) / <k>. The excess-degree form is
    exact for SIR-type processes on configuration-model networks (Newman, 2002); with remission and
    relapse the process is SIS-like, so this is an approximation to be checked by simulation."""
    k1, k2 = degree_moments(G)
    return career_transmissibility(p, beta) * (k2 - k1) / k1


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
