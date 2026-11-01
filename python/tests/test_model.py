"""Mathematical checks of the v2 model (run: python -m pytest python/tests)."""
import numpy as np
from dpm.model import Params, R0, R0_ngm, equilibrium, rhs, simulate, stable, relapse_multiplier

rng = np.random.default_rng(20261031)

def random_params(alpha_zero=False):
    return Params(b=rng.uniform(0.2, 6), eta=rng.uniform(0, 1), alpha=0.0 if alpha_zero else rng.uniform(1e-4, 0.2),
                  sigma=rng.uniform(0.05, 1), gamma1=rng.uniform(0.05, 1), gamma2=rng.uniform(0.1, 1.5),
                  delta=rng.uniform(0, 0.5), nu=rng.uniform(0.03, 0.1))

def test_R0_closed_form_equals_spectral_radius():
    for _ in range(2000):
        p = random_params(True)
        assert abs(R0(p) - R0_ngm(p)) < 1e-9 * max(1, R0(p))

def test_R0_factorises_into_no_relapse_term_and_relapse_multiplier():
    for _ in range(500):
        p = random_params(True)
        a = p.sigma + p.gamma1 + p.nu; c = p.gamma2 + p.nu
        r_nr = p.b * (p.eta * c + p.sigma) / (a * c)
        assert abs(R0(p) - r_nr * relapse_multiplier(p)) < 1e-9 * R0(p)

def test_equilibrium_is_a_fixed_point_in_the_simplex_and_stable():
    for alpha_zero in (True, False):
        for _ in range(1000):
            p = random_params(alpha_zero)
            x = equilibrium(p)
            assert np.all(x >= -1e-12) and abs(x.sum() - 1) < 1e-12
            assert max(abs(np.array(rhs(0, x, p)))) < 1e-10
            assert stable(p, x)

def test_endemic_S_equals_inverse_R0_when_alpha_zero():
    for _ in range(1000):
        p = random_params(True)
        if R0(p) > 1:
            assert abs(equilibrium(p)[0] - 1 / R0(p)) < 1e-10

def test_trajectories_stay_in_simplex_and_converge():
    for _ in range(60):
        p = random_params(rng.random() < 0.5)
        if p.alpha == 0 and abs(R0(p) - 1) < 0.3:
            continue  # critical slowing down near the threshold; covered by the stability test
        x0 = rng.dirichlet(np.ones(4))
        t, y = simulate(p, x0, 400, t_eval=np.linspace(0, 400, 401))
        assert y.min() > -1e-9 and np.allclose(y.sum(0), 1, atol=1e-8)
        assert np.allclose(y[:, -1], equilibrium(p), atol=1e-5)

def test_bifurcation_vectors_and_coefficient_sign():
    from dpm.model import jacobian
    for _ in range(300):
        p = random_params(True)
        a = p.sigma + p.gamma1 + p.nu; c = p.gamma2 + p.nu; d = p.delta + p.nu
        D = a * c * d - p.sigma * p.gamma2 * p.delta
        bstar = D / (d * (p.eta * c + p.sigma))
        q = p.with_(b=bstar)
        assert abs(R0(q) - 1) < 1e-10
        J = jacobian(q, np.array([1.0, 0, 0, 0]))
        w = np.array([(p.gamma1 - bstar * (p.eta + p.sigma / c)) / p.nu, 1.0, p.sigma / c, p.sigma * p.gamma2 / (c * d)])
        assert np.allclose(J @ w, 0, atol=1e-9)
        v4 = p.delta / d; v3 = (a - bstar * p.eta) / p.sigma
        v = np.array([0.0, 1.0, v3, v4])
        assert np.allclose(v @ J, 0, atol=1e-9)
        assert w[0] < 0   # coefficient a < 0: forward bifurcation

def test_sensitivity_of_R0_to_delta_closed_form():
    from dpm.model import elasticity
    for _ in range(200):
        p = random_params(True)
        a = p.sigma + p.gamma1 + p.nu; c = p.gamma2 + p.nu; d = p.delta + p.nu
        D = a * c * d - p.sigma * p.gamma2 * p.delta
        closed = p.sigma * p.gamma2 * p.delta * p.nu / (d * D)
        if p.delta > 1e-3:
            assert abs(elasticity(R0, p, "delta") - closed) < 1e-5

def test_floor_when_social_exposure_removed():
    for _ in range(300):
        p = random_params(False).with_(b=0.0)
        a = p.sigma + p.gamma1 + p.nu; c = p.gamma2 + p.nu; d = p.delta + p.nu
        D = a * c * d - p.sigma * p.gamma2 * p.delta
        m = 1 + p.sigma / c + p.sigma * p.gamma2 / (c * d); q = D / (c * d)
        assert abs(equilibrium(p)[2] - p.sigma / c * p.alpha / (p.alpha * m + q)) < 1e-12


def test_career_transmissibility_is_a_probability_and_reduces_without_relapse():
    from dpm import network as nw
    for _ in range(500):
        p = random_params(True)
        beta = float(np.random.uniform(1e-4, 0.3))
        Tc = nw.career_transmissibility(p, beta)
        assert 0.0 <= Tc <= 1.0
        assert Tc >= nw.transmissibility(p, beta) - 1e-12
        q = p.with_(delta=0.0)
        assert abs(nw.career_transmissibility(q, beta) - nw.transmissibility(q, beta)) < 1e-12


def test_structured_model_conserves_groups_and_has_fixed_point():
    from dpm import agesex as AS
    for _ in range(5):
        p = random_params(False)
        st = AS.Struct(phi=np.random.uniform(0.5, 3), psi=np.random.uniform(0.5, 3),
                       eps_a=np.random.uniform(0, 1), eps_s=np.random.uniform(0, 1))
        x = AS.equilibrium(p, st)
        X = x.reshape(AS.G, AS.K)
        assert np.allclose(X.sum(1), AS.N_G, atol=1e-8)
        assert np.all(X > -1e-12)
        f = AS.rhs_factory(p, st)
        assert np.max(np.abs(f(0, x))) < 1e-8


def test_structured_mixing_is_row_stochastic():
    from dpm import agesex as AS
    C = AS.mixing(AS.Struct(eps_a=0.7, eps_s=0.4))
    assert np.allclose(C.sum(1), 1.0)


def test_branching_extinction_certain_below_threshold():
    import ext_math as em
    for _ in range(50):
        p = random_params(True)
        q = em.extinction_probs(p)
        if R0(p) < 0.95:
            assert np.allclose(q, 1.0, atol=1e-6)
        elif R0(p) > 1.05:
            assert q[0] < 1 - 1e-4


def test_theorem_5_9_routh_hurwitz_at_interior_equilibrium():
    import ext_math as em
    rng = np.random.default_rng(59)
    for i in range(300):
        alpha = 0.0 if i % 2 == 0 else rng.uniform(1e-4, 0.3)
        p = Params(b=rng.uniform(0.1, 10), eta=rng.uniform(0, 1), alpha=alpha, sigma=rng.uniform(0.01, 2),
                   gamma1=rng.uniform(0.01, 2), gamma2=rng.uniform(0.05, 3), delta=rng.uniform(0, 1),
                   nu=rng.uniform(0.02, 0.2))
        if alpha == 0.0 and R0(p) <= 1.01:
            continue
        x = equilibrium(p)
        a1, a2, a3, ok = em.routh_hurwitz(em.reduced_jac(p, x))
        assert ok
        S, A, I, R = x
        lam = p.b * (I + p.eta * A) + p.alpha
        a, c, d = p.sigma + p.gamma1 + p.nu, p.gamma2 + p.nu, p.delta + p.nu
        D = a * c * d - p.sigma * p.gamma2 * p.delta
        m = 1 + p.sigma / c + p.sigma * p.gamma2 / (c * d)
        assert np.isclose(a3, D * p.alpha / lam + lam * c * d * m, rtol=1e-8)


def test_optimal_control_sweep_limits_and_adjoint_gradient():
    """Chapter 9 optimal control: huge weights give the uncontrolled path; the adjoint gives the gradient."""
    import ext_control as ec
    rng = np.random.default_rng(91)
    n = 400
    for _ in range(3):
        p = random_params()
        x0 = rng.dirichlet(np.ones(4))
        z = np.zeros(n + 1)
        PY0 = ec.trap(ec.forward(x0, z, z, p, n)[:, 2], n)
        s = ec.sweep(x0, p, 1e8, 1e8, n=n)
        assert abs(s["PY"] - PY0) < 1e-6 * PY0 and s["u1"].max() < 1e-6 and s["u2"].max() < 1e-6
        _, _, P1, P2, _ = ec.marginal(x0, p, n=n)
        eps = 1e-5
        for (c1, c2), P in (((eps, 0.0), P1), ((0.0, eps), P2)):
            fd = (PY0 - ec.prev_years(x0, p, c1, c2, n=n)[0]) / eps
            assert abs(fd - P) < 1e-3 * abs(P) + 1e-9


def test_proposition_5_10_general_durations():
    import ext_durations as ed
    rng = np.random.default_rng(510)
    for i in range(40):
        p = Params(b=rng.uniform(0.2, 6), eta=rng.uniform(0, 1), alpha=(0.0 if i % 2 == 0 else rng.uniform(1e-3, 0.1)),
                   sigma=rng.uniform(0.05, 1.5), gamma1=rng.uniform(0.05, 1.5), gamma2=rng.uniform(0.1, 2),
                   delta=rng.uniform(0.01, 0.6), nu=rng.uniform(0.03, 0.1))
        assert np.isclose(ed.R0_dur(p, 1, 1, 1), R0(p), rtol=1e-10)
        assert np.isclose(ed.I_eq_dur(p, 1, 1, 1), equilibrium(p)[2], atol=1e-12)
        k = tuple(int(v) for v in rng.integers(1, 6, size=3))
        ch = ed.Chain(*k)
        assert np.max(np.abs(ch.rhs(p)(0, ch.eq_state(p)))) < 1e-10


def test_sensitivity_indices_of_R0_all_closed_forms():
    """Proposition 1 of Paper 1: every normalised index of R0 against central differences; sign of gamma2 index."""
    names = ["b", "eta", "sigma", "gamma1", "gamma2", "delta", "nu"]
    for _ in range(300):
        p = random_params(True)
        a = p.sigma + p.gamma1 + p.nu; c = p.gamma2 + p.nu; d = p.delta + p.nu
        D = a * c * d - p.sigma * p.gamma2 * p.delta; e = p.eta * c + p.sigma
        closed = {"b": 1.0, "eta": p.eta * c / e,
                  "sigma": p.sigma * (1 / e - (c * d - p.gamma2 * p.delta) / D),
                  "gamma1": -p.gamma1 * c * d / D,
                  "gamma2": p.gamma2 * (p.eta / e - (a * d - p.sigma * p.delta) / D),
                  "delta": p.sigma * p.gamma2 * p.delta * p.nu / (d * D),
                  "nu": p.nu * (1 / d + p.eta / e - (c * d + a * d + a * c) / D)}
        for k in names:
            th = getattr(p, k)
            if th < 1e-3:
                continue
            h = 1e-6 * th
            num = (R0(p.with_(**{k: th + h})) - R0(p.with_(**{k: th - h}))) / (2 * h) * th / R0(p)
            assert abs(num - closed[k]) < 1e-5
        assert closed["gamma2"] < 0 and closed["gamma1"] < 0


def test_floor_monotone_in_rates():
    """Corollary 1 of Paper 1: the floor rises with alpha, sigma, delta and falls with gamma1, gamma2."""
    signs = {"alpha": 1, "sigma": 1, "delta": 1, "gamma1": -1, "gamma2": -1}
    for _ in range(300):
        p = random_params(False).with_(b=0.0)
        f0 = equilibrium(p)[2]
        for k, s in signs.items():
            th = getattr(p, k)
            f1 = equilibrium(p.with_(**{k: th * 1.01 + 1e-4}))[2]
            assert s * (f1 - f0) > 0
