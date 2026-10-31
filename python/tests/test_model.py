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
