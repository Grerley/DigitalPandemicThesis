# Model Specification

Full mathematical specification of the S–A–I–R model of digital addiction. Equation numbers match the
thesis. All rates are per month unless stated; annual rates are converted with
P = 1 − exp(−rate·Δt) (Eq 6.2).

## States
- **S** Susceptible — not problematically engaged, at risk of exposure
- **A** At-Risk — heavy/pre-clinical use, elevated risk, partially transmitting
- **I** Addicted — meets problematic-use thresholds (operational construct, not a settled diagnosis)
- **R** Recovered — reduced use; **not absorbing** — subject to relapse

## Transition matrix (Eq 4.1), state-dependent
```
                to S           to A     to I     to R
from S   [ 1 - P_SA           P_SA       0        0    ]
from A   [ P_AS        1-P_AS-P_AI      P_AI      0    ]
from I   [ 0                  0      1 - P_IR   P_IR   ]
from R   [ 0                 P_RA        0     1-P_RA  ]
```
Each row sums to 1. The relapse entry `P_RA` (mirroring δ in the ODE) is what makes the epidemic
self-sustaining; omitting it makes R absorbing and forces I → 0, which cannot reproduce the observed
endemic plateau.

## Mass-action force of infection (Eq 4.3)
```
P_SA(x_t) = 1 − exp[ −β · k · (I_t + η·A_t) ] + μ
```
Because `P_SA` depends on the current state `x_t`, the map is **nonlinear**. This is why the model is a
nonlinear difference-equation (compartmental) system, not a linear Markov chain: a fixed transition
matrix would converge monotonically to a stationary distribution and could not produce a
rise–peak–decline epidemic curve.

## State update (Eq 4.2)
```
x(t+1) = P(x_t)^T · x(t),     x(t) = [S(t), A(t), I(t), R(t)]^T
```

## Reed–Frost network transmission (Eq 4.4)
For a susceptible node with `m` addicted contacts:
```
P_i(S → A | m) = 1 − (1 − β)^m + μ
```

## Continuous-time ODE analogue (Appendix A.1)
```
dS/dt = Λ − β·S·(I + η·A) + γ1·A − μ·S
dA/dt = β·S·(I + η·A) − (σ + γ1 + μ)·A + δ·R
dI/dt = σ·A − (γ2 + μ)·I
dR/dt = γ2·I − (δ + μ)·R
```

## Basic reproduction number (Appendix A.4, next-generation matrix)
Infected subsystem {A, I}, evaluated at the disease-free equilibrium S₀:
```
F = [ β·η·S0   β·S0 ;   0   0 ]
V = [ σ+γ1+μ    0   ;  −σ   γ2+μ ]
R0 = ρ(F·V^{-1}) = β·S0·[ η·(γ2+μ) + σ ] / [ (σ+γ1+μ)·(γ2+μ) ] ≈ 2.5
```
A crude well-mixed estimate `β·k·D` (D = mean duration) gives the same order; an earlier per-capita
formula that multiplied a per-contact β by the contact number **and** divided by the recovery rate
double-counted contacts (yielding a spurious R0 ≈ 26) and is not used.

## Networked R0 (Appendix A.7)
```
R0_net = R0 · (⟨k²⟩ − ⟨k⟩) / ⟨k⟩
```
Degree heterogeneity (scale-free ⟨k²⟩ ≫ ⟨k⟩²) sharply amplifies spread and drives the epidemic
threshold toward zero.

## Equilibria and stability (Appendix A.5–A.6)
- Susceptible fraction at endemic equilibrium: `S* = S0 / R0`
- Positive endemic equilibrium exists iff `R0 > 1`, with `I* = [σ/(γ2+μ)]·A*`
- The system undergoes a **transcritical bifurcation** at `R0 = 1`: the disease-free equilibrium is
  stable for `R0 < 1` and loses stability to the endemic equilibrium for `R0 > 1`. Global stability of
  the endemic equilibrium follows from a Goh–Volterra Lyapunov function.

## Interventions (Chapter 8)
Interventions modify parameters (β for transmission-reducing structural/behavioural measures; μ for
age-gating; γ2, δ for treatment). Assumed efficacies — not estimated — so scenario outputs are
comparative, not predictive. Coverage ramps in via a logistic scale-up (Eq 8.1).
