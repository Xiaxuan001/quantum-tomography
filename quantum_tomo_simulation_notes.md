# Quantum tomography simulation implementation notes (CVXR1121.R aligned)

This file specifies the simulation setup to match the functions in `CVXR1121.R`.
It is written to be directly actionable in R.

## 1. Scope and deliverables

We simulate adaptive measurement design using the SU(N) Bloch-basis
parameterization implemented in `CVXR1121.R` (default N=2 for a qubit).

We compare three policy families:

- Uniform (sequential cycle over settings)
- Adaptive exact (one-step lookahead using D or A optimality)
- Adaptive first-order (D or A first-order proxy)

Deliverables in the current code are PNG plots of:

- MLE distance vs step
- D-optimal and A-optimal criteria vs cumulative samples

## 2. Bloch-basis parameterization

Use the Kimura SU(N) basis built by `build_suN_basis()`:

- Basis matrices `sigma_j` satisfy `Tr(sigma_i sigma_j) = 2 delta_ij`
- Dimension `d = N^2 - 1`

Parameterization:

- `rho(theta) = I/N + 0.5 * sum_j theta_j sigma_j`

For N=2, `theta` equals the usual Bloch vector `r`.

Constraint: positivity is enforced by the PSD constraint in the CVXR MLE;
there is no explicit eigenvalue floor or ball constraint in `CVXR1121.R`.

## 3. Measurement settings and outcomes

Settings are built directly from the basis:

- For each `a = 1..d`, set `B_a = sigma_a`
- Compute spectral projectors of `B_a` with `spectral_projectors()`
- Each setting has `r_a` outcomes (number of distinct eigenvalues)

Use `build_measurements_from_basis()` to obtain:

- `Q_list[[a]][[b]]`: projector for outcome `b` of setting `a`
- `ab_row` and `ab_df`: mapping from (a,b) to a row index

For N=2, each setting has two outcomes (Pauli PVM).

## 4. Born probabilities and affine model

Born probabilities:

- `p_{a,b}(theta) = Tr(rho(theta) Q_{a,b})`

Affine representation used in the code:

- `p = c_ab + 0.5 * S_ab %*% theta`
- `S_ab[m,j] = Tr(sigma_j Q_{a,b})`
- `c_ab[m] = Tr(Q_{a,b}) / N`

Use `build_Sab_cab()` to precompute `S_ab` and `c_ab`.

## 5. Counts and schedules

Maintain a single counts vector `Nab` over all (a,b) rows:

- Update with `increment_ab_count()`
- Aggregate per-setting counts with `counts_by_setting()`

Uniform schedule:

- `sequential_actions(d, n)` cycles settings `1..d` deterministically

Sampling outcomes:

- Use `born_probs_list()` for the chosen state
- Sample `b` with `sample.int(length(probs_a), prob = probs_a)`

## 6. CVXR MLE for theta

MLE is solved by `fit_theta_cvxr()`:

- Variables: `theta` and real-embedded `rho` (PSD)
- Constraints: `rho == real_embed(rho(theta))`, `p >= 0`
- Objective: `-sum(Nab * log(p + eps_log))`

Solver fallback chain: MOSEK -> ECOS -> SCS.
Use `eps_log = 1e-12` for numerical stability.

## 7. Fisher information and criteria

Single-setting Fisher information (Eq. 6.4 in code):

- `I_a(theta) = 0.25 * t(S_a) * diag(1 / p_a) * S_a`

Total information:

- `I_total = sum_a n_a I_a`

Design criteria:

- D-optimal: `logdet(I_total + I_a)`
- A-optimal: `-tr((I_total + I_a)^{-1})`

First-order proxies used in code:

- D: `tr(Minv * I_a)` where `Minv = solve(I_total + ridge * I)`
- A: `tr(Minv %*% I_a %*% Minv)`

Use `safe_logdet()` and `safe_trace_inverse()` with a small ridge for stability.

## 8. Policy selection

`select_next_setting()` implements:

- `criterion = "D"` or `"A"`
- `method = "exact"` or `"first-order"`

If exact scoring is not finite, the function falls back to first-order.

## 9. Simulation flow (aligned to CVXR1121.R)

Typical flow:

1. `run_demo()` builds the basis, measurements, affine model, and a true state.
2. Use `sequential_actions()` for an initial block of samples per setting.
3. Run one of:
   - `uniform_design_sequence()` for uniform/sequential sampling
   - `adaptive_design_sequence()` for exact or first-order D/A selection
4. Repeat with `average_strategy_performance()` to average curves over replications.

The code uses `random_density()` or user-specified `rho_true` as the ground truth.

## 10. Plotting and outputs

Plots produced by existing helpers:

- `plot_mle_distance_curves()` or `plot_avg_distance_curves()` -> `mle_distance_curves.png`
- `plot_criterion_vs_samples()` or `plot_avg_criterion_curves()` -> `d_optimal_vs_samples.png` and `a_optimal_vs_samples.png`

These plots use cumulative samples on the x-axis and criterion values or MLE
distance on the y-axis.

## 11. Reproducibility

- Fix RNG seeds for random state generation and sampling.
- Keep `N`, `n_per_d`, `steps_adapt`, `solver`, and `ridge` fixed across runs.
