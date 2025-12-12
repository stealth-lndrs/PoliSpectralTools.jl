---
layout: default
title: "Implementation Status"
---

# Implementation Status

This page describes what remains to be implemented for the extended usage
examples and Julia tests defined in `EXAMPLES_TEST_PLAN.md`.  Contributions
should reference the cited slides inside `class_slides/` and follow the guidance
below when preparing scripts, documentation, and automated tests.

## Usage Examples

### Completed

| # | Scenario | Artifacts | Key requirements satisfied | Slide references |
| --- | --- | --- | --- | --- |
| 1 | Chebyshev linear BVP | `examples/bvp_linear.jl`, `bvp_linear_solution.*` | Variable diffusivity `-(1+x)y'' = sin(πx)` solved on Chebyshev grid, residuals benchmarked against high-N reference. | `PE_Aula_05_N.pdf` |
| 2 | Diffusion decay | `examples/diffusion.jl`, `diffusion_decay.*` | Heat equation with analytic solution; MOL + RK4 integration, GIF of decay, error reported at \(t=0.05\). | `PE_Aula_06_N.pdf` |
| 3 | Legendre vs Chebyshev | `examples/bvp_legendre.jl`, `bvp_legendre_error.png` | Legendre Lobatto solution compared to Chebyshev reference via `Generic.Bary_Interp`, log–log convergence study. | `PE_Aula_07_N.pdf`, `PE_Aula_08_N.pdf` |
| 4 | Nonlinear BVP (`y'' = \sin y`) | `examples/bvp_nonlinear.jl`, `bvp_nonlinear_solution.png` | Newton iteration with analytic derivatives, residual norm < 1e-13, iteration count recorded. | `PE_Aula_06_N.pdf` |
| 5 | Wave mixed BCs | `examples/wave_mixed_bc.jl`, `wave_mode.*` | Leapfrog with Neumann/Dirichlet boundaries, energy drift analysis, boundary flux verification. | `PE_Aula_10_N.pdf` |

### Missing

| # | Scenario | Requirements | Slide references |
| --- | --- | --- | --- |
| 6 | Wave propagation with traveling pulse | Implement `examples/wave_pulse.jl` launching a Gaussian pulse `u(x,0) = exp(-50(x+0.5)^2)` with initial velocity `u_t = -c u_x`, Dirichlet BCs on both ends, dt study showing stable vs unstable CFL choices, and optional sponge forcing. Plot snapshots and CFL error table. | `PE_Aula_10_N.pdf` (pp. 2–3) |
| 7 | Diffusion with source forcing | Create `examples/diffusion_forced.jl` solving `u_t = 0.1 u_xx + e^{-t} cos(πx)` on `[-1,1]`, Dirichlet zero walls, initial state zero. Track total “mass” (`∫ u dx`) vs injected energy and compare against a high-resolution finite difference reference. | `PE_Aula_06_N.pdf` (forcing discussion) |
| 8 | Poisson on square | Add `examples/poisson_square.jl` reproducing the manufactured solution `sin(π(x+1)/2) sin(π(y+1)/2)` with Dirichlet BCs. Report convergence when increasing `Nx = Ny`. Include surface/heatmap plots. | `PE_Aula_09_N.pdf` (Programa 16) |
| 9 | Poisson on rectangular domain | Implement `examples/poisson_rectangle.jl` on `x ∈ [-1, 2]`, `y ∈ [0, 1]` with a manufactured solution (e.g., `sin(π(x+1)/3) sinh(π y)`), emphasizing how scaling impacts `Dx²`, `Dy²`. Include anisotropic grid study. | `PE_Aula_09_N.pdf` (pp. 8–13) |
| 10 | Mapping prototype (optional) | Create `examples/mapping_bvp.jl` once `Mapping.jl` exists. Demonstrate mapping `ξ ∈ [-1,1]` to a clustered physical grid (e.g., `x = sin(πξ/2)`), solving a boundary-layer BVP, and comparing errors vs. un-mapped grid. | `PE_Aula_09_N.pdf` (pp. 13–19) |

### Deliverables per example

1. Script under `examples/` (named as suggested above).
2. Corresponding figures/GIFs saved to `docs/assets/`.
3. Documentation snippet added to `docs/USAGE_EXAMPLES.md` and the HTML gallery
   (`docs/web/examples.html`) summarizing the new scenario.

## Julia Tests

### Completed

| # | Test description | Location | What it verifies | Linked example |
| --- | --- | --- | --- | --- |
| 1 | Chebyshev BVP residual | `test/runtests.jl` – “Chebyshev BVP residual” | Interior residual < 1e-9 and BC enforcement for Example 1. | Example 1 |
| 2 | Legendre grid properties | “Legendre grid properties” | Symmetric nodes & `D₁` annihilates constants for Example 3. | Example 3 |
| 3 | Nonlinear BVP convergence | “Nonlinear BVP convergence” | Newton converges ≤ 8 iterations, residual < 1e-8. | Example 4 |
| 4 | Diffusion analytic comparison | “Diffusion analytic comparison” | MOL + RK4 error vs analytic solution + time-step refinement ratio. | Example 2 |
| 5 | Wave energy mixed BCs | “Wave energy mixed BCs” | Energy drift bound and Neumann flux agreement for Example 5. | Example 5 |

### Missing

| # | Test description | Requirements | Linked example |
| --- | --- | --- | --- |
| 6 | Traveling pulse CFL study | Extend `test/runtests.jl` to run the pulse example twice (stable vs unstable dt) and assert the stable run remains bounded while the unstable run diverges rapidly. | Usage Example 6 |
| 7 | Diffusion forcing consistency | Verify the forced diffusion script conserves net integral according to the analytic injection minus diffusion losses. Include dt refinement check. | Usage Example 7 |
| 8 | Poisson square accuracy | Ensure the manufactured square solution achieves `< 1e-6` max error at `Nx=Ny=24` and improves ≥8× when using 32 nodes. | Usage Example 8 |
| 9 | Poisson rectangular scaling | Confirm `Dx²`/`Dy²` matrices reflect domain scaling (`(2/(b-a))^2` factors) and that the rectangular solve hits `< 5e-5` max error. | Usage Example 9 |
| 10 | Mapping impact (optional) | Once mapping exists, add a regression test comparing mapped vs. unmapped errors, expecting ≥3× improvement in the boundary layer region. Skip gracefully if the mapping module is absent. | Usage Example 10 |

### Testing guidance

- Add new `@testset`s beneath the existing “Usage-driven Tests” block in
  `test/runtests.jl`.
- Keep run times under ~5 seconds per test; prefer smaller grids when a
  refinement/ratio check conveys the same message.
- Reuse helper functions (e.g., `build_grid`, `Generic.Bary_Interp`) rather than
  duplicating logic.

## Workflow Checklist

1. **Implementation** – finish the script + figures, update documentation, and
   ensure `examples/<new>.jl` runs.
2. **Testing** – extend `test/runtests.jl`, run
   `julia --project=. -e 'include("test/runtests.jl")'`.
3. **Reporting** – regenerate the full report
   (`julia --project=. scripts/run_full_report.jl`) so the new testsets appear in
   `reports/test_report_<timestamp>.md`.
4. **Docs** – update `docs/USAGE_EXAMPLES.md`, `docs/web/examples.html`, and
   link any new assets. Ensure MathJax renders formulas where needed.

Questions? Post an issue or chat on the project’s discussion board; be sure to
mention which example/test number you are tackling so efforts don’t overlap.
