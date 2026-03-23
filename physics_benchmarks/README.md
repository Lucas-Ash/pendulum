# Physics Benchmarks

This directory contains benchmark configurations centered on **specialized integrator schemes** and their adapted problems. Each scheme is designed for a specific equation structure; these benchmarks compare how well they perform on their native problems vs general-purpose methods.

## Benchmarks

### 1. Damped oscillator (DEN3-adapted)

- **Problem**: θ̈ + 2γθ̇ + ω₀²θ = g(t,θ,θ̇)
- **Specialized scheme**: DEN3 (Discrete Exponential Nyström) – solves the linear damped part exactly via matrix exponential
- **Config**: `damped_den3/damped_den3_benchmark.yaml`
- **Convergence data**: from `SerialIntegrationDampedPendulumConvergesDen3AdaptedBenchmark`
- **Plot**: `tests/artifacts/convergence/physics_benchmark_damped_den3.png`

### 2. Driven damped oscillator (DEN3-adapted)

- **Problem**: Same structure as damped, plus sinusoidal driving
- **Specialized scheme**: DEN3
- **Config**: `driven_den3/driven_den3_benchmark.yaml`
- **Convergence data**: from `SerialIntegrationDrivenPendulumConvergesDen3AdaptedBenchmark`
- **Plot**: `tests/artifacts/convergence/physics_benchmark_driven_den3.png`

### 3. Conservative pendulum (symplectic & position-only)

- **Problem**: θ̈ = -(g/L) sin θ (no damping; Hamiltonian)
- **Specialized schemes**: Symplectic (semi_implicit_euler, leapfrog, ruth4), position-only (velocity_verlet, runge_kutta_nystrom, numerov)
- **Convergence data**: from `SerialIntegrationSimplePendulumConvergesWithIncreasingResolutionAcrossIntegrators`
- **Plots**: `convergence_results.png`, `performance_by_order.png`, `physics_benchmark_conservative.png`

## Generating the plots

Run the C++ tests with convergence plotting:

```bash
./tests/run_cpp_tests.sh --plot-convergence
```

This produces all convergence CSVs and runs `tests/plot_physics_benchmarks.py` after the tests complete.
