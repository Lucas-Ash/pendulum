# Module Structure

Modules are grouped by responsibility. Include paths use `modules/<category>/<file>`.

## Directory Layout

| Directory | Purpose | Contents |
|-----------|---------|----------|
| **core/** | Shared types, utilities, analysis | `config_utils`, `restoring_force`, `damping_force`, `error_reference`, `simulation_result`, `error_analysis`, `plotting_utils` |
| **config/** | YAML config loaders | `experiment_config` (simple), `damped_config`, `driven_config`, `coupled_config` |
| **integrators/** | Numerical integration | `rk_integrators`, `simulation_runner`, `jacobi_elliptic` |
| **simple/** | Simple (undamped) pendulum | `pendulum_simulator`, `simple_io`, `reporting`, `plotting` |
| **damped/** | Damped pendulum | `damped_simulator`, `damped_io`, `damped_plotting`, `damped_reporting` |
| **driven/** | Driven damped pendulum | `driven_simulator`, `driven_io`, `driven_plotting`, `driven_reporting`, `driven_sweep_result`, `frequency_estimation` |
| **coupled/** | Coupled oscillators | `coupled_simulator`, `coupled_state` |
| **runtime/** | Runtime environment | `runtime_env` |
| **add-ons/** | Third-party headers | `matplotlibcpp.h` |

## Workflow: YAML → Data Output

```
┌─────────────────┐     ┌──────────────────┐     ┌─────────────────┐
│ config/*.yaml   │────▶│ detect type      │────▶│ config/         │
│ (user input)   │     │ (integrated_     │     │ load_X_config()  │
└─────────────────┘     │  pendulum)       │     └────────┬────────┘
                        └──────────────────┘              │
                                                           ▼
┌─────────────────┐     ┌──────────────────┐     ┌─────────────────┐
│ data_file       │◀────│ *_io             │◀────│ *_simulator     │
│ (CSV / DAT)     │     │ write_X_data()   │     │ .simulate()      │
└─────────────────┘     └──────────────────┘     └────────┬────────┘
                                                           │
                        ┌──────────────────┐               │
                        │ integrators/     │◀───────────────┘
                        │ simulation_runner│  calls advance_state()
                        │ rk_integrators   │  (rk3, rk4, rk5, den3, etc.)
                        └──────────────────┘
```

### Step-by-step Pipeline

1. **Entry**: `integrated_pendulum <config.yaml>`
   - Reads YAML to detect type (Simple / Damped / Driven / Coupled)

2. **Config load** (from `config/`):
   - **Simple**: `load_config_from_yaml` → `ExperimentConfig`
   - **Damped**: `load_damped_config_from_yaml` → `DampedConfig`
   - **Driven**: `load_driven_config_from_yaml` → `DrivenConfig`
   - **Coupled**: `load_coupled_config_from_yaml` → `CoupledConfig`

3. **Simulate** (from `simple/`, `damped/`, `driven/`, `coupled/`):
   - Simulator builds derivatives / residuals
   - Calls `simulation_runner::run` (or equivalent) which uses `rk_integrators`
   - Returns `SimulationResult` (or `DrivenSweepResult` / `CoupledSimulationResult`)

4. **Write output** (from `*_io`):
   - **Simple**: `write_simple_data_file` → CSV
   - **Damped**: `write_damped_data_file` → DAT
   - **Driven**: `write_driven_data_file` → CSV (or sweep CSV)
   - **Coupled**: `write_coupled_data_file` → DAT

5. **Post-process** (optional):
   - `reporting`: print accuracy / summary
   - `plotting`: render via Python/matplotlib scripts

### Dependency Flow

```
core (config_utils, simulation_result, error_analysis, ...)
  ↑
config (experiment_config, damped_config, driven_config, coupled_config)
  ↑
integrators (rk_integrators, simulation_runner, jacobi_elliptic)
  ↑
simple | damped | driven | coupled (simulators, io, plotting, reporting)
```
