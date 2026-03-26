# Pendulum Simulation Configuration Guide

This guide details the available parameters and configuration options for writing a simulation deck (`.yaml` file) in the Integrated Pendulum Simulator. 

The simulator auto-detects the physics engine required based on the presence of specific keywords in your YAML file. Depending on the simulation environment (**Simple**, **Damped**, **Driven**, or **Coupled**), different physical and simulation parameters become available. 

---

## 1. Simulation Environments Auto-Detection

When you execute:
```bash
./integrated_pendulum config.yaml
```
The program parses the YAML file to classify your simulation into one of four systems:

1. **Coupled** System: Detected if `omega_1`, `alpha_11`, or `q1_0` are present.
2. **Driven** System: Detected if `omega_drive`, `drive_force`, `sweep`, or `system_model: duffing` are present.
3. **Damped** System: Detected if `gamma`, `damping`, `damping_model`, `damping_linear`, or `damping_cubic` are present.
4. **Simple** System: Detected if none of the above keywords are found.

You can organize your configuration keys logically via namespaces (e.g. `physical.g` vs `g`) or simply float them in the global namespace. The YAML parser evaluates both dot-notational groups and flattened keys identically.

---

## 2. Shared Settings Variables

Across most simulators (especially Damped, Driven, and Coupled), the simulation parameters and general program settings share identical namespacing:

```yaml
simulation:
  t_start: 0.0           # Starting time of the integration
  t_end: 100.0           # Ending time of the integration
  dt: 0.01               # Fixed simulation timestep precision
  output_every: 1        # Outputs metrics to CSV every N steps 

settings:
  integrator: "rk4"                        # Integrator to use (e.g., rk4, rk5, rk3, rk23, rkf45, velocity_verlet)
  data_file: "outputs/sim_data.csv"        # Path to dump CSV timeseries data
  output_png: "outputs/plot.png"           # Output path for the matplotlib plot
  python_script: "outputs/plot_script.py"  # Output path for the exported python graphing script
  run_plotter: true                        # If true, runs the generated python plot script
  show_plot: false                         # Opens an interactive Matplotlib GUI if true
  save_png: true                           # Exports the plot to `output_png` 
```

---

## 3. Simple Pendulum 
*A frictionless pendular system.* 

**Physical Properties:**
- `length` (or `L`) **[Required]**: Pendulum length.
- `gravity` (or `g`) **[Required]**: Gravitational acceleration.
- `theta0`: Initial angular displacement.
- `omega0`: Initial angular velocity.
- `restoring_force_model`: Set to `sine` for accurate physics, or `polynomial` for Taylor expansion approximations.
- `restoring_force_linear` / `restoring_force_cubic`: Coefficients mapped to polynomial restoring force.

**Special Settings:**
- `error_analysis` or `error_mode`: Set to `analytical`, `none`, or `hd_reference` to conduct integration error analytics against analytical solves.
- `additional_terms.*`: optional extra ODE terms shared with the damped/driven configs:
  - `inverse_cubic_enabled` / `inverse_cubic_strength`: adds a `+k / theta^3` Ermakov-Pinney term
  - `exponential_enabled` / `exponential_strength` / `exponential_scale` / `exponential_subtract_equilibrium`: adds a Toda-style exponential term `-a (exp(b*theta) - delta)`
  - `state_power_enabled` / `state_power_strength` / `state_power_exponent` / `state_power_mode`: adds a power-law restoring term `-c * theta^n`
  - `time_damping_enabled` / `time_damping_coefficient` / `time_damping_power` / `time_damping_shift`: adds a coordinate-dependent damping term `-(c / (t+shift)^p) * theta_dot`

---

## 4. Damped Pendulum
*Includes linear or cubic frictional forces.*

**Physical Properties:**
- `physical.L`: Pendulum length.
- `physical.g`: Gravitational acceleration.
- `physical.theta0`: Initial angular displacement.
- `physical.theta_dot0`: Initial angular velocity.
- `physical.damping_model`: Switch between `linear` or `polynomial` drag models.
- `physical.gamma` (or `damping` / `damping_linear`): Linear coefficient of damping.
- `physical.damping_cubic`: Cubic damping factor applied if polynomial model is active.
- `physical.additional_terms.*`: same extra-term block described above; this is how Lane-Emden style `2/t * theta_dot` and `theta^n` terms are represented.

---

## 5. Driven Systems & Duffing Oscillators
*For studying forcing, resonance cascades, chaos representations, and advanced van der Pol / Mathieu implementations.*

**Models**
- `physical.system_model` (or `mode`): Set to `pendulum` or `duffing`.

**Physical Properties**:
- `physical.mass`: Mass of the object (Defaults to 1.0 for Duffing modes).
- `physical.damping`: Primary environmental damping constant.
- `physical.linear_stiffness` / `physical.cubic_stiffness`: Stiffness coefficients specifically for the Duffing model.
- `physical.drive_force` (or `F`): External drive periodic force magnitude.
- `physical.omega_drive` (or `omega`): Angular frequency of the external force.
- `physical.drive_phase`: Phase offset added to the driving periodic component.
- `physical.parametric_amplitude` / `physical.parametric_frequency`: Time-variant amplitudes and frequencies used for Mathieu equation parameter excitations.

**Advanced Modules**:
- **Noise Engine (`noise.enabled`):**
  - `noise.force_stddev`: Standard deviation mapping for continuous random forcing.
  - `noise.correlation_time`: Color timeline for generated noise distributions.
- **Mass Events (`mass_event.enabled`):**
  - `mass_event.jump_time` / `mass_event.delta_mass`: Allows for simulating discrete dynamic structural failures mid-solve by shedding mass.
- **Frequency Sweeping (`sweep.enabled`):**
  - Automates sweeping across forced frequencies to build bifurcation or resonance charts!
  - Options: `sweep.omega_start`, `sweep.omega_end`, `sweep.points`, `sweep.direction` (`ascending`/`descending`), and `sweep.settle_time` (integration time skipped per-point to eliminate transients).
- `physical.additional_terms.*`: the same shared extra-term block is also available here, although analytical references remain limited to the linear pendulum branch.

**Analytical Model Names Added**
- `ermakov_pinney`: exact reference for the linear-plus-inverse-cubic Ermakov-Pinney case.
- `lane_emden`: exact spherical Lane-Emden reference for indices `n = 0`, `1`, or `5` when encoded via `additional_terms`.
- `toda` / `toda_liouville`: exact reference for the pure exponential Toda/Liouville branch (`exponential_subtract_equilibrium: false`).
- `toda_soliton` / `toda_soliton_continuum`: exact `sech^2` one-soliton profile for the continuum reduction of the Toda lattice.

---

## 6. Coupled Oscillators (2-DOF)
*2 Degree of Freedom systems, allowing energy transfer visualizations between interconnected harmonic/Duffing pendulums.* 

**Equations modelled:**
1) `q_1'' + 2*zeta_1*omega_1*q_1' + omega_1^2*q_1 + alpha_11*q_1^3 + alpha_12*q_1*q_2^2 = F*cos(omega_drive*t)`
2) `q_2'' + 2*zeta_2*omega_2*q_2' + omega_2^2*q_2 + alpha_22*q_2^3 + alpha_21*q_1^2*q_2 = 0`

**Physical Properties:**
- `omega_1` / `zeta_1`: Natural frequency and damping parameter for body 1.
- `alpha_11`: Primary nonlinear spring/cubic resistance for body 1.
- `alpha_12`: Transferred cross-coupling constraint inflicted by body 2 onto body 1.
- `omega_2` / `zeta_2` / `alpha_22` / `alpha_21`: Symmetrical configuration offsets for body 2.
- `F` / `omega_drive`: Forcing functions applied to driving body 1.

**Initial Conditions:**
- `q1_0`, `omega1_0` (or `q1_dot0`): Starting parameters for body 1.
- `q2_0`, `omega2_0` (or `q2_dot0`): Starting parameters for body 2.
