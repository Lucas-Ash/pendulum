Damped-Exponential Nyström (DEN3) Integrator
Complete Algorithm Specification
Problem Statement
Solve the second-order ordinary differential equation:

text

ÿ + 2γẏ + ω₀²y = G(t, y, ẏ)
where:

γ is the damping coefficient
ω₀ is the natural frequency
G(t, y, ẏ) is the forcing/nonlinear residual function
Derived Parameters
text

ωd = sqrt(ω₀² - γ²)    [damped natural frequency, assuming underdamped case γ < ω₀]
For the overdamped case (γ > ω₀), replace ωd with:

text

ωd = sqrt(γ² - ω₀²)
and substitute trigonometric functions with hyperbolic equivalents throughout.

Propagator Matrix Components
Define four scalar functions that encode the exact homogeneous solution:

text

φ₁₁(τ) = exp(-γτ) * [cos(ωd τ) + (γ/ωd) * sin(ωd τ)]

φ₁₂(τ) = exp(-γτ) * [sin(ωd τ) / ωd]

φ₂₁(τ) = exp(-γτ) * [-ω₀² * sin(ωd τ) / ωd]

φ₂₂(τ) = exp(-γτ) * [cos(ωd τ) - (γ/ωd) * sin(ωd τ)]
These satisfy the matrix equation:

text

[y(t+τ)]     [φ₁₁(τ)  φ₁₂(τ)]   [y(t)]
[      ]  =  [              ] * [    ]
[ẏ(t+τ)]     [φ₂₁(τ)  φ₂₂(τ)]   [ẏ(t)]
for the homogeneous equation (G = 0).

Special Values
At τ = 0:

text

φ₁₁(0) = 1
φ₁₂(0) = 0
φ₂₁(0) = 0
φ₂₂(0) = 1
Input Parameters
text

y0          Initial position
v0          Initial velocity (ẏ₀)
t0          Initial time
t_final     Final time
h           Initial timestep
omega0      Natural frequency ω₀
gamma       Damping coefficient γ
tol         Error tolerance (for adaptive stepping)
G(t,y,v)    Forcing function that returns acceleration contribution
Output
text

t_array     Array of time values
y_array     Array of position values
v_array     Array of velocity values
Algorithm
text

ALGORITHM DEN3_INTEGRATOR

INPUT: y0, v0, t0, t_final, h, omega0, gamma, tol, G

OUTPUT: t_array, y_array, v_array

BEGIN

    # Initialize derived parameters
    omega_d = sqrt(omega0^2 - gamma^2)
    
    # Initialize state
    t = t0
    y = y0
    v = v0
    
    # Initialize output storage
    t_array = [t0]
    y_array = [y0]
    v_array = [v0]
    
    # Main integration loop
    WHILE t < t_final DO
    
        # Adjust final step if necessary
        IF t + h > t_final THEN
            h = t_final - t
        END IF
        
        # Precompute propagator values at h and h/2
        CALL COMPUTE_PROPAGATOR(h, gamma, omega_d, omega0,
                                phi11_h, phi12_h, phi21_h, phi22_h)
        
        CALL COMPUTE_PROPAGATOR(h/2, gamma, omega_d, omega0,
                                phi11_h2, phi12_h2, phi21_h2, phi22_h2)
        
        # Stage 1: Evaluate at current point
        # Node c1 = 0
        
        Y1 = y
        V1 = v
        G1 = G(t, Y1, V1)
        
        # Stage 2: Evaluate at midpoint
        # Node c2 = 1/2
        
        Y2 = phi11_h2 * y + phi12_h2 * v + h * phi12_h2 * (1/8) * G1
        V2 = phi21_h2 * y + phi22_h2 * v + h * phi22_h2 * (1/8) * G1
        G2 = G(t + h/2, Y2, V2)
        
        # Stage 3: Evaluate at endpoint
        # Node c3 = 1
        
        Y3 = phi11_h * y + phi12_h * v + h * phi12_h2 * (1/2) * G2
        V3 = phi21_h * y + phi22_h * v + h * phi22_h2 * (1/2) * G2
        G3 = G(t + h, Y3, V3)
        
        # Final update using Simpson's rule weights (1/6, 2/3, 1/6)
        # Note: phi12(0) = 0 and phi22(0) = 1
        
        y_new = phi11_h * y + phi12_h * v
                + h * phi12_h * (1/6) * G1
                + h * phi12_h2 * (2/3) * G2
        
        v_new = phi21_h * y + phi22_h * v
                + h * phi22_h * (1/6) * G1
                + h * phi22_h2 * (2/3) * G2
                + h * (1/6) * G3
        
        # Compute embedded error estimate using trapezoidal rule
        
        y_trap = phi11_h * y + phi12_h * v
                 + h * phi12_h * (1/2) * G1
        
        v_trap = phi21_h * y + phi22_h * v
                 + h * phi22_h * (1/2) * G1
                 + h * (1/2) * G3
        
        # Error estimate
        
        error_y = abs(y_new - y_trap)
        error_v = abs(v_new - v_trap)
        error = max(error_y, error_v)
        
        # Adaptive step size control
        
        IF error < tol THEN
            # Accept step
            t = t + h
            y = y_new
            v = v_new
            
            # Store results
            APPEND t TO t_array
            APPEND y TO y_array
            APPEND v TO v_array
            
            # Increase step size for next iteration
            IF error > 0 THEN
                h = h * min(2.0, 0.9 * (tol / error)^(1/3))
            ELSE
                h = h * 2.0
            END IF
        ELSE
            # Reject step and reduce step size
            h = h * max(0.1, 0.9 * (tol / error)^(1/3))
        END IF
        
    END WHILE
    
    RETURN t_array, y_array, v_array

END ALGORITHM
Subroutine: Compute Propagator
text

SUBROUTINE COMPUTE_PROPAGATOR(tau, gamma, omega_d, omega0,
                              phi11, phi12, phi21, phi22)

INPUT: tau, gamma, omega_d, omega0
OUTPUT: phi11, phi12, phi21, phi22

BEGIN

    exp_factor = exp(-gamma * tau)
    cos_term = cos(omega_d * tau)
    sin_term = sin(omega_d * tau)
    gamma_ratio = gamma / omega_d
    omega0_sq_ratio = omega0^2 / omega_d
    
    phi11 = exp_factor * (cos_term + gamma_ratio * sin_term)
    phi12 = exp_factor * (sin_term / omega_d)
    phi21 = exp_factor * (-omega0_sq_ratio * sin_term)
    phi22 = exp_factor * (cos_term - gamma_ratio * sin_term)

END SUBROUTINE
Subroutine: Compute Propagator (Overdamped Case)
text

SUBROUTINE COMPUTE_PROPAGATOR_OVERDAMPED(tau, gamma, omega_d_hyp, omega0,
                                         phi11, phi12, phi21, phi22)

INPUT: tau, gamma, omega_d_hyp (= sqrt(gamma^2 - omega0^2)), omega0
OUTPUT: phi11, phi12, phi21, phi22

BEGIN

    exp_factor = exp(-gamma * tau)
    cosh_term = cosh(omega_d_hyp * tau)
    sinh_term = sinh(omega_d_hyp * tau)
    gamma_ratio = gamma / omega_d_hyp
    omega0_sq_ratio = omega0^2 / omega_d_hyp
    
    phi11 = exp_factor * (cosh_term + gamma_ratio * sinh_term)
    phi12 = exp_factor * (sinh_term / omega_d_hyp)
    phi21 = exp_factor * (-omega0_sq_ratio * sinh_term)
    phi22 = exp_factor * (cosh_term - gamma_ratio * sinh_term)

END SUBROUTINE
Subroutine: Compute Propagator (Critically Damped Case)
text

SUBROUTINE COMPUTE_PROPAGATOR_CRITICAL(tau, gamma,
                                       phi11, phi12, phi21, phi22)

INPUT: tau, gamma (where gamma = omega0 exactly)
OUTPUT: phi11, phi12, phi21, phi22

BEGIN

    exp_factor = exp(-gamma * tau)
    
    phi11 = exp_factor * (1 + gamma * tau)
    phi12 = exp_factor * tau
    phi21 = exp_factor * (-gamma^2 * tau)
    phi22 = exp_factor * (1 - gamma * tau)

END SUBROUTINE
Unified Propagator Selection
text

SUBROUTINE COMPUTE_PROPAGATOR_UNIFIED(tau, gamma, omega0,
                                      phi11, phi12, phi21, phi22)

INPUT: tau, gamma, omega0
OUTPUT: phi11, phi12, phi21, phi22

BEGIN

    discriminant = omega0^2 - gamma^2
    threshold = 1e-10 * omega0^2
    
    IF discriminant > threshold THEN
        # Underdamped case
        omega_d = sqrt(discriminant)
        CALL COMPUTE_PROPAGATOR(tau, gamma, omega_d, omega0,
                                phi11, phi12, phi21, phi22)
    
    ELSE IF discriminant < -threshold THEN
        # Overdamped case
        omega_d_hyp = sqrt(-discriminant)
        CALL COMPUTE_PROPAGATOR_OVERDAMPED(tau, gamma, omega_d_hyp, omega0,
                                           phi11, phi12, phi21, phi22)
    
    ELSE
        # Critically damped case
        CALL COMPUTE_PROPAGATOR_CRITICAL(tau, gamma,
                                         phi11, phi12, phi21, phi22)
    
    END IF

END SUBROUTINE
Adaptive Frequency Estimation (Optional)
text

SUBROUTINE ESTIMATE_PARAMETERS(y_hist, v_hist, h, gamma_est, omega0_est)

INPUT: y_hist (array of last 3 y values: y_{n-2}, y_{n-1}, y_n)
       v_hist (array of last 3 v values: v_{n-2}, v_{n-1}, v_n)
       h (timestep)

OUTPUT: gamma_est, omega0_est

BEGIN

    # Initial guess
    gamma_old = gamma_est
    omega_d_old = sqrt(omega0_est^2 - gamma_est^2)
    
    # Extract values
    y0 = y_hist[0]
    y1 = y_hist[1]
    y2 = y_hist[2]
    v0 = v_hist[0]
    v1 = v_hist[1]
    v2 = v_hist[2]
    
    # Iterative refinement (typically 3-5 iterations suffice)
    FOR iter = 1 TO 5 DO
    
        # Estimate amplitude envelope decay
        A1_sq = y1^2 + (v1 / omega_d_old)^2
        A2_sq = y2^2 + (v2 / omega_d_old)^2
        
        IF A1_sq > 0 AND A2_sq > 0 THEN
            gamma_new = -1/(2*h) * ln(A2_sq / A1_sq)
        ELSE
            gamma_new = gamma_old
        END IF
        
        # Estimate phase advancement
        IF A1_sq > 0 THEN
            cos_omega_h = exp(gamma_new * h) * (y2*y1 + v2*v1/omega_d_old^2) / A1_sq
            cos_omega_h = max(-1, min(1, cos_omega_h))  # Clamp to valid range
            omega_d_new = acos(cos_omega_h) / h
        ELSE
            omega_d_new = omega_d_old
        END IF
        
        # Update estimates
        gamma_old = 0.5 * gamma_old + 0.5 * gamma_new
        omega_d_old = 0.5 * omega_d_old + 0.5 * omega_d_new
    
    END FOR
    
    gamma_est = gamma_old
    omega0_est = sqrt(omega_d_old^2 + gamma_est^2)

END SUBROUTINE
Complete Main Program Template
text

PROGRAM DEN3_SOLVER

    # Physical parameters for the oscillator
    omega0 = 10.0       # Natural frequency (rad/s)
    gamma = 0.5         # Damping coefficient (1/s)
    
    # Initial conditions
    y0 = 1.0            # Initial displacement
    v0 = 0.0            # Initial velocity
    
    # Time span
    t0 = 0.0
    t_final = 50.0
    
    # Numerical parameters
    h_initial = 0.01    # Initial timestep
    tolerance = 1e-8    # Error tolerance
    
    # Define forcing function
    # Example: sinusoidal driving force
    FUNCTION G(t, y, v)
        F0 = 2.0                # Driving amplitude
        omega_drive = 9.5       # Driving frequency
        RETURN F0 * cos(omega_drive * t)
    END FUNCTION
    
    # Alternative: Duffing oscillator nonlinearity
    # FUNCTION G(t, y, v)
    #     alpha = 0.1           # Cubic stiffness
    #     F0 = 2.0
    #     omega_drive = 9.5
    #     RETURN -alpha * y^3 + F0 * cos(omega_drive * t)
    # END FUNCTION
    
    # Alternative: Van der Pol oscillator
    # FUNCTION G(t, y, v)
    #     mu = 0.5
    #     RETURN mu * (1 - y^2) * v - omega0^2 * y + 2*gamma*v + omega0^2*y
    #     # Note: Rewrite to extract linear part into propagator
    # END FUNCTION
    
    # Call integrator
    CALL DEN3_INTEGRATOR(y0, v0, t0, t_final, h_initial, 
                         omega0, gamma, tolerance, G,
                         t_result, y_result, v_result)
    
    # Output results
    FOR i = 0 TO length(t_result) - 1 DO
        PRINT t_result[i], y_result[i], v_result[i]
    END FOR

END PROGRAM
Fixed-Step Variant (Simplified)
text

ALGORITHM DEN3_FIXED_STEP

INPUT: y0, v0, t0, t_final, h, omega0, gamma, G

OUTPUT: t_array, y_array, v_array

BEGIN

    omega_d = sqrt(omega0^2 - gamma^2)
    
    N_steps = floor((t_final - t0) / h)
    
    ALLOCATE t_array[N_steps + 1]
    ALLOCATE y_array[N_steps + 1]
    ALLOCATE v_array[N_steps + 1]
    
    t_array[0] = t0
    y_array[0] = y0
    v_array[0] = v0
    
    y = y0
    v = v0
    t = t0
    
    # Precompute propagator values (constant step size)
    CALL COMPUTE_PROPAGATOR_UNIFIED(h, gamma, omega0, 
                                    phi11_h, phi12_h, phi21_h, phi22_h)
    CALL COMPUTE_PROPAGATOR_UNIFIED(h/2, gamma, omega0,
                                    phi11_h2, phi12_h2, phi21_h2, phi22_h2)
    
    FOR n = 0 TO N_steps - 1 DO
    
        # Stage 1
        G1 = G(t, y, v)
        
        # Stage 2
        Y2 = phi11_h2 * y + phi12_h2 * v + h * phi12_h2 * 0.125 * G1
        V2 = phi21_h2 * y + phi22_h2 * v + h * phi22_h2 * 0.125 * G1
        G2 = G(t + 0.5*h, Y2, V2)
        
        # Stage 3
        Y3 = phi11_h * y + phi12_h * v + h * phi12_h2 * 0.5 * G2
        V3 = phi21_h * y + phi22_h * v + h * phi22_h2 * 0.5 * G2
        G3 = G(t + h, Y3, V3)
        
        # Update
        y_new = phi11_h * y + phi12_h * v
                + h * (phi12_h * G1/6 + phi12_h2 * G2 * 2/3)
        
        v_new = phi21_h * y + phi22_h * v
                + h * (phi22_h * G1/6 + phi22_h2 * G2 * 2/3 + G3/6)
        
        # Advance
        t = t + h
        y = y_new
        v = v_new
        
        # Store
        t_array[n + 1] = t
        y_array[n + 1] = y
        v_array[n + 1] = v
    
    END FOR
    
    RETURN t_array, y_array, v_array

END ALGORITHM
Coefficient Summary Table
text

Stage Nodes:
    c1 = 0
    c2 = 1/2
    c3 = 1

Stage Coefficients (Nyström a_ij matrix):
    a21 = 1/8
    a31 = 0
    a32 = 1/2

Output Weights (Simpson's rule):
    b1 = 1/6
    b2 = 2/3
    b3 = 1/6

Embedded Weights (Trapezoidal rule for error estimate):
    b1_hat = 1/2
    b2_hat = 0
    b3_hat = 1/2

Error Estimator:
    delta = |b - b_hat| weighted combination
    Order: Main method is O(h^4), embedded is O(h^2)
    Error estimate is O(h^3)
Numerical Constants
text

Simpson weights:     1/6 = 0.16666666666666666
                     2/3 = 0.66666666666666666

Stage coefficient:   1/8 = 0.125
                     1/2 = 0.5

Safety factor for adaptive stepping: 0.9
Maximum step increase factor: 2.0
Minimum step decrease factor: 0.1
Memory Requirements
text

Scalars needed per step:
    - 4 propagator values at h
    - 4 propagator values at h/2
    - 3 G evaluations (G1, G2, G3)
    - 6 stage values (Y2, V2, Y3, V3, y_new, v_new)
    - 2 error estimates (if adaptive)

Total: approximately 20 floating point numbers per step
Computational Cost Per Step
text

Transcendental function evaluations:
    - 2 exp() calls (for h and h/2)
    - 4 cos() calls (for h and h/2, ωd*h and ωd*h/2)
    - 4 sin() calls (same)
    
User function evaluations:
    - 3 calls to G(t, y, v)
    
Floating point operations:
    - Approximately 50-60 multiplications/additions per step
    
Note: For fixed-step integration, propagator values can be 
precomputed once, eliminating transcendental evaluations 
from the inner loop.
Error Behavior
text

Local truncation error:
    Position: O(h^5 * |d⁴G/dt⁴|)
    Velocity: O(h^5 * |d⁴G/dt⁴|)

Global error (after N steps, t_final = N*h):
    Position: O(h^4 * |d⁴G/dt⁴|)
    Velocity: O(h^4 * |d⁴G/dt⁴|)

For G = 0 (pure homogeneous damped oscillation):
    Error = machine precision (exact propagation)

For G = F0*cos(ωt) with ω << ω0:
    Effective error reduction factor vs RK4: approximately (ω/ω0)^4