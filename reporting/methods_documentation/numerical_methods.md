# Numerical Integration Methods

This document describes the numerical integration methods implemented in
`modules/integrators/rk_integrators.h`. The methods fall into four families:

1. **Explicit Runge--Kutta methods** (RK3, RK4, RK5, RK23, RKF45)
2. **Symplectic integrators** (Semi-Implicit Euler, Leapfrog, Ruth4)
3. **Second-order ODE specialists** (Velocity Verlet, Runge--Kutta--Nyström, Numerov)
4. **Exponential integrators** (DEN3)

Throughout this document the ODE system is either the first-order form

\[
\frac{d}{dt}\mathbf{y} = \mathbf{f}(t,\,\mathbf{y}),
\qquad
\mathbf{y} = \begin{pmatrix}\theta \\ \omega\end{pmatrix},
\]

or the direct second-order form

\[
\ddot\theta = a(t,\,\theta).
\]

The time step is \(h\) (code variable `dt`), and subscript \(n\) denotes the
value at time \(t_n\).

---

## 1  Classical Runge--Kutta (4th Order) — `rk4_step`

### Description

The classical four-stage, fourth-order explicit Runge--Kutta method. It is the
most widely used single-step ODE integrator and provides a good balance between
accuracy and computational cost.

### Formulation

\[
\begin{aligned}
\mathbf{k}_1 &= \mathbf{f}(t_n,\;\mathbf{y}_n) \\
\mathbf{k}_2 &= \mathbf{f}\!\left(t_n + \tfrac{h}{2},\;\mathbf{y}_n + \tfrac{h}{2}\,\mathbf{k}_1\right) \\
\mathbf{k}_3 &= \mathbf{f}\!\left(t_n + \tfrac{h}{2},\;\mathbf{y}_n + \tfrac{h}{2}\,\mathbf{k}_2\right) \\
\mathbf{k}_4 &= \mathbf{f}\!\left(t_n + h,\;\mathbf{y}_n + h\,\mathbf{k}_3\right) \\[6pt]
\mathbf{y}_{n+1} &= \mathbf{y}_n + \frac{h}{6}\bigl(\mathbf{k}_1 + 2\,\mathbf{k}_2 + 2\,\mathbf{k}_3 + \mathbf{k}_4\bigr)
\end{aligned}
\]

### Butcher Tableau

```
  0   |
  1/2 | 1/2
  1/2 | 0    1/2
  1   | 0    0    1
  ----|------------------
      | 1/6  1/3  1/3  1/6
```

### Properties

| Property          | Value                              |
|-------------------|------------------------------------|
| Order             | 4                                  |
| Stages            | 4                                  |
| Implicit          | No                                 |
| Symplectic        | No                                 |
| Local error       | \(\mathcal{O}(h^5)\)              |

---

## 2  Standard Runge--Kutta (3rd Order) — `rk3_step`

### Description

A three-stage, third-order explicit Runge--Kutta method. Lower cost per step
than RK4 but also lower accuracy; useful for convergence-order studies.

### Formulation

\[
\begin{aligned}
\mathbf{k}_1 &= \mathbf{f}(t_n,\;\mathbf{y}_n) \\
\mathbf{k}_2 &= \mathbf{f}\!\left(t_n + \tfrac{h}{2},\;\mathbf{y}_n + \tfrac{h}{2}\,\mathbf{k}_1\right) \\
\mathbf{k}_3 &= \mathbf{f}\!\left(t_n + h,\;\mathbf{y}_n - h\,\mathbf{k}_1 + 2h\,\mathbf{k}_2\right) \\[6pt]
\mathbf{y}_{n+1} &= \mathbf{y}_n + \frac{h}{6}\bigl(\mathbf{k}_1 + 4\,\mathbf{k}_2 + \mathbf{k}_3\bigr)
\end{aligned}
\]

### Butcher Tableau

```
  0   |
  1/2 | 1/2
  1   | -1   2
  ----|------------
      | 1/6  2/3  1/6
```

### Properties

| Property          | Value                              |
|-------------------|------------------------------------|
| Order             | 3                                  |
| Stages            | 3                                  |
| Implicit          | No                                 |
| Symplectic        | No                                 |
| Local error       | \(\mathcal{O}(h^4)\)              |

---

## 3  Butcher's Runge--Kutta (5th Order) — `rk5_step`

### Description

A six-stage, fifth-order explicit Runge--Kutta method due to Butcher (1964).
Higher accuracy per step than RK4 at the cost of two additional derivative
evaluations.

### Formulation

\[
\begin{aligned}
\mathbf{k}_1 &= \mathbf{f}(t_n,\;\mathbf{y}_n) \\
\mathbf{k}_2 &= \mathbf{f}\!\left(t_n + \tfrac{h}{4},\;\mathbf{y}_n + \tfrac{h}{4}\,\mathbf{k}_1\right) \\
\mathbf{k}_3 &= \mathbf{f}\!\left(t_n + \tfrac{h}{4},\;\mathbf{y}_n + \tfrac{h}{8}\,\mathbf{k}_1 + \tfrac{h}{8}\,\mathbf{k}_2\right) \\
\mathbf{k}_4 &= \mathbf{f}\!\left(t_n + \tfrac{h}{2},\;\mathbf{y}_n - \tfrac{h}{2}\,\mathbf{k}_2 + h\,\mathbf{k}_3\right) \\
\mathbf{k}_5 &= \mathbf{f}\!\left(t_n + \tfrac{3h}{4},\;\mathbf{y}_n + \tfrac{3h}{16}\,\mathbf{k}_1 + \tfrac{9h}{16}\,\mathbf{k}_4\right) \\
\mathbf{k}_6 &= \mathbf{f}\!\left(t_n + h,\;\mathbf{y}_n - \tfrac{3h}{7}\,\mathbf{k}_1 + \tfrac{2h}{7}\,\mathbf{k}_2 + \tfrac{12h}{7}\,\mathbf{k}_3 - \tfrac{12h}{7}\,\mathbf{k}_4 + \tfrac{8h}{7}\,\mathbf{k}_5\right) \\[6pt]
\mathbf{y}_{n+1} &= \mathbf{y}_n + \frac{h}{90}\bigl(7\,\mathbf{k}_1 + 32\,\mathbf{k}_3 + 12\,\mathbf{k}_4 + 32\,\mathbf{k}_5 + 7\,\mathbf{k}_6\bigr)
\end{aligned}
\]

### Butcher Tableau

```
  0    |
  1/4  | 1/4
  1/4  | 1/8     1/8
  1/2  | 0      -1/2     1
  3/4  | 3/16    0       0      9/16
  1    | -3/7    2/7     12/7  -12/7    8/7
  -----|------------------------------------------
       | 7/90    0      32/90  12/90   32/90   7/90
```

### Properties

| Property          | Value                              |
|-------------------|------------------------------------|
| Order             | 5                                  |
| Stages            | 6                                  |
| Implicit          | No                                 |
| Symplectic        | No                                 |
| Local error       | \(\mathcal{O}(h^6)\)              |

---

## 4  Bogacki--Shampine (RK23) — `rk23_step`

### Description

The Bogacki--Shampine method is a three-stage embedded Runge--Kutta pair of
orders 2 and 3. In the implementation here the third-order solution is used to
advance the state. The embedded pair structure makes it well suited for adaptive
step-size control (only the fixed-step advance is implemented).

### Formulation

\[
\begin{aligned}
\mathbf{k}_1 &= \mathbf{f}(t_n,\;\mathbf{y}_n) \\
\mathbf{k}_2 &= \mathbf{f}\!\left(t_n + \tfrac{h}{2},\;\mathbf{y}_n + \tfrac{h}{2}\,\mathbf{k}_1\right) \\
\mathbf{k}_3 &= \mathbf{f}\!\left(t_n + \tfrac{3h}{4},\;\mathbf{y}_n + \tfrac{3h}{4}\,\mathbf{k}_2\right) \\[6pt]
\mathbf{y}_{n+1} &= \mathbf{y}_n + h\left(\tfrac{2}{9}\,\mathbf{k}_1 + \tfrac{1}{3}\,\mathbf{k}_2 + \tfrac{4}{9}\,\mathbf{k}_3\right)
\end{aligned}
\]

### Butcher Tableau

```
  0    |
  1/2  | 1/2
  3/4  | 0    3/4
  -----|------------
       | 2/9  1/3  4/9
```

### Properties

| Property          | Value                              |
|-------------------|------------------------------------|
| Order (advance)   | 3                                  |
| Order (embedded)  | 2                                  |
| Stages            | 3                                  |
| Implicit          | No                                 |
| Symplectic        | No                                 |
| Local error       | \(\mathcal{O}(h^4)\)              |

---

## 5  Runge--Kutta--Fehlberg (RKF45) — `rkf45_step`

### Description

The Runge--Kutta--Fehlberg method is a six-stage embedded pair of orders 4
and 5. The fifth-order solution is used to advance the state, making it an
efficient high-accuracy integrator. Originally designed for adaptive step-size
control; the implementation here uses fixed steps.

### Formulation

\[
\begin{aligned}
\mathbf{k}_1 &= \mathbf{f}(t_n,\;\mathbf{y}_n) \\
\mathbf{k}_2 &= \mathbf{f}\!\left(t_n + \tfrac{h}{4},\;\mathbf{y}_n + \tfrac{h}{4}\,\mathbf{k}_1\right) \\
\mathbf{k}_3 &= \mathbf{f}\!\left(t_n + \tfrac{3h}{8},\;\mathbf{y}_n + \tfrac{3h}{32}\,\mathbf{k}_1 + \tfrac{9h}{32}\,\mathbf{k}_2\right) \\
\mathbf{k}_4 &= \mathbf{f}\!\left(t_n + \tfrac{12h}{13},\;\mathbf{y}_n + \tfrac{1932h}{2197}\,\mathbf{k}_1 - \tfrac{7200h}{2197}\,\mathbf{k}_2 + \tfrac{7296h}{2197}\,\mathbf{k}_3\right) \\
\mathbf{k}_5 &= \mathbf{f}\!\left(t_n + h,\;\mathbf{y}_n + \tfrac{439h}{216}\,\mathbf{k}_1 - 8h\,\mathbf{k}_2 + \tfrac{3680h}{513}\,\mathbf{k}_3 - \tfrac{845h}{4104}\,\mathbf{k}_4\right) \\
\mathbf{k}_6 &= \mathbf{f}\!\left(t_n + \tfrac{h}{2},\;\mathbf{y}_n - \tfrac{8h}{27}\,\mathbf{k}_1 + 2h\,\mathbf{k}_2 - \tfrac{3544h}{2565}\,\mathbf{k}_3 + \tfrac{1859h}{4104}\,\mathbf{k}_4 - \tfrac{11h}{40}\,\mathbf{k}_5\right) \\[6pt]
\mathbf{y}_{n+1} &= \mathbf{y}_n + h\left(\tfrac{16}{135}\,\mathbf{k}_1 + \tfrac{6656}{12825}\,\mathbf{k}_3 + \tfrac{28561}{56430}\,\mathbf{k}_4 - \tfrac{9}{50}\,\mathbf{k}_5 + \tfrac{2}{55}\,\mathbf{k}_6\right)
\end{aligned}
\]

### Butcher Tableau

```
  0      |
  1/4    | 1/4
  3/8    | 3/32        9/32
  12/13  | 1932/2197  -7200/2197   7296/2197
  1      | 439/216    -8           3680/513   -845/4104
  1/2    | -8/27       2          -3544/2565   1859/4104  -11/40
  -------|-------------------------------------------------------
  5th    | 16/135      0           6656/12825  28561/56430 -9/50  2/55
```

### Properties

| Property          | Value                              |
|-------------------|------------------------------------|
| Order (advance)   | 5                                  |
| Order (embedded)  | 4                                  |
| Stages            | 6                                  |
| Implicit          | No                                 |
| Symplectic        | No                                 |
| Local error       | \(\mathcal{O}(h^6)\)              |

---

## 6  Semi-Implicit Euler (Symplectic 1st Order) — `semi_implicit_euler_step`

### Description

The semi-implicit (symplectic) Euler method updates velocity first using the
current-state derivative, then advances position using the *new* velocity. This
asymmetric coupling makes the method symplectic, preserving the phase-space
structure of Hamiltonian systems and providing long-term energy stability that
the standard explicit Euler lacks.

### Formulation

\[
\begin{aligned}
\omega_{n+1} &= \omega_n + h\,f_\omega(t_n,\;\theta_n,\;\omega_n) \\
\theta_{n+1} &= \theta_n + h\,\omega_{n+1}
\end{aligned}
\]

where \(f_\omega\) is the angular-acceleration component of the derivative.

### Properties

| Property          | Value                              |
|-------------------|------------------------------------|
| Order             | 1                                  |
| Stages            | 1                                  |
| Implicit          | Semi-implicit                      |
| Symplectic        | Yes                                |
| Local error       | \(\mathcal{O}(h^2)\)              |

---

## 7  Leapfrog / Störmer--Verlet (Symplectic 2nd Order) — `leapfrog_step`

### Description

The leapfrog method (also known as the kick-drift-kick form of Störmer--Verlet)
is a second-order symplectic integrator. It advances position by a half-step,
evaluates the force at the midpoint, updates velocity over the full step, then
completes the position half-step. It is time-reversible and preserves
phase-space volume.

### Formulation

\[
\begin{aligned}
\theta_{n+1/2} &= \theta_n + \tfrac{h}{2}\,\omega_n \\
\omega_{n+1}   &= \omega_n + h\,f_\omega\!\left(t_n + \tfrac{h}{2},\;\theta_{n+1/2}\right) \\
\theta_{n+1}   &= \theta_{n+1/2} + \tfrac{h}{2}\,\omega_{n+1}
\end{aligned}
\]

### Properties

| Property          | Value                              |
|-------------------|------------------------------------|
| Order             | 2                                  |
| Stages            | 1 force evaluation                 |
| Implicit          | No (explicit for separable \(H\))  |
| Symplectic        | Yes                                |
| Time-reversible   | Yes                                |
| Local error       | \(\mathcal{O}(h^3)\)              |

---

## 8  Ruth's 4th-Order Symplectic Integrator — `ruth4_step`

### Description

A fourth-order symplectic integrator based on the Yoshida/Forest--Ruth
decomposition. It splits the Hamiltonian flow into alternating position kicks
and momentum kicks with specially chosen coefficients that cancel lower-order
error terms while preserving the symplectic structure exactly.

### Coefficients

Define \(w = 2^{1/3}\) and \(\alpha = \dfrac{1}{2 - w}\). The position
coefficients \(c_i\) and momentum coefficients \(d_i\) are:

\[
\begin{aligned}
c_1 &= \frac{\alpha}{2}, &\quad d_1 &= \alpha, \\
c_2 &= \frac{(1-w)\,\alpha}{2}, &\quad d_2 &= -w\,\alpha, \\
c_3 &= c_2, &\quad d_3 &= \alpha, \\
c_4 &= c_1.
\end{aligned}
\]

### Formulation

The integrator applies four position drifts interleaved with three momentum
kicks:

\[
\begin{aligned}
\theta &\leftarrow \theta + c_1\,h\,\omega \\
\omega &\leftarrow \omega + d_1\,h\,f_\omega(\theta) \\
\theta &\leftarrow \theta + c_2\,h\,\omega \\
\omega &\leftarrow \omega + d_2\,h\,f_\omega(\theta) \\
\theta &\leftarrow \theta + c_3\,h\,\omega \\
\omega &\leftarrow \omega + d_3\,h\,f_\omega(\theta) \\
\theta &\leftarrow \theta + c_4\,h\,\omega
\end{aligned}
\]

### Properties

| Property          | Value                                  |
|-------------------|----------------------------------------|
| Order             | 4                                      |
| Stages            | 3 force evaluations, 4 drift steps     |
| Implicit          | No (for separable Hamiltonians)        |
| Symplectic        | Yes                                    |
| Time-reversible   | Yes                                    |
| Local error       | \(\mathcal{O}(h^5)\)                  |

---

## 9  Velocity Verlet — `velocity_verlet_step`

### Description

The Velocity Verlet algorithm is a second-order method designed for the
direct second-order ODE \(\ddot\theta = a(t,\,\theta)\). It is algebraically
equivalent to the Störmer--Verlet method but is formulated to give
synchronous position and velocity output at each step. For conservative
(velocity-independent) forces it is symplectic.

### Formulation

\[
\begin{aligned}
a_n         &= a(t_n,\;\theta_n) \\
\theta_{n+1} &= \theta_n + h\,\omega_n + \tfrac{h^2}{2}\,a_n \\
a_{n+1}     &= a(t_{n+1},\;\theta_{n+1}) \\
\omega_{n+1} &= \omega_n + \tfrac{h}{2}\bigl(a_n + a_{n+1}\bigr)
\end{aligned}
\]

### Properties

| Property          | Value                              |
|-------------------|------------------------------------|
| Order             | 2                                  |
| Force evaluations | 2 per step (1 reusable)            |
| Implicit          | No                                 |
| Symplectic        | Yes (conservative forces)          |
| Local error       | \(\mathcal{O}(h^3)\)              |

---

## 10  Classical Runge--Kutta--Nyström (4th Order) — `runge_kutta_nystrom_step`

### Description

The Runge--Kutta--Nyström (RKN) method is a fourth-order integrator
specialised for second-order ODEs of the form \(\ddot\theta = a(t,\,\theta)\).
By exploiting the second-order structure directly it avoids reducing the
problem to a first-order system, and produces both \(\theta\) and \(\omega\)
to fourth-order accuracy.

### Formulation

\[
\begin{aligned}
k_{1,\theta} &= \omega_n, &\quad k_{1,\omega} &= a(t_n,\;\theta_n) \\[4pt]
k_{2,\theta} &= \omega_n + \tfrac{h}{2}\,k_{1,\omega}, &\quad
k_{2,\omega} &= a\!\left(t_n + \tfrac{h}{2},\;\theta_n + \tfrac{h}{2}\,k_{1,\theta}\right) \\[4pt]
k_{3,\theta} &= \omega_n + \tfrac{h}{2}\,k_{2,\omega}, &\quad
k_{3,\omega} &= a\!\left(t_n + \tfrac{h}{2},\;\theta_n + \tfrac{h}{2}\,k_{2,\theta}\right) \\[4pt]
k_{4,\theta} &= \omega_n + h\,k_{3,\omega}, &\quad
k_{4,\omega} &= a\!\left(t_n + h,\;\theta_n + h\,k_{3,\theta}\right) \\[6pt]
\theta_{n+1} &= \theta_n + \frac{h}{6}\bigl(k_{1,\theta} + 2\,k_{2,\theta} + 2\,k_{3,\theta} + k_{4,\theta}\bigr) \\[4pt]
\omega_{n+1} &= \omega_n + \frac{h}{6}\bigl(k_{1,\omega} + 2\,k_{2,\omega} + 2\,k_{3,\omega} + k_{4,\omega}\bigr)
\end{aligned}
\]

### Properties

| Property          | Value                              |
|-------------------|------------------------------------|
| Order             | 4                                  |
| Force evaluations | 4                                  |
| Implicit          | No                                 |
| Local error       | \(\mathcal{O}(h^5)\)              |

---

## 11  Numerov Method (Implicit, 4th Order) — `numerov_next_theta`

### Description

Numerov's method is a high-accuracy multistep scheme for second-order ODEs
with no first-derivative term, \(\ddot\theta = a(t,\,\theta)\). It achieves
fourth-order accuracy using only three consecutive position values and their
accelerations. The method is implicit because \(a_{n+1}\) depends on the
unknown \(\theta_{n+1}\); the implementation resolves this with Newton
iteration using a numerically estimated Jacobian.

### Recurrence

\[
\theta_{n+1} - 2\,\theta_n + \theta_{n-1}
= \frac{h^2}{12}\bigl(a_{n-1} + 10\,a_n + a_{n+1}\bigr)
\]

where \(a_k = a(t_k,\,\theta_k)\).

### Newton Iteration

Since \(a_{n+1} = a(t_{n+1},\,\theta_{n+1})\) depends on the unknown, the
recurrence is rearranged into the root-finding problem

\[
F(\theta_{n+1}) = \theta_{n+1}
  - \underbrace{\left[2\,\theta_n - \theta_{n-1} + \frac{h^2}{12}\bigl(a_{n-1} + 10\,a_n\bigr)\right]}_{\text{RHS}}
  - \frac{h^2}{12}\,a(t_{n+1},\,\theta_{n+1}) = 0.
\]

The Jacobian is approximated via central finite differences:

\[
\frac{\partial a}{\partial\theta}\bigg|_{\theta_{n+1}}
\approx \frac{a(t_{n+1},\,\theta_{n+1}+\varepsilon) - a(t_{n+1},\,\theta_{n+1}-\varepsilon)}{2\varepsilon}
\]

and the Newton update is

\[
\theta_{n+1} \leftarrow \theta_{n+1} - \frac{F(\theta_{n+1})}{1 - \tfrac{h^2}{12}\,\partial a/\partial\theta}.
\]

### Bootstrap and Velocity Recovery

Numerov requires two prior position values. The implementation bootstraps
\(\theta_1\) and \(\theta_2\) using the 4th-order Runge--Kutta--Nyström method.
All arithmetic during the Numerov phase is performed in `long double` for
enhanced precision.

Velocities are recovered *a posteriori* from the position trajectory using the
corrected central-difference formula:

\[
\omega_i = \frac{\theta_{i+1} - \theta_{i-1}}{2h}
           - \frac{h}{12}\bigl(a_{i+1} - a_{i-1}\bigr),
\qquad 1 \le i \le N-1,
\]

and the endpoint velocity uses a five-point backward-difference stencil:

\[
\omega_N = \frac{25\,\theta_N - 48\,\theta_{N-1} + 36\,\theta_{N-2} - 16\,\theta_{N-3} + 3\,\theta_{N-4}}{12h}.
\]

### Properties

| Property          | Value                                          |
|-------------------|------------------------------------------------|
| Order             | 4 (global)                                     |
| Type              | Multistep (requires two prior positions)       |
| Implicit          | Yes (Newton iteration per step)                |
| Precision         | `long double` arithmetic                       |
| Local error       | \(\mathcal{O}(h^6)\)                          |

---

## 12  Discrete Exponential Nyström (DEN3) — `den3_step`

### Description

The DEN3 method is a third-order exponential integrator designed for the
damped nonlinear oscillator equation

\[
\ddot\theta + 2\gamma\,\dot\theta + \omega_0^2\,\theta = g(t,\,\theta,\,\dot\theta),
\]

where \(\gamma\) is the damping coefficient, \(\omega_0\) is the natural
frequency, and \(g\) is the nonlinear residual. The key idea is to solve the
linear part \(\ddot\theta + 2\gamma\dot\theta + \omega_0^2\theta = 0\)
*exactly* via a matrix exponential propagator, and treat the nonlinear residual
with a three-stage quadrature of Simpson type.

### Propagator Matrix

The state vector
\(\mathbf{x} = (\theta,\;\omega)^T\)
satisfies the homogeneous equation
\(\dot{\mathbf{x}} = A\,\mathbf{x}\)
with

\[
A = \begin{pmatrix} 0 & 1 \\ -\omega_0^2 & -2\gamma \end{pmatrix}.
\]

The matrix exponential \(\Phi(\tau) = e^{A\tau}\) is computed analytically in
three regimes depending on the discriminant
\(\Delta = \omega_0^2 - \gamma^2\):

**Underdamped** (\(\Delta > 0\)): let
\(\omega_d = \sqrt{\omega_0^2 - \gamma^2}\).

\[
\Phi(\tau) = e^{-\gamma\tau}
\begin{pmatrix}
\cos\omega_d\tau + \frac{\gamma}{\omega_d}\sin\omega_d\tau &
\frac{1}{\omega_d}\sin\omega_d\tau \\[4pt]
-\frac{\omega_0^2}{\omega_d}\sin\omega_d\tau &
\cos\omega_d\tau - \frac{\gamma}{\omega_d}\sin\omega_d\tau
\end{pmatrix}
\]

**Overdamped** (\(\Delta < 0\)): let
\(\omega_h = \sqrt{\gamma^2 - \omega_0^2}\).

\[
\Phi(\tau) = e^{-\gamma\tau}
\begin{pmatrix}
\cosh\omega_h\tau + \frac{\gamma}{\omega_h}\sinh\omega_h\tau &
\frac{1}{\omega_h}\sinh\omega_h\tau \\[4pt]
-\frac{\omega_0^2}{\omega_h}\sinh\omega_h\tau &
\cosh\omega_h\tau - \frac{\gamma}{\omega_h}\sinh\omega_h\tau
\end{pmatrix}
\]

**Critically damped** (\(\Delta = 0\)):

\[
\Phi(\tau) = e^{-\gamma\tau}
\begin{pmatrix}
1 + \gamma\tau & \tau \\
-\omega_0^2\,\tau & 1 - \gamma\tau
\end{pmatrix}
\]

### Constant-Residual Response

The particular solution for a constant unit forcing \(g = 1\) from zero initial
conditions is:

\[
\mathbf{r}(\tau) =
\begin{pmatrix}
r_\theta(\tau) \\ r_\omega(\tau)
\end{pmatrix}
=
\begin{pmatrix}
\displaystyle\frac{1 - \Phi_{22}(\tau) - 2\gamma\,\Phi_{12}(\tau)}{\omega_0^2} \\[6pt]
\Phi_{12}(\tau)
\end{pmatrix}
\]

which follows from the identity
\(\int_0^\tau \Phi_{12}(s)\,ds = r_\theta(\tau)\)
and the relation
\(r_\omega(\tau) = \dot{r}_\theta(\tau) = \Phi_{12}(\tau)\).

### Three-Stage Stepping Procedure

Given state \((\theta_n,\,\omega_n)\) at time \(t_n\):

1. **Stage 1**: evaluate residual at the current state
   \[g_1 = g(t_n,\;\theta_n,\;\omega_n)\]

2. **Stage 2**: propagate a half-step with constant forcing \(g_1\)
   \[
   \begin{pmatrix}\theta_2 \\ \omega_2\end{pmatrix}
   = \Phi\!\left(\tfrac{h}{2}\right)\begin{pmatrix}\theta_n \\ \omega_n\end{pmatrix}
   + g_1\,\mathbf{r}\!\left(\tfrac{h}{2}\right)
   \]
   then evaluate \(g_2 = g(t_n + h/2,\;\theta_2,\;\omega_2)\).

3. **Stage 3**: propagate a full step with constant forcing \(g_2\)
   \[
   \begin{pmatrix}\theta_3 \\ \omega_3\end{pmatrix}
   = \Phi(h)\begin{pmatrix}\theta_n \\ \omega_n\end{pmatrix}
   + g_2\,\mathbf{r}(h)
   \]
   then evaluate \(g_3 = g(t_n + h,\;\theta_3,\;\omega_3)\).

4. **Final combination** using Simpson-type weights \((1/6,\;2/3,\;1/6)\):
   \[
   \begin{aligned}
   \theta_{n+1} &= \Phi_{11}\,\theta_n + \Phi_{12}\,\omega_n
     + h\left[\Phi_{12}(h)\,\frac{g_1}{6}
       + \Phi_{12}\!\left(\tfrac{h}{2}\right)\frac{2g_2}{3}\right] \\[4pt]
   \omega_{n+1} &= \Phi_{21}\,\theta_n + \Phi_{22}\,\omega_n
     + h\left[\Phi_{22}(h)\,\frac{g_1}{6}
       + \Phi_{22}\!\left(\tfrac{h}{2}\right)\frac{2g_2}{3}
       + \frac{g_3}{6}\right]
   \end{aligned}
   \]

### Properties

| Property          | Value                                      |
|-------------------|--------------------------------------------|
| Order             | 3                                          |
| Stages            | 3 residual evaluations                     |
| Implicit          | No                                         |
| Exponential       | Yes (exact linear propagation)             |
| Suited for        | Damped oscillators with nonlinear forcing   |
| Local error       | \(\mathcal{O}(h^4)\)                      |

---

## Summary Table

| Method                   | Function                      | Order | Symplectic | Implicit | Family            |
|--------------------------|-------------------------------|:-----:|:----------:|:--------:|-------------------|
| RK3                      | `rk3_step`                    |   3   |     No     |    No    | Runge--Kutta      |
| RK4                      | `rk4_step`                    |   4   |     No     |    No    | Runge--Kutta      |
| RK5 (Butcher)            | `rk5_step`                    |   5   |     No     |    No    | Runge--Kutta      |
| Bogacki--Shampine (RK23) | `rk23_step`                   |   3   |     No     |    No    | Embedded RK       |
| RKF45                    | `rkf45_step`                  |   5   |     No     |    No    | Embedded RK       |
| Semi-Implicit Euler      | `semi_implicit_euler_step`    |   1   |    Yes     |  Semi    | Symplectic        |
| Leapfrog                 | `leapfrog_step`               |   2   |    Yes     |    No    | Symplectic        |
| Ruth4                    | `ruth4_step`                  |   4   |    Yes     |    No    | Symplectic        |
| Velocity Verlet          | `velocity_verlet_step`        |   2   |    Yes     |    No    | Verlet            |
| Runge--Kutta--Nyström    | `runge_kutta_nystrom_step`    |   4   |     No     |    No    | Nyström           |
| Numerov                  | `numerov_next_theta`          |   4   |     No     |   Yes    | Multistep         |
| DEN3                     | `den3_step`                   |   3   |     No     |    No    | Exponential       |
