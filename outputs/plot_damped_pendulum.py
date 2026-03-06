#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('outputs/damped_pendulum_1.dat', comments='#')
t = data[:,0]
theta_exact = data[:,1]
theta_nl = data[:,2]
omega_nl = data[:,3]
err_nl = data[:,4]
err_omega = data[:,5]
energy_nl = data[:,6]

fig, axes = plt.subplots(2, 2, figsize=(13, 10))
fig.suptitle('Damped Pendulum (theta0=0.3 rad, gamma=0.7)', fontsize=14)

ax_theta = axes[0, 0]
ax_err = axes[0, 1]
ax_energy = axes[1, 0]
ax_phase = axes[1, 1]

ax_theta.plot(t, theta_exact, label='Analytical (linearized)', linewidth=2.2)
ax_theta.plot(t, theta_nl, '-.', label='RK4 nonlinear (sin theta)', linewidth=1.8)
ax_theta.set_xlabel('Time t (s)')
ax_theta.set_ylabel('theta (rad)')
ax_theta.set_title('Angular Displacement')
ax_theta.grid(True, alpha=0.3)
ax_theta.legend(loc='best')

ax_err.semilogy(t, err_nl, label='|RK4 theta nonlinear - exp theta exact|', linewidth=1.6)
ax_err.set_xlabel('Time t (s)')
ax_err.set_ylabel('|Error| (rad)')
ax_err.set_title('Absolute Error vs Analytical Solution')
ax_err.grid(True, alpha=0.3)
ax_err.legend(loc='best')

ax_energy.plot(t, energy_nl, label='E(t) nonlinear', linewidth=1.6)
ax_energy.set_xlabel('Time t (s)')
ax_energy.set_ylabel('E / (m L^2)')
ax_energy.set_title('Mechanical Energy (nonlinear pendulum)')
ax_energy.grid(True, alpha=0.3)
ax_energy.legend(loc='best')

ax_phase.plot(theta_nl, omega_nl, label='Phase trajectory', linewidth=1.4)
ax_phase.plot([0.0], [0.0], 'ko', label='Equilibrium')
ax_phase.set_xlabel('theta (rad)')
ax_phase.set_ylabel('omega (rad/s)')
ax_phase.set_title('Phase Map (theta, omega)')
ax_phase.grid(True, alpha=0.3)
ax_phase.legend(loc='best')
ax_phase.set_aspect('equal', adjustable='box')

plt.tight_layout()
plt.savefig('outputs/damped_pendulum_1.png', dpi=150, bbox_inches='tight')
print('Plot saved to outputs/damped_pendulum_1.png')
plt.close(fig)
