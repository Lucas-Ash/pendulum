#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('damped_pendulum_data.dat', comments='#')
t = data[:,0]
theta_exact = data[:,1]
theta_lin = data[:,2]
theta_nl = data[:,3]
err_lin = data[:,4]
err_nl = data[:,5]
energy_nl = data[:,6]

fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
fig.suptitle('Damped Pendulum (theta0=0.3 rad, gamma=0.7, omega0=3.13209)', fontsize=14)

axes[0].plot(t, theta_exact, label='Analytical (linearized)', linewidth=2.2)
axes[0].plot(t, theta_lin, '--', label='RK4 linearized', linewidth=1.8)
axes[0].plot(t, theta_nl, '-.', label='RK4 nonlinear (sin theta)', linewidth=1.8)
axes[0].set_ylabel('theta (rad)')
axes[0].set_title('Angular Displacement')
axes[0].grid(True, alpha=0.3)
axes[0].legend(loc='best')

axes[1].semilogy(t, err_lin, label='|RK4 linear - analytical|', linewidth=1.6)
axes[1].semilogy(t, err_nl, label='|RK4 nonlinear - analytical|', linewidth=1.6)
axes[1].set_ylabel('|Error| (rad)')
axes[1].set_title('Absolute Error vs Analytical Solution')
axes[1].grid(True, alpha=0.3)
axes[1].legend(loc='best')

axes[2].plot(t, energy_nl, label='E(t) nonlinear', linewidth=1.6)
axes[2].set_xlabel('Time t (s)')
axes[2].set_ylabel('E / (m L^2)')
axes[2].set_title('Mechanical Energy (nonlinear pendulum)')
axes[2].grid(True, alpha=0.3)
axes[2].legend(loc='best')

plt.tight_layout()
plt.savefig('damped_pendulum.png', dpi=150, bbox_inches='tight')
print('Plot saved to damped_pendulum.png')
plt.show()
