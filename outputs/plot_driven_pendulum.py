#!/usr/bin/env python3
"""
Pendulum Simulation Visualization
Compares Numerical (RK4) vs Analytical (Small-Angle) Solutions
"""

import numpy as np
import matplotlib.pyplot as plt

# Load data from CSV
data = np.loadtxt('outputs/driven_pendulum_data.csv', delimiter=',', skiprows=1)
time = data[:, 0]
theta_numerical = data[:, 2]
omega_numerical = data[:, 3]
theta_analytical = data[:, 1]
difference = data[:, 4]

# Create figure with subplots
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 12))

# Plot 1: Theta vs Time
ax1.plot(time, theta_numerical, 'b-', linewidth=1.5, label='Numerical (RK4)', alpha=0.8)
ax1.plot(time, theta_analytical, 'r--', linewidth=1.5, label='Analytical (Small-Angle)', alpha=0.8)
ax1.set_xlabel('Time (s)', fontsize=12)
ax1.set_ylabel('Angle theta (rad)', fontsize=12)
ax1.set_title('Damped Driven Pendulum: Numerical vs Analytical Solution', fontsize=14, fontweight='bold')
ax1.legend(fontsize=11)
ax1.grid(True, alpha=0.3)
ax1.set_xlim([0, 50])

# Plot 2: Difference between solutions
ax2.plot(time, difference, 'g-', linewidth=1, label='|Numerical - Analytical|')
ax2.set_xlabel('Time (s)', fontsize=12)
ax2.set_ylabel('Absolute Difference (rad)', fontsize=12)
ax2.set_title('Difference Between Numerical and Analytical Solutions', fontsize=14, fontweight='bold')
ax2.legend(fontsize=11)
ax2.grid(True, alpha=0.3)
ax2.set_xlim([0, 50])

# Plot 3: Phase map
ax3.plot(theta_numerical, omega_numerical, 'm-', linewidth=1.2, label='Phase trajectory')
ax3.plot([0.0], [0.0], 'ko', label='Equilibrium')
ax3.set_xlabel('theta (rad)', fontsize=12)
ax3.set_ylabel('omega (rad/s)', fontsize=12)
ax3.set_title('Phase Map (theta, omega)', fontsize=14, fontweight='bold')
ax3.legend(fontsize=11)
ax3.grid(True, alpha=0.3)
ax3.set_aspect('equal', adjustable='box')

# Add parameter annotation
params_text = 'Parameters:\ng = 9.81 m/s^2\nL = 1 m\ndamping = 0.5\nA = 0.5\nomega_drive = 1.2 rad/s\ntheta0 = 0.1 rad\nomega0 = 0 rad/s'
fig.text(0.02, 0.02, params_text, fontsize=9, verticalalignment='bottom',
         fontfamily='monospace', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
plt.subplots_adjust(bottom=0.15)

# Save figure
plt.savefig('outputs/driven_pendulum_comparison.png', dpi=150, bbox_inches='tight')
