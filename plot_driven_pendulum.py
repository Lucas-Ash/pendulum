#!/usr/bin/env python3
"""
Pendulum Simulation Visualization
Compares Numerical (RK4) vs Analytical (Small-Angle) Solutions
"""

import numpy as np
import matplotlib.pyplot as plt

# Load data from CSV
data = np.loadtxt('driven_pendulum_data.csv', delimiter=',', skiprows=1)
time = data[:, 0]
theta_numerical = data[:, 1]
theta_analytical = data[:, 2]
difference = data[:, 3]

# Create figure with subplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

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

# Add parameter annotation
params_text = 'Parameters:\ng = 9.81 m/s^2\nL = 1 m\ndamping = 0.5\nA = 0.5\nomega_drive = 1.2 rad/s\ntheta0 = 0.1 rad\nomega0 = 0 rad/s'
fig.text(0.02, 0.02, params_text, fontsize=9, verticalalignment='bottom',
         fontfamily='monospace', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
plt.subplots_adjust(bottom=0.15)

# Save figure
plt.savefig('driven_pendulum_comparison.png', dpi=150, bbox_inches='tight')
