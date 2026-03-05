#!/usr/bin/env python3
"""
Pendulum Simulation Video Generator

This script analyzes pendulum simulation data from .dat and .csv files in QA/**/outputs
directories and generates animated plot videos. The videos play in real-time seconds,
with duration constraints of 5-20 seconds for simulations outside this range.

Usage:
    python analyze_and_plot_videos.py
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pathlib import Path
from typing import Tuple, Dict, List


class PendulumDataAnalyzer:
    """Analyzes and visualizes pendulum simulation data."""
    
    def __init__(self, qa_dir: str = "QA"):
        self.qa_dir = qa_dir
        self.output_video_dir = os.path.join(qa_dir, "video_outputs")
        os.makedirs(self.output_video_dir, exist_ok=True)
    
    def find_data_files(self) -> List[str]:
        """Find all .dat and .csv files in QA/**/outputs directories."""
        patterns = [
            os.path.join(self.qa_dir, "**", "outputs", "*.dat"),
            os.path.join(self.qa_dir, "**", "outputs", "*.csv")
        ]
        files = []
        for pattern in patterns:
            files.extend(glob.glob(pattern, recursive=True))
        return sorted(files)
    
    def load_data(self, filepath: str) -> Tuple[np.ndarray | None, np.ndarray | None, np.ndarray | None, str | None]:
        """
        Load data from .dat or .csv file.
        
        Returns:
            Tuple of (time, theta, omega, pendulum_type)
        """
        ext = os.path.splitext(filepath)[1]
        
        # Determine pendulum type from path
        if "simple" in filepath.lower():
            pendulum_type = "Simple Pendulum"
        elif "damped" in filepath.lower():
            pendulum_type = "Damped Pendulum"
        elif "driven" in filepath.lower():
            pendulum_type = "Driven Pendulum"
        else:
            pendulum_type = "Pendulum"
        
        try:
            if ext == ".dat":
                # Load .dat file (space-separated, may have comments)
                data = np.loadtxt(filepath, comments="#")
                time = data[:, 0]
                theta = data[:, 2]  # Column 2 is theta (numerical)
                omega = data[:, 3]  # Column 3 is omega
            else:  # .csv
                # Load .csv file
                df = pd.read_csv(filepath)
                time = df['Time'].values
                theta = df['Theta'].values
                omega = df['Omega'].values
            
            return time, theta, omega, pendulum_type
        
        except Exception as e:
            print(f"Error loading {filepath}: {e}")
            return None, None, None, None
    
    def calculate_video_params(self, sim_duration: float) -> Tuple[float, int, float]:
        """
        Calculate video duration and FPS based on simulation duration.
        
        Args:
            sim_duration: Total simulation time in seconds
        
        Returns:
            Tuple of (video_duration, fps)
        """
        # Target: real-time playback, but constrain to 5-20 seconds
        if sim_duration < 5:
            video_duration = 5.0
        elif sim_duration > 20:
            video_duration = 20.0
        else:
            video_duration = sim_duration
        
        # Calculate playback speed factor
        speed_factor = sim_duration / video_duration
        
        # Use 30 FPS for smooth animation
        fps = 30
        
        return video_duration, fps, speed_factor
    
    def create_animation(self, filepath: str):
        """Create and save animation for a data file."""
        print(f"\nProcessing: {filepath}")
        
        # Load data
        time, theta, omega, pendulum_type = self.load_data(filepath)
        
        if time is None:
            print(f"Skipping {filepath} due to loading error")
            return
        
        sim_duration = time[-1] - time[0]
        print(f"  Simulation duration: {sim_duration:.2f} seconds")
        
        # Calculate video parameters
        video_duration, fps, speed_factor = self.calculate_video_params(sim_duration)
        print(f"  Video duration: {video_duration:.2f} seconds")
        print(f"  Playback speed: {speed_factor:.2f}x")
        print(f"  FPS: {fps}")
        
        # Calculate total frames
        total_frames = int(video_duration * fps)
        
        # Create figure with subplots
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12))
        fig.suptitle(f"{pendulum_type} - {os.path.basename(filepath)}", fontsize=14, fontweight='bold')
        
        # Initialize plots
        line1, = ax1.plot([], [], 'b-', linewidth=2, label='Theta')
        ax1.set_xlabel('Time (s)')
        ax1.set_ylabel('Theta (rad)')
        ax1.set_title('Angular Position vs Time')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        
        line2, = ax2.plot([], [], 'r-', linewidth=2, label='Omega')
        ax2.set_xlabel('Time (s)')
        ax2.set_ylabel('Omega (rad/s)')
        ax2.set_title('Angular Velocity vs Time')
        ax2.grid(True, alpha=0.3)
        ax2.legend()
        
        line3, = ax3.plot([], [], 'g-', linewidth=2)
        ax3.set_xlabel('Theta (rad)')
        ax3.set_ylabel('Omega (rad/s)')
        ax3.set_title('Phase Space (Omega vs Theta)')
        ax3.grid(True, alpha=0.3)
        
        # Set axis limits
        ax1.set_xlim(time[0], time[-1])
        ax1.set_ylim(np.min(theta) * 1.1, np.max(theta) * 1.1)
        
        ax2.set_xlim(time[0], time[-1])
        ax2.set_ylim(np.min(omega) * 1.1, np.max(omega) * 1.1)
        
        ax3.set_xlim(np.min(theta) * 1.1, np.max(theta) * 1.1)
        ax3.set_ylim(np.min(omega) * 1.1, np.max(omega) * 1.1)
        
        # Time text
        time_text = ax1.text(0.02, 0.95, '', transform=ax1.transAxes,
                            verticalalignment='top', fontsize=10,
                            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        plt.tight_layout()
        
        def init():
            """Initialize animation."""
            line1.set_data([], [])
            line2.set_data([], [])
            line3.set_data([], [])
            time_text.set_text('')
            return line1, line2, line3, time_text
        
        def animate(frame):
            """Animation function."""
            # Calculate which data index to show based on frame
            data_idx = int((frame / total_frames) * len(time))
            data_idx = min(data_idx, len(time) - 1)
            
            # Update plots with data up to current index
            line1.set_data(time[:data_idx+1], theta[:data_idx+1])
            line2.set_data(time[:data_idx+1], omega[:data_idx+1])
            line3.set_data(theta[:data_idx+1], omega[:data_idx+1])
            
            # Update time text
            current_time = time[data_idx]
            time_text.set_text(f'Time: {current_time:.3f} s\nSpeed: {speed_factor:.2f}x')
            
            return line1, line2, line3, time_text
        
        # Create animation
        anim = animation.FuncAnimation(fig, animate, init_func=init,
                                      frames=total_frames, interval=1000/fps,
                                      blit=True, repeat=True)
        
        # Generate output filename
        base_name = os.path.splitext(os.path.basename(filepath))[0]
        output_path = os.path.join(self.output_video_dir, f"{base_name}_animation.mp4")
        
        # Save animation
        print(f"  Saving video to: {output_path}")
        try:
            Writer = animation.writers['ffmpeg']
            writer = Writer(fps=fps, metadata=dict(artist='Pendulum Analyzer'), bitrate=1800)
            anim.save(output_path, writer=writer)
            print(f"  ✓ Video saved successfully!")
        except Exception as e:
            print(f"  ✗ Error saving video: {e}")
            print(f"  Note: Make sure ffmpeg is installed (sudo apt-get install ffmpeg)")
        
        plt.close(fig)
    
    def process_all_files(self):
        """Process all data files and create animations."""
        files = self.find_data_files()
        
        if not files:
            print("No data files found in QA/**/outputs directories")
            return
        
        print(f"Found {len(files)} data file(s):")
        for f in files:
            print(f"  - {f}")
        
        print(f"\nOutput directory: {self.output_video_dir}")
        print("=" * 60)
        
        for filepath in files:
            self.create_animation(filepath)
        
        print("\n" + "=" * 60)
        print(f"Processing complete! Videos saved to: {self.output_video_dir}")


def main():
    """Main entry point."""
    analyzer = PendulumDataAnalyzer()
    analyzer.process_all_files()


if __name__ == "__main__":
    main()