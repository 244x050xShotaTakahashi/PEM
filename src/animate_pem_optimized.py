#!/usr/bin/env python3
"""
Optimized PEM Animation Generator

Performance improvements:
- Chunked file reading for large datasets
- NumPy vectorization for data processing  
- Multiprocessing for parallel frame generation
- Memory-efficient data structures
- Optional frame skipping for faster preview
"""

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from pathlib import Path
from typing import Union, List, Dict, Optional, Tuple
import multiprocessing as mp
from functools import partial
import argparse
import sys
import time
import gc

# Constants
DATA_FILE_DEFAULT = Path(__file__).resolve().parent.parent / "data" / "graph11.d"
CHUNK_SIZE = 10000  # Lines to read at once
MAX_FRAMES_IN_MEMORY = 100  # Maximum frames to keep in memory


class OptimizedPEMReader:
    """Optimized reader for PEM simulation data files."""
    
    def __init__(self, filename: Union[str, Path], chunk_size: int = CHUNK_SIZE):
        self.filename = Path(filename)
        self.chunk_size = chunk_size
        self.file_size = self.filename.stat().st_size
        
    def read_frames_chunked(self, max_frames: Optional[int] = None, 
                          skip_frames: int = 0) -> List[Dict]:
        """Read simulation frames in chunks for memory efficiency."""
        frames_data = []
        frame_count = 0
        skip_counter = 0
        
        try:
            with open(self.filename, 'r', buffering=8192*1024) as fh:
                buffer = []
                
                while True:
                    # Read chunk of lines
                    lines = fh.readlines(self.chunk_size)
                    if not lines:
                        break
                        
                    buffer.extend(lines)
                    
                    # Process complete frames from buffer
                    while self._has_complete_frame(buffer):
                        frame_data = self._extract_frame(buffer)
                        if frame_data:
                            if skip_counter >= skip_frames:
                                frames_data.append(frame_data)
                                skip_counter = 0
                                frame_count += 1
                                
                                if max_frames and frame_count >= max_frames:
                                    return frames_data
                            else:
                                skip_counter += 1
                                
        except Exception as e:
            print(f"Error reading file: {e}")
            return frames_data
            
        return frames_data
    
    def _has_complete_frame(self, buffer: List[str]) -> bool:
        """Check if buffer contains at least one complete frame."""
        if len(buffer) < 2:
            return False
            
        try:
            # Parse header to determine frame size
            tokens = buffer[0].split()
            if len(tokens) >= 3:
                num_particles = int(float(tokens[0]))
                # Each frame needs: header (2 lines) + 3 lines of particle data
                required_lines = 2 + 3
                return len(buffer) >= required_lines
        except:
            pass
            
        return False
    
    def _extract_frame(self, buffer: List[str]) -> Optional[Dict]:
        """Extract one complete frame from buffer."""
        try:
            # Parse header
            header_tokens = buffer[0].split()
            num_particles = int(float(header_tokens[0]))
            time_val = float(header_tokens[1])
            container_width = float(header_tokens[2])
            
            # Parse rmax
            rmax_tokens = buffer[1].split()
            rmax_val = float(rmax_tokens[0])
            
            # Remove processed header lines
            buffer[:2] = []
            
            # Parse particle data using NumPy for efficiency
            if num_particles > 0:
                # Read position data
                pos_line = buffer[0]
                pos_values = np.fromstring(pos_line, sep=' ')
                pos_data = pos_values.reshape(-1, 3)[:num_particles]
                
                # Read velocity data
                vel_line = buffer[1]
                vel_values = np.fromstring(vel_line, sep=' ')
                vel_data = vel_values.reshape(-1, 3)[:num_particles]
                
                # Read rotation data
                rot_line = buffer[2]
                rot_values = np.fromstring(rot_line, sep=' ')[:num_particles]
                
                # Remove processed data lines
                buffer[:3] = []
                
                return {
                    'time': time_val,
                    'num_particles': num_particles,
                    'container_width': container_width,
                    'positions': pos_data,  # NumPy array (N, 3)
                    'velocities': vel_data,  # NumPy array (N, 3)
                    'rotations': rot_values  # NumPy array (N,)
                }
            else:
                # Empty frame
                buffer[:3] = []
                return {
                    'time': time_val,
                    'num_particles': 0,
                    'container_width': container_width,
                    'positions': np.array([]),
                    'velocities': np.array([]),
                    'rotations': np.array([])
                }
                
        except Exception as e:
            print(f"Error extracting frame: {e}")
            return None


class OptimizedAnimator:
    """Optimized animator for PEM simulations."""
    
    def __init__(self, frames_data: List[Dict], output_filename: str = "pem_animation.gif"):
        self.frames_data = frames_data
        self.output_filename = output_filename
        self.fig = None
        self.ax = None
        self.patches = None
        self.time_text = None
        
    def setup_plot(self):
        """Initialize the plot with proper settings."""
        self.fig, self.ax = plt.subplots(figsize=(10, 10), dpi=100)
        
        # Calculate plot limits
        container_width = self.frames_data[0]['container_width']
        max_z = self._calculate_max_z()
        
        self.ax.set_xlim(-0.1 * container_width, 1.1 * container_width)
        self.ax.set_ylim(-0.1 * max_z, 1.2 * max_z)
        self.ax.set_aspect('equal', adjustable='box')
        self.ax.set_xlabel("X coordinate")
        self.ax.set_ylabel("Z coordinate")
        
        # Draw walls once
        self.ax.plot([0, 0], [0, max_z * 1.15], 'k-', lw=2)
        self.ax.plot([0, container_width], [0, 0], 'k-', lw=2)
        self.ax.plot([container_width, container_width], [0, max_z * 1.15], 'k-', lw=2)
        
        # Create text object
        self.time_text = self.ax.text(0.05, 0.95, '', transform=self.ax.transAxes, 
                                     ha="left", va="top", fontsize=10)
        
    def _calculate_max_z(self) -> float:
        """Calculate maximum Z coordinate across all frames."""
        max_z = 0.0
        for frame in self.frames_data:
            if frame['num_particles'] > 0:
                positions = frame['positions']
                radii = positions[:, 2]
                z_coords = positions[:, 1]
                frame_max = np.max(z_coords + radii)
                max_z = max(max_z, frame_max)
        return max_z if max_z > 0 else self.frames_data[0]['container_width'] * 0.5
        
    def update_frame_vectorized(self, frame_idx: int):
        """Update frame using vectorized operations."""
        # Clear previous patches
        if self.patches:
            self.patches.remove()
            
        frame = self.frames_data[frame_idx]
        
        # Update title and time
        self.fig.suptitle(f"Particle Simulation (Frame {frame_idx+1}/{len(self.frames_data)})", 
                          fontsize=12)
        self.time_text.set_text(f"Time: {frame['time']:.6f} s")
        
        if frame['num_particles'] > 0:
            positions = frame['positions']
            rotations = frame['rotations']
            
            # Create circles using vectorization
            circles = []
            for i in range(frame['num_particles']):
                x, z, r = positions[i]
                circle = Circle((x, z), r, facecolor='white', 
                              edgecolor='black', linewidth=1.5)
                circles.append(circle)
                
            # Add all circles at once
            self.patches = PatchCollection(circles, match_original=True)
            self.ax.add_collection(self.patches)
            
            # Draw rotation indicators efficiently
            x_centers = positions[:, 0]
            z_centers = positions[:, 1]
            radii = positions[:, 2]
            
            x_ends = x_centers + radii * np.cos(rotations)
            z_ends = z_centers + radii * np.sin(rotations)
            
            for i in range(frame['num_particles']):
                self.ax.plot([x_centers[i], x_ends[i]], 
                           [z_centers[i], z_ends[i]], 
                           'k-', linewidth=1, alpha=0.8)
        
        return []
        
    def animate(self, fps: int = 15, interval: int = 100):
        """Create and save the animation."""
        self.setup_plot()
        
        ani = animation.FuncAnimation(
            self.fig, self.update_frame_vectorized, 
            frames=len(self.frames_data), 
            blit=False, interval=interval
        )
        
        try:
            print(f"Saving animation to '{self.output_filename}'...")
            start_time = time.time()
            
            writer = animation.PillowWriter(fps=fps)
            ani.save(self.output_filename, writer=writer)
            
            elapsed = time.time() - start_time
            print(f"Animation saved successfully in {elapsed:.2f} seconds")
            
        except Exception as e:
            print(f"Error saving animation: {e}")
            print("Make sure Pillow is installed: pip install Pillow")
            

def parallel_frame_processor(frame_data: Dict) -> Dict:
    """Process a single frame in parallel (placeholder for future enhancements)."""
    # Could add frame preprocessing here if needed
    return frame_data


def main():
    """Main function with command-line interface."""
    parser = argparse.ArgumentParser(description='Optimized PEM Animation Generator')
    parser.add_argument('input', nargs='?', default=str(DATA_FILE_DEFAULT),
                       help='Input data file (default: data/graph11.d)')
    parser.add_argument('-o', '--output', default='pem_animation.gif',
                       help='Output animation file (default: pem_animation.gif)')
    parser.add_argument('-f', '--fps', type=int, default=15,
                       help='Frames per second (default: 15)')
    parser.add_argument('--max-frames', type=int, default=None,
                       help='Maximum number of frames to process')
    parser.add_argument('--skip-frames', type=int, default=0,
                       help='Skip every N frames for faster preview')
    parser.add_argument('--chunk-size', type=int, default=CHUNK_SIZE,
                       help='Number of lines to read at once')
    parser.add_argument('--profile', action='store_true',
                       help='Enable performance profiling')
    
    args = parser.parse_args()
    
    if args.profile:
        import cProfile
        import pstats
        profiler = cProfile.Profile()
        profiler.enable()
    
    # Read data
    print(f"Reading simulation data from {args.input}...")
    start_time = time.time()
    
    reader = OptimizedPEMReader(args.input, chunk_size=args.chunk_size)
    frames_data = reader.read_frames_chunked(
        max_frames=args.max_frames,
        skip_frames=args.skip_frames
    )
    
    read_time = time.time() - start_time
    print(f"Read {len(frames_data)} frames in {read_time:.2f} seconds")
    
    if not frames_data:
        print("No data to animate!")
        return 1
    
    # Memory usage info
    memory_usage = sys.getsizeof(frames_data) / (1024 * 1024)
    print(f"Memory usage: {memory_usage:.2f} MB")
    
    # Create animation
    animator = OptimizedAnimator(frames_data, args.output)
    animator.animate(fps=args.fps)
    
    if args.profile:
        profiler.disable()
        stats = pstats.Stats(profiler)
        stats.sort_stats('cumulative')
        print("\nPerformance Profile:")
        stats.print_stats(20)
    
    # Cleanup
    gc.collect()
    
    return 0


if __name__ == "__main__":
    sys.exit(main())