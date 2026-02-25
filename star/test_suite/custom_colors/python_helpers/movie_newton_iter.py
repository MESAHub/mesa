#!/usr/bin/env python3
"""
plot_newton_iter_video.py — Video plotter for MESA Colors per-iteration output

Creates an animated video showing Newton iterations accumulating over time.
Points appear one by one, sorted by model number then iteration number.
History file points appear incrementally as each timestep completes.

Author: Niall Miller (2025)
"""

import argparse
import os
import re
import sys
from typing import List, Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter
from matplotlib.colors import Normalize
from matplotlib import cm

# Import shared functionality from static plotter
from plot_newton_iter import (
    # Terminal UI
    use_color, term_width,
    BOLD, DIM, CYAN, YELL, GREEN, RED, RESET,
    print_header, print_subheader, print_success, print_error, print_info,
    # Data loading
    load_iteration_data,
    load_history_data,
    resolve_history_axis,
    MESA_READER_AVAILABLE,
    # Expression parsing
    is_expression, parse_expression, resolve_axis,
    # Prompting
    prompt_yes_no,
)


# ============================================================================
#                              SIGMA CLIPPING
# ============================================================================

def sigma_clip_mask(
    data_arrays: List[np.ndarray],
    sigma: float = 3.0,
    max_iter: int = 5
) -> np.ndarray:
    """
    Create a mask identifying outliers using iterative sigma clipping.
    
    Clips based on ALL provided arrays - a point is kept only if it's
    within sigma of the median in ALL arrays.
    """
    n_points = len(data_arrays[0])
    mask = np.ones(n_points, dtype=bool)
    
    for iteration in range(max_iter):
        prev_mask = mask.copy()
        
        for arr in data_arrays:
            valid_data = arr[mask]
            if len(valid_data) < 3:
                continue
            
            median = np.median(valid_data)
            mad = np.median(np.abs(valid_data - median))
            std_equiv = mad * 1.4826
            
            if std_equiv > 0:
                deviation = np.abs(arr - median)
                mask &= (deviation < sigma * std_equiv)
        
        if np.array_equal(mask, prev_mask):
            break
    
    return mask


# ============================================================================
#                           VIDEO-SPECIFIC DATA HANDLING
# ============================================================================

def sort_by_model_and_iter(
    data: np.ndarray, 
    column_names: List[str]
) -> Tuple[np.ndarray, np.ndarray]:
    """Sort data by model number, then by iteration number."""
    model_idx = column_names.index('model')
    iter_idx = column_names.index('iter')
    sort_indices = np.lexsort((data[:, iter_idx], data[:, model_idx]))
    return data[sort_indices], sort_indices


def get_model_completion_indices(model_data: np.ndarray) -> np.ndarray:
    """
    Find the index where each model's iterations end.
    Returns array of indices - the last point for each unique model.
    """
    indices = []
    unique_models = np.unique(model_data)
    
    for model in unique_models:
        model_mask = model_data == model
        model_indices = np.where(model_mask)[0]
        # Last index for this model
        indices.append(model_indices[-1])
    
    return np.array(sorted(indices))


# ============================================================================
#                           VIDEO CREATION
# ============================================================================

def create_video(
    x_data: np.ndarray,
    y_data: np.ndarray,
    color_data: np.ndarray,
    model_data: np.ndarray,
    iter_data: np.ndarray,
    x_label: str,
    y_label: str,
    color_label: str,
    output_path: str,
    duration: float = 30.0,
    fps: int = 30,
    cmap: str = "viridis",
    point_size: int = 20,
    alpha: float = 0.7,
    flip_y: bool = False,
    dpi: int = 150,
    sigma_clip: float = 3.0,
    history_x: Optional[np.ndarray] = None,
    history_y: Optional[np.ndarray] = None,
) -> None:
    """
    Create animated video of points appearing one by one.
    
    History points appear incrementally as each model's iterations complete.
    Axis limits dynamically expand with smooth easing.
    """
    
    # =========================================
    # SIGMA CLIPPING - Remove outliers
    # =========================================
    print_info(f"Applying {sigma_clip}-sigma clipping to remove outliers...")
    
    valid_mask = sigma_clip_mask([x_data, y_data], sigma=sigma_clip)
    n_removed = np.sum(~valid_mask)
    
    if n_removed > 0:
        print_info(f"Removed {n_removed} outlier points ({100*n_removed/len(x_data):.1f}%)")
        x_data = x_data[valid_mask]
        y_data = y_data[valid_mask]
        color_data = color_data[valid_mask]
        model_data = model_data[valid_mask]
        iter_data = iter_data[valid_mask]
    else:
        print_info("No outliers detected")
    
    n_points = len(x_data)
    
    if n_points == 0:
        print_error("No points remaining after sigma clipping!")
        return
    
    # Find where each model ends (for incremental history plotting)
    model_completion_indices = get_model_completion_indices(model_data)
    n_models = len(model_completion_indices)
    
    # Check history data availability
    has_history = (history_x is not None and history_y is not None 
                   and len(history_x) >= n_models)
    
    if has_history:
        print_info(f"History data: {len(history_x)} points available, {n_models} will be shown incrementally")
        # Truncate history to match number of models in iteration data
        history_x = history_x[:n_models]
        history_y = history_y[:n_models]
    else:
        print_info("No matching history data for incremental overlay")
    
    # Calculate frames and timing
    total_frames = int(duration * fps)
    
    if n_points <= total_frames:
        points_per_frame = 1
        total_frames = n_points
    else:
        points_per_frame = max(1, n_points // total_frames)
        total_frames = (n_points + points_per_frame - 1) // points_per_frame
    
    actual_duration = total_frames / fps
    
    print_info(f"Total points: {n_points}")
    print_info(f"Points per frame: {points_per_frame}")
    print_info(f"Total frames: {total_frames}")
    print_info(f"Video duration: {actual_duration:.1f}s at {fps} fps")
    
    # =========================================
    # PRECOMPUTE AXIS LIMITS
    # =========================================
    print_info("Precomputing dynamic axis limits...")
    
    # Final limits including history data
    all_x = [x_data]
    all_y = [y_data]
    if has_history:
        all_x.append(history_x)
        all_y.append(history_y)
    
    all_x = np.concatenate(all_x)
    all_y = np.concatenate(all_y)
    
    x_min_final, x_max_final = np.min(all_x), np.max(all_x)
    y_min_final, y_max_final = np.min(all_y), np.max(all_y)
    
    x_range_final = x_max_final - x_min_final
    y_range_final = y_max_final - y_min_final
    
    pad_frac = 0.08
    x_min_final -= x_range_final * pad_frac
    x_max_final += x_range_final * pad_frac
    y_min_final -= y_range_final * pad_frac
    y_max_final += y_range_final * pad_frac
    
    # =========================================
    # SET UP FIGURE
    # =========================================
    fig, ax = plt.subplots(figsize=(10, 8))
    
    ax.set_xlabel(x_label, fontsize=12)
    ax.set_ylabel(y_label, fontsize=12)
    ax.grid(True, alpha=0.3)
    
    norm = Normalize(vmin=np.min(color_data), vmax=np.max(color_data))
    colormap = cm.get_cmap(cmap)
    
    sm = cm.ScalarMappable(cmap=colormap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label(color_label, fontsize=12)
    
    title = ax.set_title("", fontsize=14, fontweight='bold')
    
    # Initialize scatter plots
    scatter = ax.scatter([], [], c=[], cmap=cmap, norm=norm, s=point_size, alpha=alpha)
    scatter_history = ax.scatter([], [], c='black', marker='x', s=point_size*0.5, 
                                  linewidths=2, label='History', zorder=10)
    ax.legend(loc='upper right')
    
    plt.tight_layout()
    
    def get_dynamic_limits(n_show: int, n_history_show: int) -> Tuple[float, float, float, float]:
        """Compute axis limits that smoothly expand to contain all visible points."""
        x_visible = list(x_data[:n_show])
        y_visible = list(y_data[:n_show])
        
        if has_history and n_history_show > 0:
            x_visible.extend(history_x[:n_history_show])
            y_visible.extend(history_y[:n_history_show])
        
        x_visible = np.array(x_visible)
        y_visible = np.array(y_visible)
        
        x_min_curr = np.min(x_visible)
        x_max_curr = np.max(x_visible)
        y_min_curr = np.min(y_visible)
        y_max_curr = np.max(y_visible)
        
        x_range_curr = max(x_max_curr - x_min_curr, x_range_final * 0.01)
        y_range_curr = max(y_max_curr - y_min_curr, y_range_final * 0.01)
        
        x_min_padded = x_min_curr - x_range_curr * pad_frac
        x_max_padded = x_max_curr + x_range_curr * pad_frac
        y_min_padded = y_min_curr - y_range_curr * pad_frac
        y_max_padded = y_max_curr + y_range_curr * pad_frac
        
        # Progress and easing
        progress = n_show / n_points
        ease = 1 - (1 - progress) ** 2
        
        x_min = x_min_padded + (x_min_final - x_min_padded) * ease
        x_max = x_max_padded + (x_max_final - x_max_padded) * ease
        y_min = y_min_padded + (y_min_final - y_min_padded) * ease
        y_max = y_max_padded + (y_max_final - y_max_padded) * ease
        
        return x_min, x_max, y_min, y_max
    
    def init():
        scatter.set_offsets(np.empty((0, 2)))
        scatter.set_array(np.array([]))
        scatter_history.set_offsets(np.empty((0, 2)))
        ax.set_xlim(x_min_final, x_max_final)
        if flip_y:
            ax.set_ylim(y_max_final, y_min_final)
        else:
            ax.set_ylim(y_min_final, y_max_final)
        return scatter, scatter_history, title
    
    def update(frame):
        n_show = min((frame + 1) * points_per_frame, n_points)
        
        # Update Newton iteration scatter
        x_show = x_data[:n_show]
        y_show = y_data[:n_show]
        c_show = color_data[:n_show]
        
        scatter.set_offsets(np.column_stack([x_show, y_show]))
        scatter.set_array(c_show)
        
        # Count how many models have completed (for incremental history)
        n_history_show = 0
        if has_history:
            # A model is complete when we've shown PAST its completion index
            # (history X appears one point after the model's last iteration)
            n_history_show = np.sum(model_completion_indices < n_show - 1)
            
            if n_history_show > 0:
                scatter_history.set_offsets(np.column_stack([
                    history_x[:n_history_show],
                    history_y[:n_history_show]
                ]))
            else:
                scatter_history.set_offsets(np.empty((0, 2)))
        
        # Update dynamic axis limits
        x_min, x_max, y_min, y_max = get_dynamic_limits(n_show, n_history_show)
        ax.set_xlim(x_min, x_max)
        if flip_y:
            ax.set_ylim(y_max, y_min)
        else:
            ax.set_ylim(y_min, y_max)
        
        # Update title
        current_model = int(model_data[n_show - 1])
        current_iter = int(iter_data[n_show - 1])
        title.set_text(f"Newton Iteration Colors — Model {current_model}, Iter {current_iter}\n"
                      f"Points: {n_show}/{n_points} | Timesteps: {n_history_show}")
        
        return scatter, scatter_history, title
    
    # Create animation
    print_info("Generating animation frames...")
    anim = FuncAnimation(
        fig, update, frames=total_frames,
        init_func=init, blit=False, interval=1000/fps
    )
    
    # Save video
    print_info(f"Encoding video to {output_path}...")
    
    if output_path.endswith('.gif'):
        writer = PillowWriter(fps=fps)
    else:
        try:
            writer = FFMpegWriter(fps=fps, metadata={'title': 'Newton Iterations'}, bitrate=5000)
        except Exception:
            print_info("FFmpeg not available, falling back to GIF output")
            output_path = output_path.rsplit('.', 1)[0] + '.gif'
            writer = PillowWriter(fps=fps)
    
    anim.save(output_path, writer=writer, dpi=dpi)
    plt.close(fig)
    
    print_success(f"Saved: {output_path}")


# ============================================================================
#                           INTERACTIVE PROMPTS
# ============================================================================

def prompt_axis_choice(
    column_names: List[str],
    data: np.ndarray,
    label: str,
) -> Optional[Tuple[np.ndarray, str]]:
    """Prompt user to select an axis column or expression."""
    N = len(column_names)
    
    print_subheader(f"{label} ({CYAN}{N}{RESET} columns)")
    
    col_width = max(len(s) for s in column_names) + 8
    cols = max(1, min(3, term_width() // col_width))
    
    for i, name in enumerate(column_names):
        end = '\n' if (i + 1) % cols == 0 else ''
        print(f"  [{GREEN}{i:2d}{RESET}] {name:<{col_width-8}}", end=end)
    if N % cols != 0:
        print()
    
    print(f"\n{DIM}Enter: column number | column name | expression (e.g. B-R, [9]-[13], Teff/1000){RESET}")
    
    while True:
        inp = input(f"\n{CYAN}>{RESET} ").strip()
        
        if not inp:
            continue
        
        if inp.lower() == 'q':
            return None
        
        try:
            arr, lbl = resolve_axis(inp, column_names, data)
            return arr, lbl
        except ValueError as e:
            print_error(str(e))


def prompt_duration(n_points: int, default: float = 30.0) -> Tuple[float, int]:
    """Prompt user for video duration, showing points info."""
    print_subheader("Video Duration")
    print_info(f"Total data points: {n_points}")
    
    while True:
        inp = input(f"Video duration in seconds {DIM}[{default}]{RESET}: ").strip()
        
        if not inp:
            duration = default
            break
        
        try:
            duration = float(inp)
            if duration <= 0:
                print_error("Duration must be positive")
                continue
            break
        except ValueError:
            print_error("Invalid number")
    
    default_fps = 30
    while True:
        inp = input(f"Frames per second {DIM}[{default_fps}]{RESET}: ").strip()
        
        if not inp:
            fps = default_fps
            break
        
        try:
            fps = int(inp)
            if fps <= 0:
                print_error("FPS must be positive")
                continue
            break
        except ValueError:
            print_error("Invalid number")
    
    total_frames = int(duration * fps)
    points_per_frame = max(1, n_points // total_frames)
    points_per_second = points_per_frame * fps
    
    print_info(f"≈ {points_per_second} points/second ({points_per_frame} points/frame)")
    
    return duration, fps


# ============================================================================
#                           MAIN WORKFLOWS
# ============================================================================

def run_interactive(filepath: str, history_file: str = "../LOGS/history.data") -> None:
    """Run interactive mode."""
    print_header("MESA Colors — Newton Iteration Video Maker")
    
    print_info(f"Loading: {filepath}")
    try:
        column_names, data = load_iteration_data(filepath)
    except Exception as e:
        print_error(f"Failed to load file: {e}")
        sys.exit(1)
    
    if 'model' not in column_names or 'iter' not in column_names:
        print_error("Data must have 'model' and 'iter' columns")
        sys.exit(1)
    
    data, _ = sort_by_model_and_iter(data, column_names)
    print_success(f"Loaded {data.shape[0]} data points, {data.shape[1]} columns")
    print_success("Data sorted by model number, then iteration")
    
    # Load history file
    md = load_history_data(history_file)
    
    model_idx = column_names.index('model')
    iter_idx = column_names.index('iter')
    model_data = data[:, model_idx]
    iter_data = data[:, iter_idx]
    
    # Select X axis
    print(f"\n{DIM}Select X-axis (or enter expression like 'B-R'):{RESET}")
    result = prompt_axis_choice(column_names, data, "X-axis")
    if result is None:
        return
    x_data, x_label = result
    print_success(f"X-axis: {x_label}")
    
    # Select Y axis
    print(f"\n{DIM}Select Y-axis:{RESET}")
    result = prompt_axis_choice(column_names, data, "Y-axis")
    if result is None:
        return
    y_data, y_label = result
    print_success(f"Y-axis: {y_label}")
    
    flip_y = prompt_yes_no("Flip Y axis?", default=False)
    
    # Select color axis
    print(f"\n{DIM}Select Color axis:{RESET}")
    result = prompt_axis_choice(column_names, data, "Color")
    if result is None:
        return
    color_data, color_label = result
    print_success(f"Color: {color_label}")
    
    # Prompt for duration
    duration, fps = prompt_duration(len(x_data))
    
    # Output filename
    safe_x = re.sub(r'[^\w\-]', '_', x_label)
    safe_y = re.sub(r'[^\w\-]', '_', y_label)
    default_output = f"newton_iter_{safe_y}_vs_{safe_x}.mp4"
    
    out_inp = input(f"Output filename {DIM}[{default_output}]{RESET}: ").strip()
    output_path = out_inp if out_inp else default_output
    
    # Sigma clipping
    sigma_inp = input(f"Sigma clipping threshold {DIM}[3.0]{RESET}: ").strip()
    try:
        sigma_clip = float(sigma_inp) if sigma_inp else 3.0
    except ValueError:
        sigma_clip = 3.0
        print_info("Invalid input, using default 3.0")
    
    # Get history data for overlay
    history_x, history_y = None, None
    if md is not None:
        history_x = resolve_history_axis(x_label, md)
        history_y = resolve_history_axis(y_label, md)
        
        if history_x is not None and history_y is not None:
            print_info(f"History data: {len(history_x)} points available")
        else:
            missing = []
            if history_x is None:
                missing.append(x_label)
            if history_y is None:
                missing.append(y_label)
            print_info(f"Could not find history columns for: {', '.join(missing)}")
    
    # Create video
    print_subheader("Generating Video")
    create_video(
        x_data, y_data, color_data,
        model_data, iter_data,
        x_label, y_label, color_label,
        output_path,
        duration=duration,
        fps=fps,
        flip_y=flip_y,
        sigma_clip=sigma_clip,
        history_x=history_x,
        history_y=history_y,
    )
    
    print(f"\n{GREEN}Done!{RESET}")


def run_batch(
    filepath: str,
    x_col: str,
    y_col: str,
    color_col: str,
    output: Optional[str] = None,
    duration: float = 30.0,
    fps: int = 30,
    cmap: str = "viridis",
    flip_y: bool = False,
    sigma_clip: float = 3.0,
    history_file: str = "../LOGS/history.data",
) -> None:
    """Run in batch mode."""
    column_names, data = load_iteration_data(filepath)
    
    if 'model' not in column_names or 'iter' not in column_names:
        print_error("Data must have 'model' and 'iter' columns")
        sys.exit(1)
    
    data, _ = sort_by_model_and_iter(data, column_names)
    
    md = load_history_data(history_file)
    
    model_idx = column_names.index('model')
    iter_idx = column_names.index('iter')
    model_data = data[:, model_idx]
    iter_data = data[:, iter_idx]
    
    x_data, x_label = resolve_axis(x_col, column_names, data)
    y_data, y_label = resolve_axis(y_col, column_names, data)
    color_data, color_label = resolve_axis(color_col, column_names, data)
    
    if output is None:
        safe_x = re.sub(r'[^\w\-]', '_', x_label)
        safe_y = re.sub(r'[^\w\-]', '_', y_label)
        output = f"newton_iter_{safe_y}_vs_{safe_x}.mp4"
    
    history_x, history_y = None, None
    if md is not None:
        history_x = resolve_history_axis(x_label, md)
        history_y = resolve_history_axis(y_label, md)
    
    create_video(
        x_data, y_data, color_data,
        model_data, iter_data,
        x_label, y_label, color_label,
        output,
        duration=duration,
        fps=fps,
        cmap=cmap,
        flip_y=flip_y,
        sigma_clip=sigma_clip,
        history_x=history_x,
        history_y=history_y,
    )


# ============================================================================
#                              CLI ENTRY
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Create video of MESA Colors per-iteration data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                                    # Interactive mode
  %(prog)s -f SED/iteration_colors.data       # Specify file
  %(prog)s -x Teff -y V -c model              # Batch mode
  %(prog)s -x Teff -y "V-U" -c iter           # With color index expression
  %(prog)s -x Teff -y V -c model --flip-y     # Flip Y axis (for magnitudes)
  %(prog)s -x Teff -y V -c model -d 60        # 60 second video
  %(prog)s -x Teff -y V -c model --sigma 2.5  # Stricter outlier removal
        """
    )
    
    parser.add_argument(
        "-f", "--file",
        default="SED/iteration_colors.data",
        help="Path to iteration colors data file (default: SED/iteration_colors.data)"
    )
    
    parser.add_argument(
        "--history",
        default="../LOGS/history.data",
        help="Path to MESA history file (default: ../LOGS/history.data)"
    )
    
    parser.add_argument("-x", "--x-col", help="X-axis column name or expression")
    parser.add_argument("-y", "--y-col", help="Y-axis column name or expression")
    parser.add_argument("-c", "--color-col", help="Color axis column name or expression")
    parser.add_argument("-o", "--output", help="Output video filename")
    parser.add_argument("-d", "--duration", type=float, default=30.0,
                        help="Target video duration in seconds (default: 30)")
    parser.add_argument("--fps", type=int, default=30, help="Frames per second (default: 30)")
    parser.add_argument("--cmap", default="viridis", help="Matplotlib colormap (default: viridis)")
    parser.add_argument("--flip-y", action="store_true", help="Flip Y axis")
    parser.add_argument("--sigma", type=float, default=3.0,
                        help="Sigma clipping threshold for outlier removal (default: 3.0)")
    parser.add_argument("--list-columns", action="store_true", help="List columns and exit")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.file):
        alt_path = os.path.join("..", args.file)
        if os.path.exists(alt_path):
            args.file = alt_path
        else:
            print_error(f"File not found: {args.file}")
            sys.exit(1)
    
    if args.list_columns:
        column_names, data = load_iteration_data(args.file)
        print_header("Available Columns")
        for i, name in enumerate(column_names):
            print(f"  [{GREEN}{i:2d}{RESET}] {name}")
        print(f"\n{DIM}Total: {len(column_names)} columns, {data.shape[0]} rows{RESET}")
        return
    
    if args.x_col and args.y_col and args.color_col:
        run_batch(
            filepath=args.file,
            x_col=args.x_col,
            y_col=args.y_col,
            color_col=args.color_col,
            output=args.output,
            duration=args.duration,
            fps=args.fps,
            cmap=args.cmap,
            flip_y=args.flip_y,
            sigma_clip=args.sigma,
            history_file=args.history,
        )
    else:
        run_interactive(args.file, args.history)


if __name__ == "__main__":
    main()