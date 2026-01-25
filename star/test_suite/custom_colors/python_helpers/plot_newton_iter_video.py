#!/usr/bin/env python3
"""
plot_newton_iter_video.py — Video plotter for MESA Colors per-iteration output

Creates an animated video showing Newton iterations accumulating over time.
Points appear one by one, sorted by model number then iteration number.

Author: Niall Miller (2025)
"""

import argparse
import os
import re
import shutil
import sys
from typing import List, Optional, Tuple, Union

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter
from matplotlib.colors import Normalize
from matplotlib import cm


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
    
    Args:
        data_arrays: List of 1D arrays to clip (e.g., [x_data, y_data])
        sigma: Number of standard deviations for clipping threshold
        max_iter: Maximum iterations for convergence
    
    Returns:
        Boolean mask where True = KEEP the point, False = outlier
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
            # Use median absolute deviation (more robust than std)
            mad = np.median(np.abs(valid_data - median))
            # Convert MAD to equivalent std (for normal distribution)
            std_equiv = mad * 1.4826
            
            if std_equiv > 0:
                deviation = np.abs(arr - median)
                mask &= (deviation < sigma * std_equiv)
        
        # Check for convergence
        if np.array_equal(mask, prev_mask):
            break
    
    return mask


def percentile_clip_mask(
    data_arrays: List[np.ndarray],
    lower_percentile: float = 1.0,
    upper_percentile: float = 99.0
) -> np.ndarray:
    """
    Create a mask removing points outside percentile bounds.
    
    Args:
        data_arrays: List of 1D arrays to clip
        lower_percentile: Lower percentile cutoff (default 1%)
        upper_percentile: Upper percentile cutoff (default 99%)
    
    Returns:
        Boolean mask where True = KEEP the point
    """
    n_points = len(data_arrays[0])
    mask = np.ones(n_points, dtype=bool)
    
    for arr in data_arrays:
        lower = np.percentile(arr, lower_percentile)
        upper = np.percentile(arr, upper_percentile)
        mask &= (arr >= lower) & (arr <= upper)
    
    return mask


# ============================================================================
#                              TERMINAL UI
# ============================================================================

def use_color() -> bool:
    """Check if we should use ANSI colors."""
    return sys.stdout.isatty() and ("NO_COLOR" not in os.environ)


# ANSI color codes
BOLD  = "\x1b[1m"  if use_color() else ""
DIM   = "\x1b[2m"  if use_color() else ""
CYAN  = "\x1b[36m" if use_color() else ""
YELL  = "\x1b[33m" if use_color() else ""
GREEN = "\x1b[32m" if use_color() else ""
RED   = "\x1b[31m" if use_color() else ""
RESET = "\x1b[0m"  if use_color() else ""


def term_width() -> int:
    """Get terminal width."""
    return shutil.get_terminal_size().columns


def print_header(title: str) -> None:
    """Print a styled header."""
    width = min(70, term_width())
    print(f"\n{BOLD}{CYAN}{'═' * width}{RESET}")
    print(f"{BOLD}{CYAN}  {title}{RESET}")
    print(f"{BOLD}{CYAN}{'═' * width}{RESET}\n")


def print_subheader(title: str) -> None:
    """Print a styled subheader."""
    width = min(70, term_width())
    print(f"\n{BOLD}{title}{RESET}")
    print(f"{DIM}{'─' * width}{RESET}")


def print_success(msg: str) -> None:
    """Print a success message."""
    print(f"{GREEN}✓{RESET} {msg}")


def print_error(msg: str) -> None:
    """Print an error message."""
    print(f"{RED}✗{RESET} {msg}")


def print_info(msg: str) -> None:
    """Print an info message."""
    print(f"{CYAN}ℹ{RESET} {msg}")


# ============================================================================
#                           DATA LOADING
# ============================================================================

def load_iteration_data(filepath: str) -> Tuple[List[str], np.ndarray]:
    """Load iteration colors data file."""
    
    # Read header to get column names
    with open(filepath, 'r') as f:
        header_line = f.readline().strip()
    
    # Parse column names (handle the # at the start)
    if header_line.startswith('#'):
        header_line = header_line[1:].strip()
    
    column_names = header_line.split()
    
    # Load the data
    data = np.loadtxt(filepath, comments='#')
    
    # Handle single-row data
    if data.ndim == 1:
        data = data.reshape(1, -1)
    
    return column_names, data


def sort_by_model_and_iter(
    data: np.ndarray, 
    column_names: List[str]
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Sort data by model number, then by iteration number.
    
    Returns:
        (sorted_data, sort_indices)
    """
    model_idx = column_names.index('model')
    iter_idx = column_names.index('iter')
    
    # Create sort key: (model, iter)
    sort_indices = np.lexsort((data[:, iter_idx], data[:, model_idx]))
    
    return data[sort_indices], sort_indices


def get_final_iteration_mask(model_data: np.ndarray, iter_data: np.ndarray) -> np.ndarray:
    """
    Identify the final Newton iteration for each timestep (model).
    
    Returns:
        Boolean mask where True indicates the final iteration for that model
    """
    mask = np.zeros(len(model_data), dtype=bool)
    
    unique_models = np.unique(model_data)
    
    for model in unique_models:
        model_mask = model_data == model
        model_indices = np.where(model_mask)[0]
        iters_for_model = iter_data[model_mask]
        max_iter_local_idx = np.argmax(iters_for_model)
        final_idx = model_indices[max_iter_local_idx]
        mask[final_idx] = True
    
    return mask


# ============================================================================
#                        EXPRESSION PARSING
# ============================================================================

def is_expression(s: str) -> bool:
    """Check if string contains math operators (is an expression)."""
    return bool(re.search(r'(?<!^)[\+\-\*/]', s.replace(' ', '')))


def parse_expression(
    expr: str, 
    column_names: List[str], 
    data: np.ndarray
) -> Tuple[np.ndarray, str]:
    """Parse and evaluate a column expression."""
    original_expr = expr
    expr = expr.strip()
    
    namespace = {}
    
    def replace_index(match):
        idx = int(match.group(1))
        if idx < 0 or idx >= len(column_names):
            raise ValueError(f"Column index [{idx}] out of range (0-{len(column_names)-1})")
        var_name = f"__col_{idx}__"
        namespace[var_name] = data[:, idx]
        return var_name
    
    expr = re.sub(r'\[(\d+)\]', replace_index, expr)
    
    sorted_names = sorted(column_names, key=len, reverse=True)
    
    for i, name in enumerate(sorted_names):
        orig_idx = column_names.index(name)
        var_name = f"__col_{orig_idx}__"
        pattern = r'\b' + re.escape(name) + r'\b'
        if re.search(pattern, expr):
            namespace[var_name] = data[:, orig_idx]
            expr = re.sub(pattern, var_name, expr)
    
    safe_funcs = {
        'abs': np.abs,
        'sqrt': np.sqrt,
        'log': np.log,
        'log10': np.log10,
        'exp': np.exp,
        'sin': np.sin,
        'cos': np.cos,
        'tan': np.tan,
        'pi': np.pi,
    }
    namespace.update(safe_funcs)
    
    if not re.match(r'^[\w\s\+\-\*/\(\)\.\,]+$', expr):
        raise ValueError(f"Invalid characters in expression: {original_expr}")
    
    try:
        result = eval(expr, {"__builtins__": {}}, namespace)
    except Exception as e:
        raise ValueError(f"Failed to evaluate expression '{original_expr}': {e}")
    
    if np.isscalar(result):
        result = np.full(data.shape[0], result)
    
    return result, original_expr


def resolve_axis(
    spec: Union[int, str],
    column_names: List[str],
    data: np.ndarray
) -> Tuple[np.ndarray, str]:
    """Resolve an axis specification to data array and label."""
    if isinstance(spec, int):
        if 0 <= spec < len(column_names):
            return data[:, spec], column_names[spec]
        raise ValueError(f"Column index {spec} out of range")
    
    spec = str(spec).strip()
    
    if spec in column_names:
        idx = column_names.index(spec)
        return data[:, idx], spec
    
    match = re.match(r'^\[?(\d+)\]?$', spec)
    if match:
        idx = int(match.group(1))
        if 0 <= idx < len(column_names):
            return data[:, idx], column_names[idx]
        raise ValueError(f"Column index {idx} out of range")
    
    for i, name in enumerate(column_names):
        if name.lower() == spec.lower():
            return data[:, i], name
    
    return parse_expression(spec, column_names, data)


# ============================================================================
#                              VIDEO CREATION
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
    cmap: str = "viridis",
    point_size: int = 20,
    alpha: float = 0.7,
    flip_y: bool = False,
    fps: int = 30,
    dpi: int = 150,
    sigma_clip: float = 3.0,
) -> None:
    """
    Create an animated video of points appearing one by one.
    
    Points appear in order: model 0 iter 0, model 0 iter 1, ..., model 1 iter 0, ...
    Video duration is fixed at `duration` seconds (default 30s) if enough points.
    Axis limits dynamically expand as points are added.
    Outliers are removed via sigma clipping.
    """
    
    # =========================================
    # SIGMA CLIPPING - Remove outliers
    # =========================================
    print_info(f"Applying {sigma_clip}-sigma clipping to remove outliers...")
    
    # Clip based on x and y data
    valid_mask = sigma_clip_mask([x_data, y_data], sigma=sigma_clip)
    n_removed = np.sum(~valid_mask)
    
    if n_removed > 0:
        print_info(f"Removed {n_removed} outlier points ({100*n_removed/len(x_data):.1f}%)")
        
        # Apply mask to all arrays
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
    
    # Calculate how many points per frame to achieve desired duration
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
    
    # Identify final iterations for marking
    final_iter_mask = get_final_iteration_mask(model_data, iter_data)
    
    # =========================================
    # PRECOMPUTE AXIS LIMITS FOR EACH FRAME
    # =========================================
    # This enables smooth dynamic zooming
    print_info("Precomputing dynamic axis limits...")
    
    # Final limits (with padding)
    x_min_final, x_max_final = np.min(x_data), np.max(x_data)
    y_min_final, y_max_final = np.min(y_data), np.max(y_data)
    
    x_range_final = x_max_final - x_min_final
    y_range_final = y_max_final - y_min_final
    
    # Add padding
    pad_frac = 0.08
    x_min_final -= x_range_final * pad_frac
    x_max_final += x_range_final * pad_frac
    y_min_final -= y_range_final * pad_frac
    y_max_final += y_range_final * pad_frac
    
    # Set up the figure
    fig, ax = plt.subplots(figsize=(10, 8))
    
    ax.set_xlabel(x_label, fontsize=12)
    ax.set_ylabel(y_label, fontsize=12)
    ax.grid(True, alpha=0.3)
    
    # Color normalization (use clipped data range)
    norm = Normalize(vmin=np.min(color_data), vmax=np.max(color_data))
    colormap = cm.get_cmap(cmap)
    
    # Add colorbar
    sm = cm.ScalarMappable(cmap=colormap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label(color_label, fontsize=12)
    
    # Title with counters (will be updated)
    title = ax.set_title("", fontsize=14, fontweight='bold')
    
    # Initialize empty scatter plots
    scatter = ax.scatter([], [], c=[], cmap=cmap, norm=norm, s=point_size, alpha=alpha)
    scatter_final = ax.scatter([], [], c='black', marker='x', s=point_size*3, 
                                linewidths=2, label='Final iteration', zorder=10)
    ax.legend(loc='upper right')
    
    plt.tight_layout()
    
    def get_dynamic_limits(n_show: int) -> Tuple[float, float, float, float]:
        """
        Compute axis limits that smoothly expand to contain all visible points.
        Starts zoomed in on first points, ends at final limits.
        """
        x_visible = x_data[:n_show]
        y_visible = y_data[:n_show]
        
        # Current data bounds
        x_min_curr = np.min(x_visible)
        x_max_curr = np.max(x_visible)
        y_min_curr = np.min(y_visible)
        y_max_curr = np.max(y_visible)
        
        # Add padding to current bounds
        x_range_curr = max(x_max_curr - x_min_curr, x_range_final * 0.01)
        y_range_curr = max(y_max_curr - y_min_curr, y_range_final * 0.01)
        
        x_min_padded = x_min_curr - x_range_curr * pad_frac
        x_max_padded = x_max_curr + x_range_curr * pad_frac
        y_min_padded = y_min_curr - y_range_curr * pad_frac
        y_max_padded = y_max_curr + y_range_curr * pad_frac
        
        # Smoothly interpolate towards final limits
        # Progress factor: how far through the video we are
        progress = n_show / n_points
        
        # Use smooth easing (ease-out cubic)
        ease = 1 - (1 - progress) ** 2
        
        # Blend current bounds with final bounds
        # Early: mostly current bounds (zoomed in)
        # Late: mostly final bounds (zoomed out)
        x_min = x_min_padded + (x_min_final - x_min_padded) * ease
        x_max = x_max_padded + (x_max_final - x_max_padded) * ease
        y_min = y_min_padded + (y_min_final - y_min_padded) * ease
        y_max = y_max_padded + (y_max_final - y_max_padded) * ease
        
        return x_min, x_max, y_min, y_max
    
    def init():
        scatter.set_offsets(np.empty((0, 2)))
        scatter.set_array(np.array([]))
        scatter_final.set_offsets(np.empty((0, 2)))
        # Set initial limits
        ax.set_xlim(x_min_final, x_max_final)
        if flip_y:
            ax.set_ylim(y_max_final, y_min_final)
        else:
            ax.set_ylim(y_min_final, y_max_final)
        return scatter, scatter_final, title
    
    def update(frame):
        # Calculate how many points to show
        n_show = min((frame + 1) * points_per_frame, n_points)
        
        # Get data up to this point
        x_show = x_data[:n_show]
        y_show = y_data[:n_show]
        c_show = color_data[:n_show]
        
        # Update main scatter
        scatter.set_offsets(np.column_stack([x_show, y_show]))
        scatter.set_array(c_show)
        
        # Update final iteration markers
        final_mask_show = final_iter_mask[:n_show]
        if np.any(final_mask_show):
            scatter_final.set_offsets(np.column_stack([
                x_show[final_mask_show], 
                y_show[final_mask_show]
            ]))
        else:
            scatter_final.set_offsets(np.empty((0, 2)))
        
        # Update dynamic axis limits
        x_min, x_max, y_min, y_max = get_dynamic_limits(n_show)
        ax.set_xlim(x_min, x_max)
        if flip_y:
            ax.set_ylim(y_max, y_min)
        else:
            ax.set_ylim(y_min, y_max)
        
        # Update title with current model/iter info
        current_model = int(model_data[n_show - 1])
        current_iter = int(iter_data[n_show - 1])
        n_final = np.sum(final_mask_show)
        title.set_text(f"Newton Iteration Colors — Model {current_model}, Iter {current_iter}\n"
                      f"Points: {n_show}/{n_points} | Timesteps: {n_final}")
        
        return scatter, scatter_final, title
    
    # Create animation - use blit=False since we're changing axis limits
    print_info("Generating animation frames...")
    anim = FuncAnimation(
        fig, update, frames=total_frames,
        init_func=init, blit=False, interval=1000/fps
    )
    
    # Save video
    print_info(f"Encoding video to {output_path}...")
    
    # Try FFmpeg first, fall back to Pillow for GIF
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
#                           INTERACTIVE MODE
# ============================================================================

def prompt_axis_choice(
    column_names: List[str], 
    data: np.ndarray,
    label: str
) -> Optional[Tuple[np.ndarray, str]]:
    """
    Prompt for axis selection - can be column index, name, or expression.
    
    Returns:
        (data_array, label) or None for quit
    """
    N = len(column_names)
    
    print_subheader(f"{label} ({CYAN}{N}{RESET} columns)")
    
    # Print in columns
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
        
        # Try to resolve as axis (handles indices, names, and expressions)
        try:
            arr, lbl = resolve_axis(inp, column_names, data)
            return arr, lbl
        except ValueError as e:
            print_error(str(e))


def run_interactive(filepath: str) -> None:
    """Run interactive mode."""
    print_header("MESA Colors — Newton Iteration Video Maker")
    
    # Load data
    print_info(f"Loading: {filepath}")
    try:
        column_names, data = load_iteration_data(filepath)
    except Exception as e:
        print_error(f"Failed to load file: {e}")
        sys.exit(1)
    
    # Verify required columns exist
    if 'model' not in column_names or 'iter' not in column_names:
        print_error("Data must have 'model' and 'iter' columns")
        sys.exit(1)
    
    # Sort by model then iter
    data, _ = sort_by_model_and_iter(data, column_names)
    print_success(f"Loaded {data.shape[0]} data points, {data.shape[1]} columns")
    print_success("Data sorted by model number, then iteration")
    
    # Get model and iter columns
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
    
    # Flip Y?
    flip_inp = input(f"Flip Y axis? {DIM}[y/N]{RESET} ").strip().lower()
    flip_y = flip_inp in ('y', 'yes')
    
    # Select color axis
    print(f"\n{DIM}Select Color axis:{RESET}")
    result = prompt_axis_choice(column_names, data, "Color")
    if result is None:
        return
    color_data, color_label = result
    print_success(f"Color: {color_label}")
    
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
    
    # Create video
    print_subheader("Generating Video")
    create_video(
        x_data, y_data, color_data,
        model_data, iter_data,
        x_label, y_label, color_label,
        output_path,
        flip_y=flip_y,
        sigma_clip=sigma_clip
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
) -> None:
    """Run in batch mode."""
    # Load data
    column_names, data = load_iteration_data(filepath)
    
    # Verify required columns
    if 'model' not in column_names or 'iter' not in column_names:
        print_error("Data must have 'model' and 'iter' columns")
        sys.exit(1)
    
    # Sort by model then iter
    data, _ = sort_by_model_and_iter(data, column_names)
    
    # Get model and iter columns
    model_idx = column_names.index('model')
    iter_idx = column_names.index('iter')
    model_data = data[:, model_idx]
    iter_data = data[:, iter_idx]
    
    # Resolve axes
    x_data, x_label = resolve_axis(x_col, column_names, data)
    y_data, y_label = resolve_axis(y_col, column_names, data)
    color_data, color_label = resolve_axis(color_col, column_names, data)
    
    # Output filename
    if output is None:
        safe_x = re.sub(r'[^\w\-]', '_', x_label)
        safe_y = re.sub(r'[^\w\-]', '_', y_label)
        output = f"newton_iter_{safe_y}_vs_{safe_x}.mp4"
    
    # Create video
    create_video(
        x_data, y_data, color_data,
        model_data, iter_data,
        x_label, y_label, color_label,
        output,
        duration=duration,
        fps=fps,
        cmap=cmap,
        flip_y=flip_y,
        sigma_clip=sigma_clip
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
    
    # Check file exists
    if not os.path.exists(args.file):
        alt_path = os.path.join("..", args.file)
        if os.path.exists(alt_path):
            args.file = alt_path
        else:
            print_error(f"File not found: {args.file}")
            sys.exit(1)
    
    # List columns mode
    if args.list_columns:
        column_names, data = load_iteration_data(args.file)
        print_header("Available Columns")
        for i, name in enumerate(column_names):
            print(f"  [{GREEN}{i:2d}{RESET}] {name}")
        print(f"\n{DIM}Total: {len(column_names)} columns, {data.shape[0]} rows{RESET}")
        return
    
    # Batch or interactive
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
        )
    else:
        run_interactive(args.file)


if __name__ == "__main__":
    main()