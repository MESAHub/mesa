#!/usr/bin/env python3
"""
plot_newton_iter.py — Interactive plotter for MESA Colors per-iteration output

Plots stellar photometry data from Newton solver iterations with an
interactive column picker and colored terminal UI.

Author: Niall Miller (2025)
"""

import argparse
import math
import os
import re
import shutil
import sys
from typing import List, Optional, Sequence, Tuple, Union

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm


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
#                           COLUMN PICKER
# ============================================================================

def prompt_choice(
    options: Sequence[str],
    label: str,
    allow_back: bool = False,
    max_cols: int = 3,
    filter_enabled: bool = True,
) -> Optional[int]:
    """
    Interactive column picker with grid display and filtering.
    
    Returns:
        Selected index (0-based), None for quit, -1 for back
    """
    if not options:
        print(f"No {label} options available.")
        return None

    labels = list(options)
    N = len(labels)
    max_label_len = max(len(s) for s in labels) + 2
    filt: Optional[Tuple[str, str]] = None  # (kind, pattern)

    def apply_filter(indices: List[int]) -> List[int]:
        if filt is None:
            return indices
        kind, patt = filt
        p = patt.lower()
        if kind == "substr":
            return [i for i in indices if p in labels[i].lower()]
        elif kind == "neg":
            return [i for i in indices if p not in labels[i].lower()]
        else:  # regex
            rx = re.compile(patt, re.I)
            return [i for i in indices if rx.search(labels[i])]

    def highlight(s: str) -> str:
        if not use_color() or filt is None:
            return s
        kind, patt = filt
        if kind != "substr" or not patt:
            return s
        rx = re.compile(re.escape(patt), re.I)
        return rx.sub(lambda m: f"{YELL}{m.group(0)}{RESET}", s)

    def grid_print(visible_ids: List[int]) -> None:
        width = max(70, term_width())
        col_w = 8 + max_label_len
        cols = max(1, min(max_cols, width // col_w))
        
        cells = []
        for i in visible_ids:
            cell = f"[{GREEN}{i:2d}{RESET}] {highlight(labels[i])}"
            cells.append(cell)
        
        # Pad to fill grid
        while len(cells) % cols:
            cells.append("")
        
        rows = [cells[k:k+cols] for k in range(0, len(cells), cols)]
        
        print_subheader(f"{label} ({CYAN}{N}{RESET} columns)")
        for row in rows:
            print("  " + "".join(cell.ljust(col_w) for cell in row))

    all_idx = list(range(N))
    
    while True:
        kept = apply_filter(all_idx)
        grid_print(kept)
        
        # Show filter status
        if filt:
            kind, patt = filt
            if kind == "substr":
                print(f"\n{DIM}Filter: /{patt}{RESET}")
            elif kind == "neg":
                print(f"\n{DIM}Filter: !{patt}{RESET}")
            else:
                print(f"\n{DIM}Filter: //{patt}{RESET}")
        
        # Show controls
        controls = f"\n{DIM}Enter column number"
        if filter_enabled:
            controls += " | /text !text //regex | clear"
        if allow_back:
            controls += " | b=back"
        controls += f" | q=quit{RESET}"
        print(controls)
        
        inp = input(f"{CYAN}>{RESET} ").strip()
        
        if not inp:
            continue
        
        # Quit
        if inp.lower() == "q":
            return None
        
        # Back
        if inp.lower() == "b" and allow_back:
            return -1
        
        # Clear filter
        if inp.lower() == "clear":
            filt = None
            continue
        
        # Substring filter: /text
        if inp.startswith("/") and not inp.startswith("//"):
            filt = ("substr", inp[1:])
            continue
        
        # Negative filter: !text
        if inp.startswith("!"):
            filt = ("neg", inp[1:])
            continue
        
        # Regex filter: //pattern
        if inp.startswith("//"):
            try:
                re.compile(inp[2:])
                filt = ("regex", inp[2:])
            except re.error:
                print_error("Invalid regex pattern")
            continue
        
        # Try to parse as number
        try:
            idx = int(inp)
            if 0 <= idx < N:
                return idx
            else:
                print_error(f"Index must be between 0 and {N-1}")
        except ValueError:
            # Try to match by name
            matches = [i for i, lbl in enumerate(labels) if inp.lower() == lbl.lower()]
            if len(matches) == 1:
                return matches[0]
            elif len(matches) > 1:
                print_error(f"Ambiguous: {len(matches)} columns match '{inp}'")
            else:
                print_error(f"Invalid input: '{inp}'")


def prompt_yes_no(prompt: str, default: bool = True) -> bool:
    """Prompt for yes/no with default."""
    suffix = "[Y/n]" if default else "[y/N]"
    inp = input(f"{prompt} {DIM}{suffix}{RESET} ").strip().lower()
    if not inp:
        return default
    return inp in ("y", "yes")


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


# ============================================================================
#                        EXPRESSION PARSING
# ============================================================================

def is_expression(s: str) -> bool:
    """Check if string contains math operators (is an expression)."""
    # Look for operators not at the start (to allow negative numbers)
    return bool(re.search(r'(?<!^)[\+\-\*/]', s.replace(' ', '')))


def parse_expression(
    expr: str, 
    column_names: List[str], 
    data: np.ndarray
) -> Tuple[np.ndarray, str]:
    """
    Parse and evaluate a column expression.
    
    Supports:
        - Column names: V, Teff, log_g
        - Column indices: [0], [14], [15]
        - Math operators: +, -, *, /
        - Parentheses: (V-U)/(B-V)
        - Constants: 2.5, 1000
    
    Examples:
        "V-U"           -> data[:, V_idx] - data[:, U_idx]
        "[15]-[14]"     -> data[:, 15] - data[:, 14]
        "Teff/1000"     -> data[:, Teff_idx] / 1000
        "2.5*log_g"     -> 2.5 * data[:, log_g_idx]
    
    Returns:
        (result_array, label_string)
    """
    original_expr = expr
    expr = expr.strip()
    
    # Build a safe evaluation namespace
    namespace = {}
    
    # First, replace [N] index references with placeholder variable names
    def replace_index(match):
        idx = int(match.group(1))
        if idx < 0 or idx >= len(column_names):
            raise ValueError(f"Column index [{idx}] out of range (0-{len(column_names)-1})")
        var_name = f"__col_{idx}__"
        namespace[var_name] = data[:, idx]
        return var_name
    
    expr = re.sub(r'\[(\d+)\]', replace_index, expr)
    
    # Sort column names by length (longest first) to avoid partial matches
    sorted_names = sorted(column_names, key=len, reverse=True)
    
    # Replace column names with placeholder variable names
    for i, name in enumerate(sorted_names):
        # Find the original index
        orig_idx = column_names.index(name)
        var_name = f"__col_{orig_idx}__"
        
        # Use word boundaries to avoid partial matches
        # But column names might have special chars, so escape them
        pattern = r'\b' + re.escape(name) + r'\b'
        if re.search(pattern, expr):
            namespace[var_name] = data[:, orig_idx]
            expr = re.sub(pattern, var_name, expr)
    
    # Add safe math functions
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
    
    # Validate expression only contains safe characters
    # Allow: digits, letters, underscores, operators, parentheses, dots, spaces
    if not re.match(r'^[\w\s\+\-\*/\(\)\.\,]+$', expr):
        raise ValueError(f"Invalid characters in expression: {original_expr}")
    
    # Evaluate
    try:
        result = eval(expr, {"__builtins__": {}}, namespace)
    except Exception as e:
        raise ValueError(f"Failed to evaluate expression '{original_expr}': {e}")
    
    # Ensure result is array
    if np.isscalar(result):
        result = np.full(data.shape[0], result)
    
    return result, original_expr


def resolve_axis(
    spec: Union[int, str],
    column_names: List[str],
    data: np.ndarray
) -> Tuple[np.ndarray, str]:
    """
    Resolve an axis specification to data array and label.
    
    Args:
        spec: Column index (int), column name (str), or expression (str)
        column_names: List of column names
        data: Data array
    
    Returns:
        (data_array, label_string)
    """
    # If it's already an integer index
    if isinstance(spec, int):
        if 0 <= spec < len(column_names):
            return data[:, spec], column_names[spec]
        raise ValueError(f"Column index {spec} out of range")
    
    spec = str(spec).strip()
    
    # Check if it's a simple column name
    if spec in column_names:
        idx = column_names.index(spec)
        return data[:, idx], spec
    
    # Check if it's a simple index like "5" or "[5]"
    match = re.match(r'^\[?(\d+)\]?$', spec)
    if match:
        idx = int(match.group(1))
        if 0 <= idx < len(column_names):
            return data[:, idx], column_names[idx]
        raise ValueError(f"Column index {idx} out of range")
    
    # Check case-insensitive column name match
    for i, name in enumerate(column_names):
        if name.lower() == spec.lower():
            return data[:, i], name
    
    # Must be an expression
    return parse_expression(spec, column_names, data)


# ============================================================================
#                              PLOTTING
# ============================================================================




def get_final_iteration_mask(model_data: np.ndarray, iter_data: np.ndarray) -> np.ndarray:
    """
    Identify the final Newton iteration for each timestep (model).
    
    Args:
        model_data: Array of model numbers
        iter_data: Array of iteration numbers
    
    Returns:
        Boolean mask where True indicates the final iteration for that model
    """
    mask = np.zeros(len(model_data), dtype=bool)
    
    # Get unique models
    unique_models = np.unique(model_data)
    
    for model in unique_models:
        # Find indices for this model
        model_mask = model_data == model
        model_indices = np.where(model_mask)[0]
        
        # Find the index with maximum iteration for this model
        iters_for_model = iter_data[model_mask]
        max_iter_local_idx = np.argmax(iters_for_model)
        final_idx = model_indices[max_iter_local_idx]
        
        mask[final_idx] = True
    
    return mask


def create_plot(
    x_data: np.ndarray,
    y_data: np.ndarray,
    color_data: np.ndarray,
    x_label: str,
    y_label: str,
    color_label: str,
    z_data: Optional[np.ndarray] = None,
    z_label: Optional[str] = None,
    cmap: str = "viridis",
    point_size: int = 20,
    alpha: float = 0.7,
    flip_y: bool = False,
    final_iter_mask: Optional[np.ndarray] = None,
) -> plt.Figure:
    """Create the plot with the provided data arrays.
    
    Args:
        final_iter_mask: Boolean mask indicating final iterations to mark with black X
    """
    
    # Create figure
    if z_data is not None:
        # 3D plot
        fig = plt.figure(figsize=(12, 9))
        ax = fig.add_subplot(111, projection='3d')
        
        scatter = ax.scatter(
            x_data, y_data, z_data, 
            c=color_data, cmap=cmap, 
            s=point_size, alpha=alpha
        )
        
        # Plot final iterations with black X markers
        if final_iter_mask is not None and np.any(final_iter_mask):
            ax.scatter(
                x_data[final_iter_mask], 
                y_data[final_iter_mask], 
                z_data[final_iter_mask],
                c='black', marker='x', s=point_size*3, 
                linewidths=2, label='Final iteration', zorder=10
            )
            ax.legend(loc='best')
        
        ax.set_xlabel(x_label, fontsize=12, labelpad=10)
        ax.set_ylabel(y_label, fontsize=12, labelpad=10)
        ax.set_zlabel(z_label, fontsize=12, labelpad=10)
        
        title = f"Newton Iteration Colors\n{z_label} vs {y_label} vs {x_label}"
        
    else:
        # 2D plot
        fig, ax = plt.subplots(figsize=(10, 8))
        
        scatter = ax.scatter(
            x_data, y_data, 
            c=color_data, cmap=cmap, 
            s=point_size, alpha=alpha
        )
        
        # Plot final iterations with black X markers
        if final_iter_mask is not None and np.any(final_iter_mask):
            ax.scatter(
                x_data[final_iter_mask], 
                y_data[final_iter_mask],
                c='black', marker='x', s=point_size*1.1, 
                linewidths=1, label='Final iteration', zorder=10
            )
            ax.legend(loc='best')
        
        ax.set_xlabel(x_label, fontsize=12)
        ax.set_ylabel(y_label, fontsize=12)
        ax.grid(True, alpha=0.3)
        
        title = f"Newton Iteration Colors\n{y_label} vs {x_label}"
    
    if flip_y:
        ax.invert_yaxis()

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax, pad=0.02)
    cbar.set_label(color_label, fontsize=12)
    
    ax.set_title(title, fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    return fig


def save_figure(fig: plt.Figure, base_name: str, formats: List[str], dpi: int = 300) -> None:
    """Save figure in multiple formats."""
    for fmt in formats:
        path = f"{base_name}.{fmt}"
        fig.savefig(path, dpi=dpi, bbox_inches='tight')
        print_success(f"Saved: {path}")


# ============================================================================
#                           MAIN WORKFLOWS
# ============================================================================

def prompt_axis(
    column_names: List[str],
    data: np.ndarray,
    label: str,
    allow_back: bool = True,
) -> Optional[Tuple[np.ndarray, str]]:
    """
    Prompt user for axis selection - can be column or expression.
    
    Returns:
        (data_array, label) or None for quit, or "back" string for back
    """
    while True:
        # Show column picker
        result = prompt_choice(column_names, label, allow_back=allow_back)
        
        if result is None:
            return None
        if result == -1:
            return "back"
        
        # User selected a column by index
        return data[:, result], column_names[result]


def prompt_axis_or_expr(
    column_names: List[str],
    data: np.ndarray,
    label: str,
    allow_back: bool = True,
) -> Optional[Tuple[np.ndarray, str]]:
    """
    Prompt user for axis - can be column selection OR custom expression.
    
    Returns:
        (data_array, label) or None for quit, or "back" string for back
    """
    N = len(column_names)
    max_label_len = max(len(s) for s in column_names) + 2
    
    def grid_print():
        width = max(70, term_width())
        col_w = 8 + max_label_len
        cols = max(1, min(3, width // col_w))
        
        cells = []
        for i, name in enumerate(column_names):
            cell = f"[{GREEN}{i:2d}{RESET}] {name}"
            cells.append(cell)
        
        while len(cells) % cols:
            cells.append("")
        
        rows = [cells[k:k+cols] for k in range(0, len(cells), cols)]
        
        print_subheader(f"{label} ({CYAN}{N}{RESET} columns)")
        for row in rows:
            print("  " + "".join(cell.ljust(col_w) for cell in row))
    
    while True:
        grid_print()
        
        # Show controls with expression hint
        print(f"\n{DIM}Enter: column number | column name | expression (e.g. V-U, [15]-[14], Teff/1000){RESET}")
        controls = f"{DIM}"
        if allow_back:
            controls += "b=back | "
        controls += f"q=quit{RESET}"
        print(controls)
        
        inp = input(f"{CYAN}>{RESET} ").strip()
        
        if not inp:
            continue
        
        # Quit
        if inp.lower() == "q":
            return None
        
        # Back
        if inp.lower() == "b" and allow_back:
            return "back"
        
        # Try to resolve as axis (column or expression)
        try:
            arr, lbl = resolve_axis(inp, column_names, data)
            return arr, lbl
        except ValueError as e:
            print_error(str(e))
            continue


def run_interactive(filepath: str) -> None:
    """Run the interactive plotting workflow."""
    
    print_header("MESA Colors — Newton Iteration Plotter")
    
    # Load data
    print_info(f"Loading: {filepath}")
    try:
        column_names, data = load_iteration_data(filepath)
    except Exception as e:
        print_error(f"Failed to load file: {e}")
        sys.exit(1)
    
    print_success(f"Loaded {data.shape[0]} data points, {data.shape[1]} columns")
    
    # Select plot type
    print_subheader("Plot Type")
    print(f"  [{GREEN}2{RESET}] 2D scatter plot (x, y)")
    print(f"  [{GREEN}3{RESET}] 3D scatter plot (x, y, z)")
    
    while True:
        inp = input(f"\n{CYAN}>{RESET} ").strip()
        if inp in ("2", "3"):
            plot_type = int(inp)
            break
        elif inp.lower() == "q":
            print("Goodbye!")
            return
        print_error("Enter 2 or 3")
    
    # Select X axis
    result = prompt_axis_or_expr(column_names, data, "X-axis", allow_back=False)
    if result is None:
        return
    x_data, x_label = result
    print_success(f"X-axis: {x_label}")
    
    # Select Y axis
    result = prompt_axis_or_expr(column_names, data, "Y-axis", allow_back=True)
    if result is None:
        return
    if result == "back":
        return run_interactive(filepath)
    y_data, y_label = result
    print_success(f"Y-axis: {y_label}")
    
    flip_y = prompt_yes_no("Flip Y axis?", default=False)


    # Select Z axis (for 3D)
    z_data, z_label = None, None
    if plot_type == 3:
        result = prompt_axis_or_expr(column_names, data, "Z-axis", allow_back=True)
        if result is None:
            return
        if result == "back":
            return run_interactive(filepath)
        z_data, z_label = result
        print_success(f"Z-axis: {z_label}")
    
    # Select color axis
    result = prompt_axis_or_expr(column_names, data, "Color axis", allow_back=True)
    if result is None:
        return
    if result == "back":
        return run_interactive(filepath)
    color_data, color_label = result
    print_success(f"Color: {color_label}")
    
    # Generate output filename
    base_name = "newton_iter"
    axes = [x_label, y_label]
    if z_label is not None:
        axes.append(z_label)
    # Sanitize labels for filename
    safe_axes = [re.sub(r'[^\w\-]', '_', a) for a in reversed(axes)]
    base_name += "_" + "_vs_".join(safe_axes)
    
    # Create and save plot
    print_subheader("Generating Plot")
    
    # Compute final iteration mask if model and iter columns exist
    final_iter_mask = None
    if 'model' in column_names and 'iter' in column_names:
        model_idx = column_names.index('model')
        iter_idx = column_names.index('iter')
        final_iter_mask = get_final_iteration_mask(data[:, model_idx], data[:, iter_idx])
        print_info(f"Identified {np.sum(final_iter_mask)} final iterations (marked with black X)")
    
    fig = create_plot(
        x_data, y_data, color_data, 
        x_label, y_label, color_label,
        z_data, z_label, flip_y=flip_y,
        final_iter_mask=final_iter_mask
    )
    save_figure(fig, base_name, ["pdf", "jpg"])
    
    # Show interactive plot
    print_info("Displaying interactive plot (close window to exit)")
    plt.show()
    
    print(f"\n{GREEN}Done!{RESET}")


def run_batch(
    filepath: str,
    x_col: str,
    y_col: str,
    color_col: str,
    z_col: Optional[str] = None,
    output: Optional[str] = None,
    formats: List[str] = ["pdf", "jpg"],
    no_show: bool = False,
    cmap: str = "viridis",
) -> None:
    """Run in batch mode with specified columns or expressions."""
    
    # Load data
    column_names, data = load_iteration_data(filepath)
    
    # Resolve each axis (supports columns and expressions)
    x_data, x_label = resolve_axis(x_col, column_names, data)
    y_data, y_label = resolve_axis(y_col, column_names, data)
    color_data, color_label = resolve_axis(color_col, column_names, data)
    
    z_data, z_label = None, None
    if z_col:
        z_data, z_label = resolve_axis(z_col, column_names, data)
    
    # Generate output name
    if output is None:
        axes = [x_label, y_label]
        if z_label is not None:
            axes.append(z_label)
        # Sanitize labels for filename
        safe_axes = [re.sub(r'[^\w\-]', '_', a) for a in reversed(axes)]
        output = "newton_iter_" + "_vs_".join(safe_axes)
    
    # Compute final iteration mask if model and iter columns exist
    final_iter_mask = None
    if 'model' in column_names and 'iter' in column_names:
        model_idx = column_names.index('model')
        iter_idx = column_names.index('iter')
        final_iter_mask = get_final_iteration_mask(data[:, model_idx], data[:, iter_idx])
    
    # Create plot
    fig = create_plot(
        x_data, y_data, color_data,
        x_label, y_label, color_label,
        z_data, z_label,
        cmap=cmap,
        final_iter_mask=final_iter_mask
    )
    save_figure(fig, output, formats)
    
    if not no_show:
        plt.show()


# ============================================================================
#                              CLI ENTRY
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Plot MESA Colors per-iteration data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                                    # Interactive mode
  %(prog)s -f SED/iteration_colors.data       # Specify file
  %(prog)s -x iter -y Teff -c model           # Batch mode with columns
  %(prog)s -x 1 -y 4 -c 0 --no-show           # Use column indices
  %(prog)s -x iter -y Teff -z R -c model      # 3D plot

Expression examples (color indices, ratios, etc.):
  %(prog)s -x iter -y "V-U" -c model          # Color index V-U
  %(prog)s -x iter -y "[9]-[13]" -c Teff      # Using column indices
  %(prog)s -x "Teff/1000" -y log_g -c model   # Scaled temperature
  %(prog)s -x iter -y "B-V" -c "U-B"          # Color-color diagram

Supported expression syntax:
  - Column names:  V, Teff, log_g, Mag_bol
  - Column indices: [0], [14], [15]
  - Operators: +, -, *, /
  - Parentheses: (V-U)/(B-V)
  - Functions: sqrt(), log10(), abs()
        """
    )
    
    parser.add_argument(
        "-f", "--file",
        default="SED/iteration_colors.data",
        help="Path to iteration colors data file (default: SED/iteration_colors.data)"
    )
    
    parser.add_argument(
        "-x", "--x-col",
        help="X-axis column name or index (enables batch mode)"
    )
    
    parser.add_argument(
        "-y", "--y-col",
        help="Y-axis column name or index"
    )
    
    parser.add_argument(
        "-z", "--z-col",
        help="Z-axis column name or index (for 3D plots)"
    )
    
    parser.add_argument(
        "-c", "--color-col",
        help="Color axis column name or index"
    )
    
    parser.add_argument(
        "-o", "--output",
        help="Output filename base (without extension)"
    )
    
    parser.add_argument(
        "--formats",
        default="pdf,jpg",
        help="Output formats, comma-separated (default: pdf,jpg)"
    )
    
    parser.add_argument(
        "--cmap",
        default="viridis",
        help="Matplotlib colormap (default: viridis)"
    )
    
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Don't display interactive plot"
    )
    
    parser.add_argument(
        "--list-columns",
        action="store_true",
        help="List available columns and exit"
    )
    
    args = parser.parse_args()
    
    # Check file exists
    if not os.path.exists(args.file):
        # Try with ../ prefix
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
    
    # Batch mode if columns specified
    if args.x_col and args.y_col and args.color_col:
        formats = [f.strip() for f in args.formats.split(",")]
        run_batch(
            filepath=args.file,
            x_col=args.x_col,
            y_col=args.y_col,
            color_col=args.color_col,
            z_col=args.z_col,
            output=args.output,
            formats=formats,
            no_show=args.no_show,
            cmap=args.cmap,
        )
    else:
        # Interactive mode
        run_interactive(args.file)


if __name__ == "__main__":
    main()