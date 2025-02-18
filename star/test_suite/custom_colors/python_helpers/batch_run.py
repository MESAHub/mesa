import os
import subprocess
import shutil
import re

# Get the directory where the script is located
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Assume the MESA work directory is one level up from the script
MESA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))  # Moves up one folder

INLIST_PATH = os.path.join(MESA_DIR, "inlist_1.0")  # Path to inlist
HISTORY_FILE = os.path.join(MESA_DIR, "LOGS/history.data")  # Default history file

# List of parameter file names corresponding to inlist options
param_options = [
    "MDwarf.params",
    "PopIII.params",
    "PopII.params"
    "OType.params",
    "BType.params",
    "AType.params",
]

def update_inlist(selected_param, reset=False):
    """Modify the inlist file:
    - If `reset=True`, re-comment all entries.
    - Otherwise, uncomment only the selected parameter file.
    """
    with open(INLIST_PATH, "r") as f:
        inlist_content = f.readlines()

    # Regex pattern to match the commented/uncommented parameter lines
    param_pattern = re.compile(r"^\s*!?\s*(extra_controls_inlist_name\(1\) = '(.+?)')")

    new_content = []
    for line in inlist_content:
        match = param_pattern.match(line)
        if match:
            param_file = match.group(2)
            if reset:
                # Always re-comment every parameter entry
                new_content.append(f"    !{match.group(1)}\n")
            elif param_file == selected_param:
                # Uncomment the selected parameter file
                new_content.append(f"    {match.group(1)}\n")
            else:
                # Keep other entries commented
                new_content.append(f"    !{match.group(1)}\n")
        else:
            new_content.append(line)

    # Write back the modified inlist
    with open(INLIST_PATH, "w") as f:
        f.writelines(new_content)

def run_mesa():
    """Runs the MESA simulation in the correct directory and handles normal termination cases."""
    original_cwd = os.getcwd()  # Save original working directory
    os.chdir(MESA_DIR)  # Change to MESA working directory
    try:
        subprocess.run(["./clean"], check=True)
        subprocess.run(["./mk"], check=True)
        result = subprocess.run(["./rn"], check=True)
        
        # Check MESA output for normal termination cases
        if "termination code: max_age" in result.stdout or "required_termination_code_string" in result.stdout or "1" in result.stdout:
            print("✅ MESA completed successfully (reached termination criterion).")
        else:
            print("❌ Unexpected MESA failure.")
            print(result.stdout)
            print(result.stderr)
            #raise subprocess.CalledProcessError(result.returncode, result.args)
    
    finally:
        os.chdir(original_cwd)  # Restore original directory

def rename_history(param_file):
    """Rename the history file to prevent overwrites."""
    param_name = os.path.splitext(param_file)[0]
    new_history_file = os.path.join(MESA_DIR, f"history_{param_name}.data")
    if os.path.exists(HISTORY_FILE):
        shutil.move(HISTORY_FILE, new_history_file)
        print(f"Saved history as {new_history_file}")
    else:
        print(f"Warning: {HISTORY_FILE} not found after simulation.")

# Re-comment all entries before running batch
update_inlist("", reset=True)

# Run batch simulations
for param_file in param_options:
    print(f"Running MESA with {param_file}...")
    update_inlist(param_file)
    try:
        run_mesa()
        rename_history(param_file)
    except subprocess.CalledProcessError as e:
        print(f"❌ MESA failed for {param_file}. Skipping to next parameter set.")
        print(e)

# Re-comment all parameter files at the end
update_inlist("", reset=True)

print("All simulations completed.")
