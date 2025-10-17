import os
import re
import math

def eyring_equation(barrier, T=423.15):
    """Computes the rate constant using the Eyring equation with correct unit handling."""
    kB = 1.380649e-23  # J/K (Boltzmann constant)
    h = 6.62607015e-34  # J·s (Planck constant)
    R = 8.314  # J/(mol·K) (Gas constant)

    # convert J/mol into kJ/mol
    barrier_kJ = barrier/1000
    # Compute rate constant
    rate_constant = ((kB * T)/h) * math.exp(-(barrier_kJ / (R*T)))
    
    return rate_constant

def extract_rate_constants(mkm_file):
    """Extracts transition state energies (in kJ/mol), and calculates rate constants."""
    print(f"🔍 Reading TS energy data from: {mkm_file}")
    rate_constants = {}

    with open(mkm_file, "r") as f:
        lines = f.readlines()

    for line in lines:
        line = line.strip()

        # Match transition state energy lines (now in kJ/mol)
        match = re.search(r";\s*(forw_\w+)\s*;\s*(back_\w+)\s*;\s*([\d\.\-e+]+)\s*;\s*([\d\.\-e+]+)\s*;", line)
        if match:
            forw_name = match.group(1)
            back_name = match.group(2)
            forw_energy_kJmol = float(match.group(3))  # Already in kJ/mol
            back_energy_kJmol = float(match.group(4))  # Already in kJ/mol

            # Compute rate constants using fixed Eyring equation
            rate_constants[forw_name] = eyring_equation(forw_energy_kJmol)
            rate_constants[back_name] = eyring_equation(back_energy_kJmol)

    if not rate_constants:
        print(f"⚠️ No TS energies found in {mkm_file}. Check the file format.")

    return rate_constants

def update_rate_constants(mkm_file, rate_constants, output_file):
    """Replaces placeholders with computed rate constants in the input.mkm file."""
    print(f"✏️ Updating rate constants in: {mkm_file}")

    with open(mkm_file, "r") as f:
        lines = f.readlines()

    updated_lines = []
    replacements_made = 0

    for line in lines:
        original_line = line  # Store the original line for debugging
        for key, rate in rate_constants.items():
            if key in line:
                line = line.replace(key, f"{rate:.6e}")
                replacements_made += 1
        
        updated_lines.append(line)

        if original_line != line:
            print(f"🔄 BEFORE: {original_line.strip()}\n   AFTER : {line.strip()}\n")

    with open(output_file, "w") as f:
        f.writelines(updated_lines)

    print(f"✅ Updated {output_file} with {replacements_made} replacements.")

    if replacements_made == 0:
        print(f"⚠️ No placeholders were replaced in {mkm_file}. Check if the rate constants match the placeholders.")

def main():
    input_folder = "input"
    mechanisms = ["lh", "l2", "mix", "mix2", "mvk"]

    # Create output directories
    output_dirs = {mech: f"Kinetics_{mech.capitalize()}" for mech in mechanisms}
    for dir_path in output_dirs.values():
        os.makedirs(dir_path, exist_ok=True)

    for mech in mechanisms:
        input_mkm = os.path.join(input_folder, f"{mech}_input.mkm")
        
        if not os.path.exists(input_mkm):
            print(f"⚠️ Warning: {input_mkm} not found. Skipping...")
            continue

        # Extract rate constants from {mech}_input.mkm
        rate_constants = extract_rate_constants(input_mkm)

        # Update the file with rate constants
        updated_mkm = os.path.join(output_dirs[mech], f"{mech}_input.mkm")
        update_rate_constants(input_mkm, rate_constants, updated_mkm)

        # Copy {mech}_input.mkm as input.mkm inside the folder
        final_mkm = os.path.join(output_dirs[mech], "input.mkm")
        with open(updated_mkm, "r") as src, open(final_mkm, "w") as dest:
            dest.writelines(src.readlines())
        
        # Remove {mech}_input.mkm after copying
        os.remove(updated_mkm)
        print(f"🗑️ Removed {updated_mkm} after copying to {final_mkm}")

if __name__ == "__main__":
    main()


