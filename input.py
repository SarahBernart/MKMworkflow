import os
import re

def normalize_species_name(species):
    """Removes extra spaces and normalizes species names for comparison."""
    return " ".join(species.split())

def parse_ped_file(ped_file):
    """Parses the ped_0.txt file to extract species energies and assign indices for each mechanism."""
    print(f"Reading energy data from: {ped_file}")
    energy_lookup = {"lh": {}, "l2": {}, "mix": {},"mix2": {}, "mvk": {}}
    index_lookup = {"lh": {}, "l2": {}, "mix": {},"mix2": {}, "mvk": {}}
    rc_values = {"lh": [], "l2": [], "mix": [],"mix2": [], "mvk": []}
    
    with open(ped_file, "r") as f:
        lines = f.readlines()
    
    current_mech = None
    current_species = None
    counter = {"lh": 0, "l2": 0, "mix": 0, "mix2": 0, "mvk": 0}
    
    for line in lines:
        line = line.strip()
        match = re.match(r"^(lh|l2|mix|mix2|mvk)\s*\|\s*(.+)", line)
        if match:
            current_mech = match.group(1)
            current_species = normalize_species_name(match.group(2).strip())
            index_lookup[current_mech][counter[current_mech]] = current_species  # Map index to species name
            rc_values[current_mech].append(counter[current_mech])
            counter[current_mech] += 1
            continue
        
        energy_match = re.search(r"=\s*(-?\d+\.\d+)\s*eV", line)
        if energy_match and current_mech is not None and current_species is not None:
            energy_lookup[current_mech][counter[current_mech] - 1] = float(energy_match.group(1))
    
    return energy_lookup, index_lookup, rc_values

def extract_reaction_energies(mech_file, energy_lookup):
    """Extracts and calculates reaction energies from the kinetics file using numerical indices."""
    print(f"Reading reaction data from: {mech_file}")
    reaction_energies = {}
    conversion_factor = 96.48533212  # Convert eV to kJ/mol
    
    with open(mech_file, "r") as f:
        lines = f.readlines()
    
    for line in lines:
        line = line.strip()
        match = re.match(r"(\w+):\s*\{(\d+)\}\s*-\s*\{(\d+)\}", line)
        if match:
            reaction_name = match.group(1).strip()
            index_A = int(match.group(2))
            index_B = int(match.group(3))
            
            energy_A = energy_lookup.get(index_A, None)
            energy_B = energy_lookup.get(index_B, None)
            
            if energy_A is not None and energy_B is not None:
                reaction_energies[reaction_name] = (energy_A - energy_B) * conversion_factor
            else:
                print(f"⚠️ Missing energy values for {reaction_name}: {index_A}({energy_A}) - {index_B}({energy_B})")
    
    print(f"Extracted reaction energies: {reaction_energies}")
    return reaction_energies

def update_placeholder_file(placeholder_file, reaction_energies, output_file):
    """Replaces placeholders in input_placeholder.mkm with computed reaction energies."""
    print(f"Reading placeholder file from: {placeholder_file}")
    
    with open(placeholder_file, "r") as f:
        lines = f.readlines()
    
    updated_lines = []
    replacements_made = 0
    
    for line in lines:
        original_line = line  # Keep original for debugging
        for reaction, energy in reaction_energies.items():
            if reaction in line:
                line = line.replace(reaction, f"{energy:.0f}e3")
                replacements_made += 1
        updated_lines.append(line)
        
        # Debugging: Print before & after if a change was made
        if original_line != line:
            print(f"BEFORE: {original_line.strip()}\nAFTER : {line.strip()}\n")
    
    print(f"Total replacements made: {replacements_made}")
    
    print(f"Writing updated file to: {output_file}")
    with open(output_file, "w") as f:
        f.writelines(updated_lines)

def main():
    ped_file = "temp/ped_0.txt"
    energy_lookup, index_lookup, rc_values = parse_ped_file(ped_file)
    
    mechanisms = ["lh", "l2", "mix", "mix2", "mvk"]
    kinetics_files = {mech: f"mech/{mech}_mech_kinetics.txt" for mech in mechanisms}
    placeholder_files = {mech: f"input/input_placeholder_{mech}.mkm" for mech in mechanisms}
    output_folder = "input"  # ✅ All output files will be saved in input/ 

    # Ensure input folder exists
    os.makedirs(output_folder, exist_ok=True)

    for mech in mechanisms:
        reaction_energies = extract_reaction_energies(kinetics_files[mech], energy_lookup[mech])

        # ✅ Write the updated file into the input folder without overwriting input.mkm
        output_mkm = os.path.join(output_folder, f"{mech}_input.mkm")
        update_placeholder_file(placeholder_files[mech], reaction_energies, output_mkm)

        print(f"✅ Updated {output_mkm}")

if __name__ == "__main__":
    main()
