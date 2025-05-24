import os
import math

def calculate_entropy(frequencies, temperature=423.15):
    """
    Calculate vibrational entropy from frequencies in eV.
    """
    k_B = 8.617333262e-5  # Boltzmann constant in eV/K
    h = 4.135667696e-15  # Planck's constant in eV*s
    c = 2.99792e10  # Speed of light in cm/s
    entropy = 0.0

    for freq in frequencies:
        freq_eV = h * c * freq  # Convert cm⁻¹ to eV
        x = freq_eV / (k_B * temperature)

        if x > 1e-6:  # Avoid numerical instability for very small x
            entropy += (x / (math.exp(x) - 1)) - math.log1p(-math.exp(-x))

    entropy *= k_B  # Convert to eV/K
    return f"{entropy:.3e}"

def extract_total_energy(selected_folders, gasphase_folders, output_file="data.txt"):
    """
    Extracts:
    - total energies from 'total.e'
    - zero-point energies (ZPE) from multiple 'zpe_*' folders
    - entropy values (S_*) from vibrational frequencies in 'zpe_*' folders
    """
    results = {}  # Total energies
    zpe_results = {}  # Zero-point energies for each zpe_*
    entropy_results = {}  # Entropies for each zpe_*
    
    data_dir = "data"
    os.makedirs(data_dir, exist_ok=True)  # Ensure the 'data' directory exists
    output_path = os.path.join(data_dir, output_file)

    all_folders = selected_folders + gasphase_folders
    found_any_energy = False  # Track if any total.e is found

    for folder in all_folders:
        if not os.path.isdir(folder):
            continue  # Skip if it's not a valid directory

        # ✅ Handling TS_* folders: Look inside _TS subfolders
        if folder.startswith("TS_"):
            print(f"🔍 Searching TS folder: {folder}")

            for subdir in os.listdir(folder):
                subpath = os.path.join(folder, subdir)
                if os.path.isdir(subpath) and subdir.endswith("_TS"):
                    print(f"  🔍 Found _TS subfolder: {subpath}")

                    # Extract total.e from the _TS subfolder
                    total_e_path = os.path.join(subpath, "total.e")
                    if os.path.isfile(total_e_path):
                        try:
                            with open(total_e_path, "r") as f:
                                total_energy = round(float(f.readline().strip()), 3)
                            results[subpath] = total_energy
                            found_any_energy = True
                            print(f"  ✅ Extracted energy from {total_e_path}: {total_energy}")
                        except Exception as e:
                            print(f"  ❌ Error reading {total_e_path}: {e}")

                    # Look for multiple zpe_* folders inside _TS subfolder
                    for zpe_subdir in os.listdir(subpath):
                        zpe_path = os.path.join(subpath, zpe_subdir)
                        if os.path.isdir(zpe_path) and zpe_subdir.startswith("zpe_"):
                            print(f"    🔍 Found ZPE folder: {zpe_path}")
                            zpe_file = os.path.join(zpe_path, "frequencies.txt")
                            
                            if os.path.isfile(zpe_file):
                                try:
                                    frequencies = []
                                    with open(zpe_file, "r") as f:
                                        for line in f:
                                            if "Zero-point energy:" in line:
                                                parts = line.split()
                                                for part in parts:
                                                    try:
                                                        zpe_value = float(part)
                                                        zpe_results.setdefault(subpath, {})[zpe_subdir] = zpe_value
                                                        print(f"    ✅ Extracted ZPE from {zpe_file}: {zpe_value}")
                                                        break
                                                    except ValueError:
                                                        continue
                                            elif len(line.split()) == 3 and line.split()[0].isdigit():
                                                try:
                                                    freq_cm = float(line.split()[2])
                                                    frequencies.append(freq_cm)
                                                except ValueError:
                                                    continue
                                    if frequencies:
                                        entropy_results.setdefault(subpath, {})[zpe_subdir] = calculate_entropy(frequencies)
                                        print(f"    ✅ Calculated Entropy for {zpe_file}: {entropy_results[subpath][zpe_subdir]}")
                                except Exception as e:
                                    print(f"    ❌ Error reading {zpe_file}: {e}")

        # ✅ Handling normal (non-TS) folders
        else:
            total_e_path = os.path.join(folder, "total.e")
            if os.path.isfile(total_e_path):
                try:
                    with open(total_e_path, "r") as f:
                        total_energy = round(float(f.readline().strip()), 3)
                    results[folder] = total_energy
                    found_any_energy = True
                    print(f"✅ Extracted energy from {total_e_path}: {total_energy}")
                except Exception as e:
                    print(f"❌ Error reading {total_e_path}: {e}")

            # Extract multiple ZPE and entropy values for each zpe_* folder
            for zpe_subdir in os.listdir(folder):
                zpe_path = os.path.join(folder, zpe_subdir)
                if os.path.isdir(zpe_path) and zpe_subdir.startswith("zpe_"):
                    print(f"🔍 Found ZPE folder: {zpe_path}")
                    zpe_file = os.path.join(zpe_path, "frequencies.txt")
                    
                    if os.path.isfile(zpe_file):
                        try:
                            frequencies = []
                            with open(zpe_file, "r") as f:
                                for line in f:
                                    if "Zero-point energy:" in line:
                                        parts = line.split()
                                        for part in parts:
                                            try:
                                                zpe_value = float(part)
                                                zpe_results.setdefault(folder, {})[zpe_subdir] = zpe_value
                                                print(f"✅ Extracted ZPE from {zpe_file}: {zpe_value}")
                                                break
                                            except ValueError:
                                                continue
                                    elif len(line.split()) == 3 and line.split()[0].isdigit():
                                        try:
                                            freq_cm = float(line.split()[2])
                                            frequencies.append(freq_cm)
                                        except ValueError:
                                            continue
                            if frequencies:
                                entropy_results.setdefault(folder, {})[zpe_subdir] = calculate_entropy(frequencies)
                                print(f"✅ Calculated Entropy for {zpe_file}: {entropy_results[folder][zpe_subdir]}")
                        except Exception as e:
                            print(f"❌ Error reading {zpe_file}: {e}")

    # Debug: If no energy files were found
    if not found_any_energy:
        print("🚨 No 'total.e' files found! Check if the folders contain energy files.")

    # ✅ WRITE OUTPUT FILE WITHOUT WRITING *_TS SUBFOLDERS
    all_zpe_keys = sorted(set(k for v in zpe_results.values() for k in v.keys()))
    
    with open(output_path, "w") as out_file:
        # Define column headers dynamically
        zpe_headers = "\t".join(all_zpe_keys)
        entropy_headers = "\t".join([f"S_{zpe}" for zpe in all_zpe_keys])
        out_file.write(f"Species\tTotal_E\t{zpe_headers}\t{entropy_headers}\n")
    
        summarized_results = {}
    
        for folder in results.keys():
            # ✅ If this is a TS_* subfolder (e.g., TS_COox_LH/TS1), map it to the main TS_* folder
            if "/" in folder and folder.split("/")[0].startswith("TS_"):
                main_folder = folder.split("/")[0]  # Extract only the TS_* main folder
            else:
                main_folder = folder  # Keep normal species as they are
    
            # ✅ Combine energies, ZPE, and entropy values under the main folder
            if main_folder not in summarized_results:
                summarized_results[main_folder] = {"total_e": results[folder], "zpe": {}, "S": {}}

            for zpe in all_zpe_keys:
                summarized_results[main_folder]["zpe"][zpe] = zpe_results.get(folder, {}).get(zpe, "-")
                summarized_results[main_folder]["S"][zpe] = entropy_results.get(folder, {}).get(zpe, "-")

        # ✅ Write the summarized data (without *_TS subfolders)
        for folder, values in summarized_results.items():
            total_e = values["total_e"]
            zpe_values = {zpe: values["zpe"].get(zpe, "-") for zpe in all_zpe_keys}
            entropy_values = {zpe: values["S"].get(zpe, "-") for zpe in all_zpe_keys}
    
            zpe_str = "\t".join(str(zpe_values[zpe]) for zpe in all_zpe_keys)
            entropy_str = "\t".join(str(entropy_values[zpe]) for zpe in all_zpe_keys)

            out_file.write(f"{folder}\t{total_e}\t{zpe_str}\t{entropy_str}\n")

    print(f"✅ Data extracted and saved to {output_path} (without *_TS subfolders).")

def read_species_file(filepath="mech/species.txt"):
    """
    Reads species.txt and extracts gasphase, mixed, lh, and mvk folder lists.
    """
    species_data = {"gasphase_folders": [], "mixed_folders": [],"mix2_folders": [], "lh_folders": [], "l2_folders": [], "mvk_folders": []}

    try:
        with open(filepath, "r") as file:
            current_key = None
            for line in file:
                line = line.strip()
                if not line or line.startswith("#"):  # Ignore empty lines & comments
                    continue

                if line in species_data:  # Identify section headers
                    current_key = line
                elif current_key:  # Append species to the current section
                    species_data[current_key].append(line)

    except FileNotFoundError:
        print(f"❌ Error: File {filepath} not found.")

    return species_data["gasphase_folders"], species_data["mixed_folders"], species_data["mix2_folders"], species_data["lh_folders"], species_data["l2_folders"], species_data["mvk_folders"]

# ✅ Read species from file
gasphase_folders, mixed_folders, mix2_folders, lh_folders, l2_folders, mvk_folders = read_species_file()

# ✅ Run the extraction process
if __name__ == "__main__":
    print("Processing Mixed Data:")
    extract_total_energy(mixed_folders, gasphase_folders, output_file="mix_data.txt")

    print("Processing Mix2 Data:")
    extract_total_energy(mix2_folders, gasphase_folders, output_file="mix2_data.txt")

    print("Processing LH Data:")
    extract_total_energy(lh_folders, gasphase_folders, output_file="lh_data.txt")

    print("Processing L2 Data:")
    extract_total_energy(l2_folders, gasphase_folders, output_file="l2_data.txt")
 
    print("Processing MvK Data:")
    extract_total_energy(mvk_folders, gasphase_folders, output_file="mvk_data.txt")
