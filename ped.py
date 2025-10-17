import os
import numpy as np

# Set temperature (modifiable by searching for "T =")
T = 0  # Temperature in Kelvin

def read_data_file(data_file):
    """Reads species energies, ZPE, and entropy * T values from a data.txt file."""
    data = {}
    with open(data_file, 'r') as f:
        lines = f.readlines()
        headers = lines[0].split()
        for line in lines[1:]:
            values = line.split()
            species = values[0].strip()
            total_e = float(values[1]) if values[1] != '-' else 0.0  # ✅ Fix: Handle missing values
            zpe_values = {headers[i]: float(values[i]) if values[i] != '-' else 0.0 for i in range(2, len(headers)) if headers[i].startswith('zpe_')}
            s_values = {headers[i]: float(values[i]) if values[i] != '-' else 0.0 for i in range(2, len(headers)) if headers[i].startswith('S_')}
            data[species] = {'total_e': total_e, 'zpe': zpe_values, 'S': s_values}
    return data

def process_mechanism(mech_file, data):
    """Processes a mechanism file to compute total energies using only the ZPE/S keys specified per line."""
    mechanism = []
    with open(mech_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            parts = line.strip().split('+')
            main_species = parts[0].strip()
            zpe_keys = []
            gas_species = []

            for p in parts[1:]:
                p = p.strip()
                if p.startswith("zpe_"):
                    zpe_keys.append(p)
                elif p.endswith("(g)") or p.endswith("(s)"):
                    gas_species.append(p.replace("(g)", "_g").replace("(s)", ""))

            main_lookup = main_species.replace("(g)", "_g").replace("(s)", "")
            if main_lookup not in data:
                print(f"⚠️ Warning: Main species '{main_lookup}' not found in data!")
                continue

            # Total energy, ZPE and entropy from main species
            total_energy = data[main_lookup]["total_e"]
            zpe_total = 0.0
            S_total = 0.0
            zpe_terms = []
            S_terms = []

            for z in zpe_keys:
                s = "S_" + z.split("zpe_")[1]
                zval = data[main_lookup]["zpe"].get(z, None)
                sval = data[main_lookup]["S"].get(s, None)
                if zval is None:
                    print(f"⚠️ Missing ZPE '{z}' for species '{main_lookup}'")
                    zval = 0.0
                if sval is None:
                    print(f"⚠️ Missing Entropy '{s}' for species '{main_lookup}'")
                    sval = 0.0
                zpe_total += zval
                S_total += sval
                zpe_terms.append(f"{z}={zval:.6f}")
                S_terms.append(f"{s}={sval:.6f}")

            # Gas-phase species
            energy_terms = [f"{total_energy:.6f} (electronic)"]
            for g in gas_species:
                multiplier = 1
                if '*' in g:
                    try:
                        multiplier_str, g = g.split('*', 1)
                        multiplier = int(multiplier_str.strip())
                    except:
                        print(f"⚠️ Could not parse multiplier in gas species '{g}'")
                        continue

                if g not in data:
                    print(f"⚠️ Gas species '{g}' not found in data!")
                    continue

                ge = data[g]["total_e"]
                gzpe = sum(data[g]["zpe"].values())
                gS = sum(data[g]["S"].values())
                total_energy += multiplier * ge
                zpe_total += multiplier * gzpe
                S_total += multiplier * gS
                energy_terms.append(f"{multiplier}×{ge:.6f} (from {g})")
                zpe_terms.append(f"{multiplier}×{gzpe:.6f} (ZPE {g})")
                S_terms.append(f"{multiplier}×{gS:.6f} (S {g})")

            # Apply ZPE and entropy correction to main energy
            total_energy += zpe_total
            total_energy -= S_total * T

            mechanism.append((line.strip(), total_energy, zpe_total, S_total, energy_terms, zpe_terms, S_terms, [main_lookup] + gas_species))
    return mechanism

def compute_ped(output_file, mechanisms):
    """Generates a potential energy diagram with computed total energies including gas-phase species."""
    with open(output_file, 'w') as f:
        f.write("Mechanism | Reaction | Total Energy (eV)\n")
        f.write("-" * 150 + "\n")
        
        for mech_name, reactions in mechanisms.items():
            reference_energy = reactions[0][1]  # First energy as reference
            ped_energies = []
            reaction_coords = []
            step_number = 0
            
            # Assign only integer steps for all mechanisms
            for reaction, energy, zpe, S, energy_terms, zpe_terms, S_terms, species in reactions:
                ped_energy = energy - reference_energy
                ped_energies.append(f"{ped_energy:.6f}")
                reaction_coords.append(f"{step_number}")
                step_number += 1
                
                f.write(f"{mech_name}\t| {reaction}\n")
                f.write(f"\t= {' + '.join(energy_terms)}\n")
                f.write(f"\t+ ZPE: {' + '.join(zpe_terms)}\n")
                f.write(f"\t- S*T: ({' + '.join(S_terms)}) * {T}\n")  # Updated line to show entropy contribution
                f.write(f"\t= {energy:.6f} eV\n\n")
                
            f.write(f"ped_energies_{mech_name} = {', '.join(ped_energies)}\n")
            f.write(f"rc1_{mech_name} = np.array([{', '.join(reaction_coords)}])\n\n")

if __name__ == "__main__":
    data_files = {
        'lh': read_data_file('data/lh_data.txt'),
        'l2': read_data_file('data/l2_data.txt'),
        'mix': read_data_file('data/mix_data.txt'),
        'mix2': read_data_file('data/mix2_data.txt'),
        'mvk': read_data_file('data/mvk_data.txt')
        }
    mechanisms = {
        'lh': process_mechanism('mech/lh_mech.txt', data_files['lh']),
        'l2': process_mechanism('mech/l2_mech.txt', data_files['l2']),
        'mix': process_mechanism('mech/mix_mech.txt', data_files['mix']),
        'mix2': process_mechanism('mech/mix2_mech.txt', data_files['mix2']),
        'mvk': process_mechanism('mech/mvk_mech.txt', data_files['mvk'])
        }
    
    compute_ped('data/ped.txt', mechanisms)



