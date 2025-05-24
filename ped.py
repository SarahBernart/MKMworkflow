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
    """Processes a mechanism file to compute total energies, including gas-phase species."""
    mechanism = []
    with open(mech_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            parts = line.split('+')
            species = [p.strip() for p in parts]
            total_energy = 0.0
            total_zpe = 0.0
            total_S = 0.0
            energy_terms = []
            zpe_terms = []
            S_terms = []
            
            for sp in species:
                multiplier = 1
                if '*' in sp:  # Handle multipliers like 2*CO
                    multiplier, sp = sp.split('*')
                    multiplier = int(multiplier.strip())
                    sp = sp.strip()
                
                sp_lookup = sp.replace("(g)", "_g") if "(g)" in sp else sp.replace("(s)", "")  # Match correctly
                if sp_lookup in data:
                    energy_terms.append(f"({data[sp_lookup]['total_e']:.6f} * {multiplier})")
                    total_energy += multiplier * data[sp_lookup]['total_e']
                    
                    zpe_terms.append(f"({sum(data[sp_lookup]['zpe'].values()):.6f} * {multiplier})")
                    total_zpe += multiplier * sum(data[sp_lookup]['zpe'].values())
                    
                    S_terms.append(f"({sum(data[sp_lookup]['S'].values()):.6f} * {multiplier})")
                    total_S += multiplier * sum(data[sp_lookup]['S'].values())
                elif sp_lookup.startswith("zpe_"):
                    continue  # Ignore zpe_ species from mech file as they are not actual species
                else:
                    print(f"Warning: Species {sp_lookup} not found in data!")
            
            total_energy += total_zpe  # Adding ZPE correction
            total_energy += total_S * T  # Subtract entropy contribution correctly

            mechanism.append((line.strip(), total_energy, total_zpe, total_S, energy_terms, zpe_terms, S_terms, species))
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

