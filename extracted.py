import os
import math

def calculate_entropy(frequencies, temperature=450):
    k_B = 8.617333262e-5
    h = 4.135667696e-15
    c = 2.99792e10
    entropy = 0.0
    for freq in frequencies:
        freq_eV = h * c * freq
        x = freq_eV / (k_B * temperature)
        if x > 1e-6:
            entropy += (x / (math.exp(x) - 1)) - math.log1p(-math.exp(-x))
    return f"{entropy * k_B:.3e}"

def read_energy(path):
    try:
        with open(path, "r") as f:
            return round(float(f.readline().strip()), 3)
    except:
        return None

def read_zpe_and_entropy(folder, zpe_results, entropy_results, is_gas):
    for sub in os.listdir(folder):
        path = os.path.join(folder, sub)
        if not os.path.isdir(path) or not sub.startswith("zpe_"):
            continue
        freq_path = os.path.join(path, "frequencies.txt")
        if not os.path.isfile(freq_path):
            continue
        freqs, zpe = [], None
        with open(freq_path) as f:
            for line in f:
                if "Zero-point energy:" in line:
                    for part in line.split():
                        try:
                            zpe = float(part)
                            break
                        except:
                            continue
                elif len(line.split()) == 3 and line.split()[0].isdigit():
                    try:
                        freqs.append(float(line.split()[2]))
                    except:
                        continue
        if zpe is not None:
            zpe_results.setdefault(folder, {})[sub] = zpe
        if is_gas:
            ent_file = os.path.join(folder, "only_entropy.txt")
            if os.path.isfile(ent_file):
                try:
                    with open(ent_file) as ef:
                        val = float(ef.readline().strip())
                        entropy_results.setdefault(folder, {})["S_thermo"] = f"{val:.3e}"
                except:
                    pass
        elif freqs:
            entropy_key = "S_" + sub.replace("zpe_", "")
            entropy_results.setdefault(folder, {})[entropy_key] = calculate_entropy(freqs)

def extract_total_energy(folders, gas_folders, out="data.txt"):
    all_folders = folders + gas_folders
    results, zpe, entropy = {}, {}, {}
    for folder in all_folders:
        if not os.path.isdir(folder):
            continue
        if folder.startswith("TS_"):
            for sub in os.listdir(folder):
                path = os.path.join(folder, sub)
                if os.path.isdir(path) and sub.endswith("_TS"):
                    e = read_energy(os.path.join(path, "total.e"))
                    if e is not None:
                        results[path] = e
                    read_zpe_and_entropy(path, zpe, entropy, False)
        else:
            e = read_energy(os.path.join(folder, "total.e"))
            if e is not None:
                results[folder] = e
            read_zpe_and_entropy(folder, zpe, entropy, folder in gas_folders)

    os.makedirs("data", exist_ok=True)
    keys = sorted(set(k for v in zpe.values() for k in v) | set("zpe_" + k[2:] for v in entropy.values() for k in v if k.startswith("S_")))
    out_path = os.path.join("data", out)
    with open(out_path, "w") as f:
        f.write("Species\tTotal_E\t" + "\t".join(keys) + "\t" + "\t".join("S_" + k.replace("zpe_", "") for k in keys) + "\n")
        combined = {}
        for k in results:
            main = k.split("/")[0] if "/" in k and k.startswith("TS_") else k
            if main not in combined:
                combined[main] = {"e": results[k], "z": {}, "s": {}}
            for z in keys:
                combined[main]["z"][z] = zpe.get(k, {}).get(z, "-")
                s_key = "S_" + z.replace("zpe_", "")
                combined[main]["s"][z] = entropy.get(k, {}).get(s_key, "-")
        for k, v in combined.items():
            zpe_str = "\t".join(str(v["z"].get(z, "-")) for z in keys)
            S_str = "\t".join(str(v["s"].get(z, "-")) for z in keys)
            f.write(f"{k}\t{v['e']}\t{zpe_str}\t{S_str}\n")
    print(f"✅ Saved to {out_path}")

def read_species_file(filepath="mech/species.txt"):
    species_data = {
        "gasphase_folders": [], "mixed_folders": [],
        "mix2_folders": [], "lh_folders": [],
        "l2_folders": [], "mvk_folders": []
    }
    try:
        with open(filepath, "r") as file:
            current_key = None
            for line in file:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                if line in species_data:
                    current_key = line
                elif current_key:
                    species_data[current_key].append(line)
    except FileNotFoundError:
        print(f"❌ Error: File {filepath} not found.")

    return (
        species_data["gasphase_folders"],
        species_data["mixed_folders"],
        species_data["mix2_folders"],
        species_data["lh_folders"],
        species_data["l2_folders"],
        species_data["mvk_folders"]
    )

if __name__ == "__main__":
    gas, mix, mix2, lh, l2, mvk = read_species_file()

    print("Processing Mixed Data:")
    extract_total_energy(mix, gas, out="mix_data.txt")

    print("Processing Mix2 Data:")
    extract_total_energy(mix2, gas, out="mix2_data.txt")

    print("Processing LH Data:")
    extract_total_energy(lh, gas, out="lh_data.txt")

    print("Processing L2 Data:")
    extract_total_energy(l2, gas, out="l2_data.txt")

    print("Processing MvK Data:")
    extract_total_energy(mvk, gas, out="mvk_data.txt")


