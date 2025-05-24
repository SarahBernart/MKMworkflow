import numpy as np
import matplotlib.pyplot as plt
import re
import os
import time

def read_ped_file(filename):
    """Reads ped.txt and extracts reaction coordinates, energies, and identifies TS positions."""
    ped_data = {}
    ts_indices = {}
    reaction_order = {}

    with open(filename, "r") as f:
        lines = f.readlines()

    for line in lines:
        if line.startswith("ped_energies_"):
            mechanism = line.split("=")[0].strip().replace("ped_energies_", "")
            ped_data[mechanism] = {
                "energies": np.array([float(x) for x in line.split("=")[1].strip().split(",")])
            }
        elif line.startswith("rc1_"):
            mechanism = line.split("=")[0].strip().replace("rc1_", "")
            match = re.search(r"\[(.*?)\]", line)
            if match:
                numbers = match.group(1).split(",")
                ped_data[mechanism]["coordinates"] = np.array([float(x.strip()) for x in numbers])
                reaction_order[mechanism] = []
            else:
                raise ValueError(f"Could not parse reaction coordinates in line: {line}")

    for mechanism in reaction_order.keys():
        index = 0
        for line in lines:
            if f"{mechanism}\t|" in line:
                reaction_order[mechanism].append((line.strip(), index))
                index += 1
        ts_indices[mechanism] = [idx for name, idx in reaction_order[mechanism] if "TS_" in name]

    return ped_data, ts_indices

temp_folder = "temp"
ped_files = sorted([f for f in os.listdir(temp_folder) if f.startswith("ped_") and f.endswith(".txt")],
                   key=lambda x: int(x.split("_")[1].split(".")[0]))

plt.ion()
fig, ax = plt.subplots(figsize=(8, 3))

while True:
    for ped_file in ped_files:
        filepath = os.path.join(temp_folder, ped_file)
        temperature = ped_file.split("_")[1].split(".")[0]
        ped_data, ts_indices = read_ped_file(filepath)

        ax.clear()
        marker_width = 0.5

        ts_lh = ts_indices.get("lh", [])
        ts_mvk = ts_indices.get("mvk", [])
        ts_mix = ts_indices.get("mix", [])

        def interpol_ts(xxx, yyy):
            if yyy[1] < yyy[0]:
                yyy[1] = yyy[2]

            xx = list(np.linspace(xxx[0], xxx[1], 50))
            prefac = yyy[1] - yyy[0]
            yy = [yyy[0] + prefac * np.sin(0.5 * np.pi * (x - xxx[0]) / (xxx[1] - xxx[0])) for x in xx]
            xx2 = list(np.linspace(xxx[1], xxx[2], 50))
            prefac = yyy[2] - yyy[1]
            yy += [yyy[2] - prefac * np.sin(0.5 * np.pi * (x - xxx[0]) / (xxx[2] - xxx[1])) for x in xx2]
            xx += xx2
            return xx, yy

        def plot_mechanism(rc, set_values, ts_indices, color, label):
            for i in range(len(rc)):
                if i not in ts_indices:
                    ax.hlines(set_values[i], rc[i] - marker_width / 2, rc[i] + marker_width / 2, color=color, linewidth=2)
            for i in range(1, len(rc)):
                if i in ts_indices:
                    xxx = [rc[i-1] + marker_width / 2, rc[i], rc[i+1] - marker_width / 2]
                    yyy = [set_values[i-1], set_values[i], set_values[i+1]]
                    xx, yy = interpol_ts(xxx, yyy)
                    ax.plot(xx, yy, color=color, linewidth=2)
                    ax.plot(rc[i], set_values[i], 'o', color=color, markersize=0.01)
                elif i - 1 in ts_indices:
                    continue
                else:
                    ax.plot([rc[i-1] + marker_width / 2, rc[i] - marker_width / 2],
                             [set_values[i-1], set_values[i]], color=color, linestyle='dashed')
    
        plot_mechanism(ped_data["lh"]["coordinates"], ped_data["lh"]["energies"], ts_lh, 'lightblue', 'LH')
        plot_mechanism(ped_data["mvk"]["coordinates"], ped_data["mvk"]["energies"], ts_mvk, 'black', 'MVK')
        plot_mechanism(ped_data["mix"]["coordinates"], ped_data["mix"]["energies"], ts_mix, 'gray', 'MIX')

        ax.set_xlabel("Reaction Coordinate")
        ax.set_ylabel("ΔG [eV]")
        ax.set_title(f"CO Oxidation (Mixed route), M1 at {temperature} K")
        ax.legend(title=f"Temperature: {temperature} K")
        ax.set_ylim(-8, 3)
        plt.pause(0.001)

plt.ioff()
plt.close(fig)
print("All temperature plots displayed successfully.")
plt.close() 
