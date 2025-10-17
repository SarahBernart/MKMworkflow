import numpy as np
import matplotlib.pyplot as plt
import re

plt.rcParams.update({
    'font.size': 12,
    'axes.titlesize': 12,
    'axes.labelsize': 12,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
})

def read_ped_file(filename):
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

# === Load data ===
ped_file = "temp/ped_450.txt"
ped_data, ts_indices = read_ped_file(ped_file)

# === Define absolute coordinate shifts ===
coordinate_shifts_absolute = {
    "mix": [(5.0, 2.0)],
    "mix2": [(5.0, 2.0), (6.0, 1.0), (7.0, 1.0)],
    "lh":   [(7.0, 2.0)],
    "l2":   [(2.0, 1.0), (4.0, 1.0), (5.0, 2.0)]
    }

def apply_absolute_shifts(ped_data, mech, shifts):
    rc_array = ped_data[mech]["coordinates"]

    for threshold_rc, shift in sorted(shifts, reverse=True):
        print(f"[{mech}] Shift all points ≥ {threshold_rc} by {shift}")
        for i, x in enumerate(rc_array):
            if x >= threshold_rc:
                rc_array[i] += shift

# === Apply shifts ===
for mech, shifts in coordinate_shifts_absolute.items():
    apply_absolute_shifts(ped_data, mech, shifts)

# === Assign variables ===
rc1_lh, set1 = ped_data["lh"]["coordinates"], ped_data["lh"]["energies"]
rc1_l2, set2 = ped_data["l2"]["coordinates"], ped_data["l2"]["energies"]
rc1_mvk, set3 = ped_data["mvk"]["coordinates"], ped_data["mvk"]["energies"]
rc1_mix, set4 = ped_data["mix"]["coordinates"], ped_data["mix"]["energies"]
rc1_mix2, set5 = ped_data["mix2"]["coordinates"], ped_data["mix2"]["energies"]

ts_lh   = ts_indices.get("lh", [])
ts_l2   = ts_indices.get("l2", [])
ts_mvk  = ts_indices.get("mvk", [])
ts_mix  = ts_indices.get("mix", [])
ts_mix2 = ts_indices.get("mix2", [])

print("Detected TS indices:")
print(f"LH: {ts_lh}")
print(f"L2: {ts_l2}")
print(f"MVK: {ts_mvk}")
print(f"MIX: {ts_mix}")
print(f"MIX2: {ts_mix2}")

# === Plotting ===
marker_width = 0.5

def interpol_ts(xxx, yyy):
    if yyy[1] < yyy[0]:  # Spontaneous
        yyy[1] = yyy[2]
    xx = list(np.linspace(xxx[0], xxx[1], 50))
    prefac = yyy[1] - yyy[0]
    yy = [yyy[0] + prefac * np.sin(0.5 * np.pi * (x - xxx[0]) / (xxx[1] - xxx[0])) for x in xx]
    xx2 = list(np.linspace(xxx[1], xxx[2], 50))
    prefac = yyy[2] - yyy[1]
    yy += [yyy[2] - prefac * np.sin(0.5 * np.pi * (x - xxx[0]) / (xxx[2] - xxx[1])) for x in xx2]
    xx += xx2
    return xx, yy

custom_offsets = {
     'mix' : {7: (0.5, 0.0), 8: (0.5, 0.0),10: (0.5, 0.0)},
     'mix2': {6: (0.5,0.0), 8: (0.5,0.0)},
     'lh'  : {5: (-0.5,0.0)},
     'mvk' : {1: (0.5,0.0),3: (-0.5, 0.0)} 
}
default_offset = 0.15

def plot_mechanism(rc, set_values, ts_indices, color, label, mechanism_name):
    for i in range(len(rc)):
        if i not in ts_indices:
            dy, dx = custom_offsets.get(mechanism_name, {}).get(i, (default_offset, 0.0))
            plt.hlines(set_values[i], rc[i] - marker_width / 2, rc[i] + marker_width / 2,
                       color=color, linewidth=2)
            plt.text(rc[i] + dx, set_values[i] + dy, f"{set_values[i]:.2f}",
                     color=color, fontsize=11, ha='center')

    for i in range(1, len(rc)):
        if i in ts_indices:
            if i - 1 < 0 or i + 1 >= len(rc):
                continue  # skip TS at edges
            dy, dx = custom_offsets.get(mechanism_name, {}).get(i, (default_offset, 0.0))
            rc_mid = (rc[i - 1] + rc[i + 1]) / 2
            xxx = [rc[i - 1] + marker_width / 2, rc_mid, rc[i + 1] - marker_width / 2]
            yyy = [set_values[i - 1], set_values[i], set_values[i + 1]]
            spontaneous = set_values[i] < set_values[i - 1]
            xx, yy = interpol_ts(xxx, yyy)
            plt.plot(xx, yy, color=color, linewidth=2)
            plt.plot(rc_mid + dx, set_values[i], 'o', color=color, markersize=0.01)
            if not spontaneous:
                plt.text(rc_mid, set_values[i] + dy - 0.1, f"{set_values[i]:.2f}",
                         color=color, fontsize=11, ha='center')
        elif i - 1 in ts_indices:
            continue
        else:
            plt.plot([rc[i - 1] + marker_width / 2, rc[i] - marker_width / 2],
                     [set_values[i - 1], set_values[i]], color=color, linestyle='dashed')

# === Final plot ===
plt.figure(figsize=(8, 6))
plot_mechanism(rc1_mix2, set5, ts_mix2, '#47D45A',  'Mix2', 'mix2')
plot_mechanism(rc1_lh,   set1, ts_lh,   '#83CBEB',  'LH, stepwise', 'lh')
plot_mechanism(rc1_l2,   set2, ts_l2,   '#4344E4',  'LH, concerted', 'l2')
plot_mechanism(rc1_mvk,  set3, ts_mvk,  '#FFA032',  'MvK', 'mvk')
plot_mechanism(rc1_mix,  set4, ts_mix,  '#7F7F7F',  'Mix', 'mix')

plt.xlabel('Reaction coordinate')
plt.ylabel(r'$\Delta G$ [eV]')
plt.title(r'$Pd_{10}/CeO_{2}(111)$')
plt.ylim(-7.5, 5)

plt.plot([], [], color='#4344E4',   label='LH, conc.', linewidth=2)
plt.plot([], [], color='#FFA032',   label='MvK', linewidth=2)
plt.plot([], [], color='#7F7F7F',   label='Mix, step.', linewidth=2)
plt.plot([], [], color='#47D45A',   label='Mix, conc.', linewidth=2)
plt.plot([], [], color='#83CBEB',   label='LH, step.', linewidth=2)
plt.legend()

plt.gca().set_xticklabels([])
plt.gca().set_xticks([])

plt.tight_layout()
#plt.savefig('data/Pd1.pdf')
plt.show()


