import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
import re

# === Plot style ===
plt.rcParams.update({
    'font.size': 12,
    'axes.titlesize': 12,
    'axes.labelsize': 12,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
})

# === Read PED file ===
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

    for mechanism in reaction_order.keys():
        index = 0
        for line in lines:
            if f"{mechanism}\t|" in line:
                reaction_order[mechanism].append((line.strip(), index))
                index += 1
        ts_indices[mechanism] = [idx for name, idx in reaction_order[mechanism] if "TS_" in name]

    return ped_data, ts_indices

# === Coordinate shifts ===
coordinate_shifts_absolute = {
    "mix": [(5.0, 1.0)],
    "mix2": [(5.0, 1.0), (6.0, 1.0), (7.0, 1.0)],
    "lh":   [(7.0, 2.0)],
    "l2":   [(2.0, 1.0), (4.0, 1.0), (5.0, 1.0)]
}

def apply_absolute_shifts(ped_data, mech, shifts):
    rc_array = ped_data[mech]["coordinates"]
    for threshold_rc, shift in sorted(shifts, reverse=True):
        for i, x in enumerate(rc_array):
            if x >= threshold_rc:
                rc_array[i] += shift

# === Sine interpolation (your version) ===
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

# === Plotting with text and exact styling ===
marker_width = 0.5
default_offset = 0.15
custom_offsets = {
    'mix2': {6: (0.25, 0.0)},
}

def plot_mechanism(ax, rc, set_values, ts_indices, color, label, mech_name):
    for i in range(len(rc)):
        if i not in ts_indices:
            dy, dx = custom_offsets.get(mech_name, {}).get(i, (default_offset, 0.0))
            ax.hlines(set_values[i], rc[i] - marker_width / 2, rc[i] + marker_width / 2, color=color, linewidth=2)
            ax.text(rc[i] + dx, set_values[i] + dy, f"{set_values[i]:.2f}", color=color,
                    fontsize=11, ha='center')

    for i in range(1, len(rc)):
        if i in ts_indices:
            if i - 1 < 0 or i + 1 >= len(rc):
                continue
            dy, dx = custom_offsets.get(mech_name, {}).get(i, (default_offset, 0.0))
            rc_mid = (rc[i - 1] + rc[i + 1]) / 2
            xxx = [rc[i - 1] + marker_width / 2, rc_mid, rc[i + 1] - marker_width / 2]
            yyy = [set_values[i - 1], set_values[i], set_values[i + 1]]
            spontaneous = set_values[i] < set_values[i - 1]
            xx, yy = interpol_ts(xxx, yyy)
            ax.plot(xx, yy, color=color, linewidth=2)
            ax.plot(rc_mid + dx, set_values[i], 'o', color=color, markersize=0.01)
            if not spontaneous:
                ax.text(rc_mid, set_values[i] + dy - 0.1, f"{set_values[i]:.2f}",
                        color=color, fontsize=11, ha='center')
        elif i - 1 in ts_indices:
            continue
        else:
            ax.plot([rc[i - 1] + marker_width / 2, rc[i] - marker_width / 2],
                    [set_values[i - 1], set_values[i]], color=color, linestyle='dashed')

# === Load temperature files ===
temp_folder = "temp"
ped_files = sorted([f for f in os.listdir(temp_folder) if f.startswith("ped_") and f.endswith(".txt")],
                   key=lambda x: int(x.split("_")[1].split(".")[0]))

fig, ax = plt.subplots(figsize=(10, 7))

def create_frame(i):
    ax.clear()
    ped_file = ped_files[i]
    temperature = ped_file.split("_")[1].split(".")[0]
    filepath = os.path.join(temp_folder, ped_file)
    ped_data, ts_indices = read_ped_file(filepath)

    # Apply shifts
    for mech, shifts in coordinate_shifts_absolute.items():
        if mech in ped_data:
            apply_absolute_shifts(ped_data, mech, shifts)

    # Plot each mechanism
    mechanisms = [
        ("mix2", "#47D45A", 'Mix2'),
        ("l2",   "#4344E4", 'LH, concerted'),
        ("lh",   "#83CBEB", 'LH, stepwise'),
        ("mvk",  "#FFA032", 'MvK'),
        ("mix",  "#7F7F7F", 'Mix'),
    ]

    for mech, color, label in mechanisms:
        if mech in ped_data:
            plot_mechanism(ax, ped_data[mech]["coordinates"], ped_data[mech]["energies"],
                           ts_indices.get(mech, []), color, label, mech)

    ax.set_xlabel("Reaction coordinate")
    ax.set_ylabel(r'$\Delta G$ [eV]')
    ax.set_title(rf'$Pd_4/CeO_2(111)$ – {temperature} K')
    ax.set_ylim(-7.5, 5)
    ax.set_xticks([])
    ax.set_xticklabels([])

    for mech, color, label in mechanisms:
        ax.plot([], [], color=color, label=label, linewidth=2)
    ax.legend()

    return fig,

# === Create animation ===
ani = animation.FuncAnimation(fig, create_frame, frames=len(ped_files), blit=False)

# === Save video ===
os.makedirs("movie", exist_ok=True)
ani.save("movie/Pd4_movie.mp4", writer='ffmpeg', fps=20, dpi=200)

plt.close(fig)
print("Movie saved successfully to movie/Pd4_movie.mp4")

