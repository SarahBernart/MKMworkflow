#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# ---------- Style only (no content changes) ----------
plt.rcParams['font.size'] = 12
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['legend.fontsize'] = 12

# ---------- Mechanism subfolders (same as your single-plot script) ----------
MECH_SUBFOLDERS = [
    "Kinetics_L2",
    "Kinetics_Lh",
    "Kinetics_Mix",
    "Kinetics_Mix2",
    "Kinetics_MvK",
    "Kinetics_All",
]

# ---------- Systems (positions only; content unchanged) ----------
SYSTEMS = [
    ("../SACs_edge/Pd1_CeO2",  r"$Pd_{1}$"),  # top-left
    ("../Cluster/Pd4_CeO2",    r"$Pd_{4}$"),  # top-middle
    ("../Cluster/Pd10_CeO2",   r"$Pd_{10}$"), # top-right
    ("../SACs_edge/Pt1_CeO2",  r"$Pt_{1}$"),  # bottom-left
    ("../Cluster/Pt4_CeO2",    r"$Pt_{4}$"),  # bottom-middle
    ("../Cluster/Pt10_CeO2",   r"$Pt_{10}$"), # bottom-right
]

# ---------- Legend labels / colors / linestyles (unchanged content) ----------
legend_labels = {
    "MIX":  "Mix, step.",
    "MIX2": "Mix, conc.",
    "LH":   "LH, step.",
    "L2":   "LH, conc.",
    "MVK":  "MvK",
    "ALL":  "Combined"
}

color_map = {
    "MIX":  "#7F7F7F",
    "MIX2": "#47D45A",
    "LH":   "#83CBEB",
    "L2":   "#4344E4",
    "MVK":  "#FFA032",
    "ALL":  "black"
}

linestyle_map = {
    "ALL": "-.",
    "MIX2": "-."
}

# ---------- YOU control draw/legend order here ----------
# Keys are the folder suffixes turned into mechanism codes below.
# "ALL" will be drawn last **regardless**, to stay on top.
MECH_ORDER = ["L2", "LH", "MVK", "MIX", "MIX2", "ALL"]

# ---------- Figure / grid ----------
fig, axs = plt.subplots(2, 3, figsize=(11.5, 7.5), sharex=True)
axs = axs.flatten()
fig.subplots_adjust(hspace=0.28, wspace=0.10)

# Global legend map: one handle per mechanism (first encountered)
mech_to_handle = {}
present_mechs = set()

for i, (root, panel_title) in enumerate(SYSTEMS):
    ax = axs[i]
    combined_df = pd.DataFrame()

    # collect data like your single plot
    folders = [os.path.join(root, m) for m in MECH_SUBFOLDERS]
    for folder in folders:
        range_path = os.path.join(folder, "range")
        if os.path.isdir(range_path):
            derivatives_file = os.path.join(range_path, "derivatives.dat")
            if os.path.exists(derivatives_file):
                df = pd.read_csv(derivatives_file, delimiter="\t")
                df.columns = [col.strip() for col in df.columns]
                df = df.apply(pd.to_numeric, errors='coerce')

                mech_code = os.path.basename(folder).replace("Kinetics_", "").upper()  # -> L2, LH, MIX, MIX2, MVK, ALL
                df["Mechanism"] = mech_code
                df["Source"] = legend_labels.get(mech_code, mech_code)

                combined_df = pd.concat([combined_df, df], ignore_index=True)

    # same content logic as your example
    if "CO2" in combined_df.columns and "Temperature" in combined_df.columns:
        combined_df = combined_df.dropna(subset=["Temperature", "CO2"])
        combined_df = combined_df.sort_values(by="Temperature")
        combined_df["Log_CO2"] = np.log10(combined_df["CO2"].replace(0, np.nan))
        combined_df["Log_CO2"] = combined_df["Log_CO2"].replace([np.inf, -np.inf], np.nan)
        combined_df = combined_df.dropna(subset=["Log_CO2"])
        combined_df = combined_df[(combined_df["Temperature"] >= 400) & (combined_df["Temperature"] <= 1000)]

        # ---- Plot order per panel ----
        mechs_here = list(combined_df["Mechanism"].unique())
        # non-ALL in requested order, then any remaining unknowns (except ALL)
        ordered_non_all = [m for m in MECH_ORDER if m != "ALL" and m in mechs_here]
        others_non_all  = [m for m in mechs_here if m not in MECH_ORDER and m != "ALL"]
        draw_order = ordered_non_all + others_non_all

        # 1) draw all non-ALL first
        for mech in draw_order:
            data = combined_df[combined_df["Mechanism"] == mech]
            label = data["Source"].iloc[0]
            color = color_map.get(mech, "black")
            linestyle = linestyle_map.get(mech, "-")
            ln, = ax.plot(
                data["Temperature"], data["Log_CO2"],
                label=label, color=color, linestyle=linestyle, linewidth=1.6, zorder=2
            )
            # remember first handle for global legend (by mechanism)
            if mech not in mech_to_handle:
                mech_to_handle[mech] = ln
                present_mechs.add(mech)

        # 2) draw ALL last and on top (if present)
        if "ALL" in mechs_here:
            data = combined_df[combined_df["Mechanism"] == "ALL"]
            label = data["Source"].iloc[0]
            color = color_map.get("ALL", "black")
            linestyle = linestyle_map.get("ALL", "-.")
            ln, = ax.plot(
                data["Temperature"], data["Log_CO2"],
                label=label, color=color, linestyle=linestyle, linewidth=1.6, zorder=10
            )
            if "ALL" not in mech_to_handle:
                mech_to_handle["ALL"] = ln
                present_mechs.add("ALL")

        # Titles/axes (content unchanged)
        ax.set_title(panel_title)
        ax.axhline(0, color='black', linewidth=0.75, linestyle='--')
        ax.set_ylim(-15, 10)
        ax.grid(linewidth=0.45)
        if i % 3 == 0:
            ax.set_ylabel("$\mathrm{TOF_{log}(CO_{2})}$ [s$\mathrm{^{-1}}$]")
        else:
            ax.tick_params(labelleft=False)
        if i >= 3:
            ax.set_xlabel("Temperature [K]")
        else:
            ax.tick_params(labelbottom=False)
    else:
        ax.set_title(panel_title)
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
        if i % 3 == 0:
            ax.set_ylabel("$\mathrm{TOF_{log}(CO_{2})}$ [s$\mathrm{^{-1}}$]")
        if i >= 3:
            ax.set_xlabel("Temperature [K]")

# ---------- Global legend in your chosen order (ALL last) ----------
present_order = [m for m in MECH_ORDER if m in present_mechs] + \
                [m for m in mech_to_handle.keys() if m not in MECH_ORDER]

handles_all = [mech_to_handle[m] for m in present_order]
labels_all  = [legend_labels.get(m, m) for m in present_order]

plt.tight_layout(rect=[0.02, 0.12, 1.00, 1.00]) 
legend_ax = fig.add_axes([0.00, 0.00, 1.00, 0.2])  #[left, bottom, width, height] 
legend_ax.axis("off") 
leg = legend_ax.legend( 
                       handles_all, 
                       labels_all, 
                       loc="center", 
                       ncol=min(6, 
                                len(labels_all)), 
                       frameon=True, # Rahmen einschalten 
                       fancybox=True, # abgerundete Ecken 
                       fontsize=12, 
                       handlelength=1.8, 
                       columnspacing=0.8, 
                       handletextpad=0.3, 
                       borderpad=0.4, # Innenabstand Box ↔ Inhalt 
                       labelspacing=0.3 )

os.makedirs("plots", exist_ok=True)
plt.savefig("plots/all_mech_grid.png", dpi=300)
plt.close()
print("Saved plots/all_mech_grid.png")

