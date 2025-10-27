#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot DRC (Degree of Rate Control) for Pd_n/Pt_n on CeO2(111)
- 2x3 Layout: oben Pd1, Pd4, Pd10; unten Pt1, Pt4, Pt10
- Legende unten (eigene, unsichtbare Achse – wird nicht abgeschnitten)

Outputs:
- plots/all.png
- plots/all.txt (DRC-Export über definierte Temperaturen)
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# ---------- kompaktere Schriftgrößen (nur Stil) ----------
plt.rcParams.update({
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
})

# === Mechanismen aus externer Datei laden ===
mechanism_file = "mech.txt"
farbe_rgb_map = {}
reaction_colors = {}
reaction_labels = {}

with open(mechanism_file, "r", encoding="utf-8") as f:
    code = f.read()
    local_vars = {}
    exec(code, {}, local_vars)
    farbe_rgb_map = local_vars.get("farbe_rgb_map", {})
    reaction_colors = local_vars.get("reaction_colors", {})
    reaction_labels = local_vars.get("reaction_labels", {})

# === Reaktionsstile automatisch ableiten ===
reaction_styles = {}
for reaction in reaction_colors:
    if "OCOO" in reaction and "CO2" in reaction:
        reaction_styles[reaction] = (0, (5, 2))        # OCOO → CO2 (lange-kurz)
    elif "OCOO" in reaction:
        reaction_styles[reaction] = ':'                # concerted Mechanismus
    elif "CO" in reaction and "O2" in reaction:
        reaction_styles[reaction] = '-.'               # CO Oxidation
    elif "CO2" in reaction:
        reaction_styles[reaction] = (0, (2, 2, 6, 2))  # CO2 Desorption
    elif "O2" in reaction and "+" not in reaction:
        reaction_styles[reaction] = ':'                # O2 Dissoziation
    elif "Ovic_Vosub" in reaction:
        reaction_styles[reaction] = (0, (1, 1))        # Rückreaktion
    else:
        reaction_styles[reaction] = '-'                # Default

# === Ordnerstruktur und Plot-Zuordnung ===
folders = [
    "~/calc/Kinetics/SACs_edge/Pd1_CeO2/Kinetics_All", 
    "~/calc/Kinetics/SACs_edge/Pt1_CeO2/Kinetics_All",
    "~/calc/Kinetics/Cluster/Pd4_CeO2/Kinetics_All",
    "~/calc/Kinetics/Cluster/Pt4_CeO2/Kinetics_All",
    "~/calc/Kinetics/Cluster/Pd10_CeO2/Kinetics_All",
    "~/calc/Kinetics/Cluster/Pt10_CeO2/Kinetics_All"
]
folders = [os.path.expanduser(f) for f in folders]

# gewünschte Anordnung: oben Pd1, Pd4, Pd10; unten Pt1, Pt4, Pt10
subplot_order = [
    "$Pd_{1}$", "$Pd_{4}$", "$Pd_{10}$",   # obere Reihe (links→rechts)
    "$Pt_{1}$", "$Pt_{4}$", "$Pt_{10}$"    # untere Reihe (links→rechts)
]

legend_labels = {
    os.path.expanduser("~/calc/Kinetics/SACs_edge/Pd1_CeO2/Kinetics_All"): "$Pd_{1}$",
    os.path.expanduser("~/calc/Kinetics/SACs_edge/Pt1_CeO2/Kinetics_All"): "$Pt_{1}$",
    os.path.expanduser("~/calc/Kinetics/Cluster/Pd4_CeO2/Kinetics_All"):   "$Pd_{4}$",
    os.path.expanduser("~/calc/Kinetics/Cluster/Pt4_CeO2/Kinetics_All"):   "$Pt_{4}$",
    os.path.expanduser("~/calc/Kinetics/Cluster/Pd10_CeO2/Kinetics_All"):  "$Pd_{10}$",
    os.path.expanduser("~/calc/Kinetics/Cluster/Pt10_CeO2/Kinetics_All"):  "$Pt_{10}$"
}

# === Daten laden ===
min_temp = 350
max_temp = 1000
min_abs_drc = 0.05   # Linien mit |DRC| < 0.05 über den gesamten Bereich ausblenden

combined_df = pd.DataFrame()

for folder in folders:
    drc_file = os.path.join(folder, "range", "drc", "drc.dat")
    if os.path.exists(drc_file):
        df = pd.read_csv(drc_file, delimiter="\t")
        df.columns = [col.strip() for col in df.columns]
        df = df.apply(pd.to_numeric, errors='coerce')
        df["Source"] = legend_labels.get(folder, folder)
        combined_df = pd.concat([combined_df, df], ignore_index=True)

# === Plot erstellen ===
if not combined_df.empty:
    combined_df = combined_df[
        (combined_df["Temperature"] >= min_temp) & (combined_df["Temperature"] <= max_temp)
    ]

    export_temps = list(range(100, 1101, 100))
    export_folder = "plots/drc_export"
    os.makedirs(export_folder, exist_ok=True)

    all_lines = []
    header_all = ["Mechanism", "Species"] + [str(T) for T in export_temps]
    all_lines.append("\t".join(header_all))

    for source in combined_df["Source"].unique():
        source_df = combined_df[combined_df["Source"] == source]
        species_columns = [col for col in source_df.columns if col not in ["Temperature", "Source"]]

        # Mechanismus-Name kürzen für Dateinamen
        mech_short = source.replace(",", "").replace(" ", "_").replace(".", "")
        
        # === Einzeldokument für Mechanismus (optional – aktuell auskommentiert unten) ===
        output_lines = []
        header = ["Species"] + [str(T) for T in export_temps]
        output_lines.append("\t".join(header))

        for species in species_columns:
            line = [species]
            all_line = [mech_short, species]

            for T in export_temps:
                row = source_df[np.isclose(source_df["Temperature"], T)]
                if not row.empty:
                    value = row[species].values[0]
                    line.append(f"{value:.6e}")
                    all_line.append(f"{value:.6e}")
                else:
                    line.append("NaN")
                    all_line.append("NaN")

            output_lines.append("\t".join(line))
            all_lines.append("\t".join(all_line))

        # # Speichern einzelner Mechanismus (falls gewünscht)
        # outpath = os.path.join(export_folder, f"DRC_{mech_short}.txt")
        # with open(outpath, "w") as f:
        #     f.write("\n".join(output_lines))
        # print(f"Exported DRC data to {outpath}")

    # === Alle Mechanismen gesammelt ===
    os.makedirs("plots", exist_ok=True)
    outpath_all = os.path.join("plots", "all.txt")
    with open(outpath_all, "w") as f:
        f.write("\n".join(all_lines))
    print(f"Combined DRC export written to {outpath_all}")

    # === Plot ===
    high_res_temp = np.linspace(min_temp, max_temp, 1000)

    # größere Figure, gleiche Inhalte
    fig, axs = plt.subplots(2, 3, figsize=(11.5, 7.5), sharex=True)
    axs = axs.flatten()
    fig.subplots_adjust(hspace=0.28, wspace=0.10)

    for i, source in enumerate(subplot_order[:6]):
        ax = axs[i]
        source_df = combined_df[combined_df["Source"] == source]

        species_columns = [col for col in source_df.columns if col not in ["Temperature", "Source"]]
        # optionale Schwelle beibehalten
        species_columns = [s for s in species_columns if source_df[s].abs().max() >= 0.01]
        species_columns = [s for s in species_columns if source_df[s].abs().max() >= min_abs_drc]

        for species in species_columns:
            label = reaction_labels.get(species, species)
            colorname = reaction_colors.get(species, "grey")
            linestyle = reaction_styles.get(species, '-')
            rgb = farbe_rgb_map.get(colorname, "#000000")

            interpolation = interp1d(
                source_df["Temperature"], source_df[species],
                kind='cubic', fill_value='extrapolate'
            )
            smooth_values = interpolation(high_res_temp)
            ax.plot(
                high_res_temp, smooth_values,
                label=label, color=rgb, linestyle=linestyle, linewidth=1.4   # etwas dünner
            )

        # === Plotgestaltung je nach Index i ===
        ax.set_title(rf"{source}/CeO$_2$(111)", fontsize=9, pad=4)  # Content unverändert
        if i in [0, 3]:
            ax.set_ylabel("DRC", fontsize=8)                         # Content unverändert
        else:
            ax.tick_params(labelleft=False)

        if i in [3, 4, 5]:
            ax.set_xlabel("Temperature [K]", fontsize=12)             # Content unverändert
            ax.tick_params(labelbottom=True)
        else:
            ax.tick_params(labelbottom=False)

        ax.grid(True, linewidth=0.5, alpha=0.6)
        ax.set_ylim(-0.05, 1.05)
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.title.set_position([.5, 1.01])

    fig.suptitle("DRC of $Pd_{n}/CeO_{2}(111)$", fontsize=9, y=0.97)  # Content unverändert

    # === Legende (robust, kleine Schrift, nicht abgeschnitten) ===
    used_labels = set()
    handles, labels = [], []
    for ax in axs:
        for line in ax.get_lines():
            lab = line.get_label()
            if lab not in used_labels and not lab.startswith('_'):
                handles.append(line)
                labels.append(lab)
                used_labels.add(lab)

    # Platz unten für Legende
    plt.tight_layout(rect=[0, 0.16, 1, 0.97])

    # eigene, unsichtbare Achse für die Legende
    legend_ax = fig.add_axes([0.0, 0.0, 1.0, 0.2])  # [left, bottom, width, height]
    legend_ax.axis("off")
    legend_ax.legend(
        handles, labels,
        loc="center",
        ncol=7,
        frameon=True,
        fancybox=True,
        fontsize=10,         # kleiner
        handlelength=1.8,   # kürzer
        columnspacing=0.7,
        handletextpad=0.3,
        borderpad=0.2, 
        labelspacing=0.2
    )
    fig.suptitle("DRC of $Pd_{n}/CeO_{2}(111)$", fontsize=12, y=0.98)
    os.makedirs("plots", exist_ok=True)
    plt.savefig("plots/all.png", format="png", dpi=600)
    plt.close()
    print("Ordered combined plot saved as plots/all.png")

else:
    print("No data found.")

