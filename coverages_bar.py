import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.interpolate import interp1d

# ---------- Style only (keine Inhaltsänderungen) ----------
plt.rcParams['font.size'] = 12
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['legend.fontsize'] = 12

# === Mechanismen, Farben, Styles aus mech.txt laden ===
mechanism_file = "mech.txt"
farbe_rgb_map = {}
reaction_colors = {}
reaction_labels = {}
reaction_styles = {}

with open(mechanism_file, "r", encoding="utf-8") as f:
    code = f.read()
    local_vars = {}
    exec(code, {}, local_vars)
    farbe_rgb_map = local_vars.get("farbe_rgb_map", {})
    reaction_colors = local_vars.get("reaction_colors", {})
    reaction_labels = local_vars.get("reaction_labels", {})
    reaction_styles = local_vars.get("reaction_styles", {})

# === Ordner und Labels ===
folders = [
    "~/calc/Kinetics/SACs_edge/Pd1_CeO2/Kinetics_All", 
    "~/calc/Kinetics/SACs_edge/Pt1_CeO2/Kinetics_All",
    "~/calc/Kinetics/Cluster/Pd4_CeO2/Kinetics_All",
    "~/calc/Kinetics/Cluster/Pt4_CeO2/Kinetics_All",
    "~/calc/Kinetics/Cluster/Pd10_CeO2/Kinetics_All",
    "~/calc/Kinetics/Cluster/Pt10_CeO2/Kinetics_All"
]
folders = [os.path.expanduser(f) for f in folders]

# NEU: 2×3-Layout wie „perfektes“ Script
subplot_order = [
    "$Pd_{1}$", "$Pd_{4}$", "$Pd_{10}$",  # obere Reihe
    "$Pt_{1}$", "$Pt_{4}$", "$Pt_{10}$"   # untere Reihe
]

legend_labels = {
    os.path.expanduser("~/calc/Kinetics/SACs_edge/Pd1_CeO2/Kinetics_All"): "$Pd_{1}$",
    os.path.expanduser("~/calc/Kinetics/SACs_edge/Pt1_CeO2/Kinetics_All"): "$Pt_{1}$",
    os.path.expanduser("~/calc/Kinetics/Cluster/Pd4_CeO2/Kinetics_All"):   "$Pd_{4}$",
    os.path.expanduser("~/calc/Kinetics/Cluster/Pt4_CeO2/Kinetics_All"):   "$Pt_{4}$",
    os.path.expanduser("~/calc/Kinetics/Cluster/Pd10_CeO2/Kinetics_All"):  "$Pd_{10}$",
    os.path.expanduser("~/calc/Kinetics/Cluster/Pt10_CeO2/Kinetics_All"):  "$Pt_{10}$"
}

# === Coverage-Daten einlesen ===
min_temp = 350
max_temp = 1000
combined_df = pd.DataFrame()
columns_to_exclude = ["CO", "CO2", "O2", "N2"]

for folder in folders:
    cov_file = os.path.join(folder, "range", "coverage.dat")
    if os.path.exists(cov_file):
        df = pd.read_csv(cov_file, delimiter="\t")
        df.columns = [col.strip() for col in df.columns]
        df = df.apply(pd.to_numeric, errors='coerce')
        df["Source"] = legend_labels.get(folder, folder)
        combined_df = pd.concat([combined_df, df], ignore_index=True)

# === Plot erstellen: gestapelte Balkendiagramme (2×3, 50 K) ===
if not combined_df.empty:
    combined_df = combined_df[(combined_df["Temperature"] >= min_temp) & (combined_df["Temperature"] <= max_temp)]

    # ---- 50-K-Raster zwischen min/max ----
    start_T = int(np.ceil(min_temp / 50.0) * 50)
    end_T   = int(np.floor(max_temp / 50.0) * 50)
    export_temps = list(range(start_T, end_T + 1, 50))

    export_folder = "plots/cov_export"
    os.makedirs(export_folder, exist_ok=True)

    # ---------- Helpers ----------
    def sp_color(species):
        cname = reaction_colors.get(species, "grey")
        return farbe_rgb_map.get(cname, "#888888")

    def sp_label(species):
        return reaction_labels.get(species, species)

    # Sammelcontainer für Gesamt-Export + Legende
    all_lines = []
    all_header = ["Mechanism", "Species"] + [str(T) for T in export_temps]
    all_lines.append("\t".join(all_header))
    legend_handles = {}  # label -> handle (einzigartig)

    # ---------- Figure/Grid im „perfekten“ Stil ----------
    fig, axs = plt.subplots(2, 3, figsize=(11.5, 7.5), sharex=True)
    axs = axs.flatten()
    fig.subplots_adjust(hspace=0.28, wspace=0.10)

    # feste Reihenfolge (falls vorhanden)
    sources_in_df = set(combined_df["Source"].unique())
    ordered_sources = [s for s in subplot_order if s in sources_in_df]

    for i, source in enumerate(ordered_sources[:6]):
        ax = axs[i]
        sdf = combined_df[combined_df["Source"] == source].copy()
        if sdf.empty:
            ax.axis("off")
            continue

        # Spalten (nur Oberflächenspezies)
        species_cols = [c for c in sdf.columns if c not in ["Temperature", "Source"] + columns_to_exclude]
        # wie zuvor: Mini-Beiträge raus
        species_cols = [s for s in species_cols if sdf[s].abs().max() >= 0.01]
        if not species_cols:
            ax.axis("off")
            continue

        # Werte-Matrix (Temps × Species)
        vals = []
        for T in export_temps:
            row = sdf[np.isclose(sdf["Temperature"], T)]
            if row.empty:
                vals.append([0.0] * len(species_cols))  # fehlende T als 0
            else:
                vals.append([float(row[s].values[0]) for s in species_cols])
        vals = np.array(vals)

        # Species nach Gesamtsumme sortieren (absteigend) → stabilere Legende
        totals = vals.sum(axis=0)
        order = np.argsort(-totals)
        species_cols = [species_cols[j] for j in order]
        vals = vals[:, order]

        # --- Gestapelte Balken ---
        x = np.arange(len(export_temps))
        bottom = np.zeros(len(export_temps))
        for j, sp in enumerate(species_cols):
            bar = ax.bar(
                x, vals[:, j], bottom=bottom,
                label=sp_label(sp), color=sp_color(sp),
                width=0.85, linewidth=0
            )
            bottom += vals[:, j]
            # Legendeneintrag einmalig merken
            lab = sp_label(sp)
            if lab not in legend_handles:
                legend_handles[lab] = bar

        # --- Titel & Achsen wie gehabt ---
        ax.set_title(f"{source}", fontsize=12, pad=5)
        row, col = divmod(i, 3)
        if col == 0:
            ax.set_ylabel("Coverage", fontsize=12)
        else:
            ax.tick_params(labelleft=False)

        if row == 1:
            ax.set_xlabel("Temperature [K]", fontsize=12)
            ax.tick_params(labelbottom=True)
        else:
            ax.tick_params(labelbottom=False)

        ax.set_ylim(0.0, 1.02)
        ax.set_xlim(-0.6, len(export_temps) - 0.4)
        # nur Labels bei 100-K-Schritten zeigen
        tick_idx = [i for i, T in enumerate(export_temps) if T % 100 == 0]
        ax.set_xticks(tick_idx)
        ax.set_xticklabels([str(export_temps[i]) for i in tick_idx])
        ax.grid(axis='y', linestyle='--', alpha=0.5)
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.title.set_position([.5, 1.01])

        # --- Mechanismus-spezifischer Export ---
        short = source.replace(",", "").replace(" ", "_").replace(".", "")
        out_lines = []
        out_lines.append("\t".join(["Species"] + [str(T) for T in export_temps]))
        for j, sp in enumerate(species_cols):
            out_lines.append("\t".join([sp] + [f"{v:.6e}" for v in vals[:, j]]))
            # auch zur Gesamtdatei hinzufügen
            all_lines.append("\t".join([short, sp] + [f"{v:.6e}" for v in vals[:, j]]))
        with open(os.path.join(export_folder, f"{short}.txt"), "w") as f:
            f.write("\n".join(out_lines))

    # Falls weniger als 6 Quellen vorhanden: übrige Panels ausblenden
    for k in range(len(ordered_sources), 6):
        axs[k].axis("off")

    # ---------- Legende unten (eigene Achse) ----------
    # eindeutige Reihenfolge: nach Label sortieren, dann dynamische Spaltenzahl
    labels_sorted = sorted(legend_handles.keys())
    handles_sorted = [legend_handles[l] for l in labels_sorted]
    plt.tight_layout(rect=[0, 0.18, 1, 0.98])
    legend_ax = fig.add_axes([0.0, 0.0, 1.0, 0.25])
    legend_ax.axis("off")
    n_items = len(labels_sorted)
    ncol = 7 if n_items <= 28 else 8 if n_items <= 40 else 9
    legend_ax.legend(
        handles_sorted, labels_sorted,
        loc="center", ncol=ncol,
        frameon=True, fancybox=True,
        fontsize=10, handlelength=1.6,
        columnspacing=0.7, handletextpad=0.3,
        borderpad=0.2, labelspacing=0.2
    )

    fig.suptitle("Coverage of $Pd_{n}/CeO_{2}(111)$", fontsize=12, y=0.98)
    os.makedirs("plots", exist_ok=True)
    plt.savefig("plots/all_bars.png", dpi=600)
    plt.savefig("plots/all_bars.pdf")
    plt.close()
    print("Ordered combined bar plots saved as plots/all_bars.(png|pdf)")

    # ---------- Gesamtdaten speichern (50 K) ----------
    outpath_all = os.path.join("plots", "all_mechs_50K.txt")
    with open(outpath_all, "w") as f:
        f.write("\n".join(all_lines))
    print(f"Combined coverage export (50 K) written to {outpath_all}")

else:
    print("No data found.")

