#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from matplotlib.ticker import MultipleLocator

# ---------- einheitliche Schriftgrößen: alles 10 ----------
plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 10,
    "axes.labelsize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
})

folders = [
    "~/calc/Kinetics/SACs_edge/Pd1_CeO2/Kinetics_All",
    "~/calc/Kinetics/Cluster/Pd4_CeO2/Kinetics_All",
    "~/calc/Kinetics/Cluster/Pd10_CeO2/Kinetics_All",
    "~/calc/Kinetics/SACs_edge/Pt1_CeO2/Kinetics_All",
    "~/calc/Kinetics/Cluster/Pt4_CeO2/Kinetics_All",
    "~/calc/Kinetics/Cluster/Pt10_CeO2/Kinetics_All",
]
folders = [os.path.expanduser(p) for p in folders]
subplot_labels = ["$Pd_{1}$", "$Pd_{4}$", "$Pd_{10}$", "$Pt_{1}$", "$Pt_{4}$", "$Pt_{10}$"]
folder_to_label = dict(zip(folders, subplot_labels))

# Farben ohne Marker
series_style = {
    "CO": {"color": "#919191"},     # schwarz
    "O2": {"color": "#0096FF"},     # brick red (firebrick)
}

def load_orders_df(folder):
    orders_dir = os.path.join(folder, "range", "orders")
    for fname in ("orders.dat", "order.dat"):
        fpath = os.path.join(orders_dir, fname)
        if os.path.exists(fpath):
            try:
                df = pd.read_csv(fpath, sep="\t")
            except Exception:
                df = pd.read_csv(fpath, delim_whitespace=True)
            df.columns = [c.strip() for c in df.columns]
            for c in df.columns:
                df[c] = pd.to_numeric(df[c], errors='coerce')
            df = df.dropna(subset=["Temperature"]).sort_values("Temperature")
            df["Source"] = folder_to_label.get(folder, folder)
            return df
    print(f"[WARN] keine orders.dat/order.dat in: {orders_dir}")
    return None

dfs = []
for folder in folders:
    df = load_orders_df(folder)
    if df is not None:
        dfs.append(df)
if not dfs:
    print("No data found."); raise SystemExit

combined = pd.concat(dfs, ignore_index=True)
min_temp, max_temp = 350, 1000
combined = combined[(combined["Temperature"] >= min_temp) & (combined["Temperature"] <= max_temp)]

export_temps = list(range(400, 1001, 50))
os.makedirs("plots", exist_ok=True)
export_path = os.path.join("plots", "orders_all.txt")

lines = []
header = ["Mechanism", "Species"] + [str(T) for T in export_temps]
lines.append("\t".join(header))
for label in subplot_labels:
    sdf = combined[combined["Source"] == label]
    if sdf.empty: 
        continue
    for specie in ["CO", "O2"]:
        if specie not in sdf.columns: 
            continue
        row = [label, specie]
        for T in export_temps:
            hit = sdf[np.isclose(sdf["Temperature"], T)]
            row.append(f"{hit[specie].values[0]:.6e}" if not hit.empty else "NaN")
        lines.append("\t".join(row))
with open(export_path, "w") as f:
    f.write("\n".join(lines))
print(f"Combined order export written to {export_path}")

# ---------------- Plot ----------------
fig, axs = plt.subplots(2, 3, figsize=(11.5, 7.5), sharex=True, sharey=True)
axs = axs.flatten()
fig.subplots_adjust(hspace=0.20, wspace=0.10)

for i, label in enumerate(subplot_labels):
    ax = axs[i]
    sdf = combined[combined["Source"] == label].sort_values("Temperature")
    if sdf.empty:
        ax.set_visible(False); continue

    T_data = sdf["Temperature"].values
    xhi = np.linspace(max(min_temp, np.nanmin(T_data)), min(max_temp, np.nanmax(T_data)), 800)

    for specie in ["CO", "O2"]:
        if specie not in sdf.columns:
            continue
        y = sdf[specie].values
        kind = "cubic" if len(sdf) >= 4 else "linear"
        f = interp1d(T_data, y, kind=kind, fill_value="extrapolate")
        y_smooth = f(xhi)
        ax.plot(xhi, y_smooth, lw=1.8, color=series_style[specie]["color"], label=specie)

    ax.set_title(rf"{label}/CeO$_2$(111)", pad=4)

    if i in [0, 3]:
        ax.set_ylabel("Reaction Order")
    else:
        ax.tick_params(labelleft=False)

    if i in [3, 4, 5]:
        ax.set_xlabel("Temperature [K]")
        ax.tick_params(labelbottom=True)
    else:
        ax.tick_params(labelbottom=False)

    ax.set_xlim(min_temp, max_temp)
    ax.set_ylim(-0.5, 1.0)
    ax.xaxis.set_major_locator(MultipleLocator(100))
    ax.yaxis.set_major_locator(MultipleLocator(0.5))   # 0.5er sichtbar
    ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)

fig.suptitle("Reaction orders of $Pd_{n}$/$Pt_{n}$ on CeO$_2$(111)", y=0.98)

# Legende unten
handles, labels = [], []
seen = set()
for ax in axs:
    for line in ax.get_lines():
        lab = line.get_label()
        if lab and not lab.startswith("_") and lab not in seen:
            handles.append(line); labels.append(lab); seen.add(lab)

plt.tight_layout(rect=[0, 0.12, 1, 0.97])
legend_ax = fig.add_axes([0.0, 0.0, 1.0, 0.2])
legend_ax.axis("off")
legend_ax.legend(handles, labels, loc="center",
                 ncol=6, frameon=True, fancybox=True,
                 handlelength=1.8, columnspacing=0.8,
                 handletextpad=0.4, borderpad=0.2, labelspacing=0.2)

#plt.savefig("plots/orders_all.png", dpi=600)
plt.savefig("plots/orders_all.pdf")
plt.close()
print("Updated plot saved as plots/orders_all.pdf")

