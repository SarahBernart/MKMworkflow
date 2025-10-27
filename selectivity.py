#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CO2 production split by pathway (no markers, 2x3 layout, example style):
- CO2_LH1, CO2_LH2, CO2_MvK1
- CO2_MvK2_int + CO2_MvK2_nn (if present), otherwise CO2_MvK2

Rohdaten: keine Interpolation/Glättung. Plot zeigt log10(Rate) auf linearer y-Achse.
Horizontale Linie bei y=0. Einheitliche y-Limits für alle Subplots.

Outputs:
- plots/CO2_prod_all.png
- plots/CO2_prod_all.txt  (Rohwerte bei 50-K-Schritten; NaN falls T nicht exakt vorhanden)
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =========================
# Settings
# =========================
MECHANISM_FILE = "mech.txt"  # optional: farbe_rgb_map, reaction_colors, reaction_labels, species_* ...

FOLDERS = [
    os.path.expanduser("~/calc/Kinetics/SACs_edge/Pd1_CeO2/Kinetics_All_dsc"),
    os.path.expanduser("~/calc/Kinetics/SACs_edge/Pt1_CeO2/Kinetics_All_dsc"),
    os.path.expanduser("~/calc/Kinetics/Cluster/Pd4_CeO2/Kinetics_All_dsc"),
    os.path.expanduser("~/calc/Kinetics/Cluster/Pt4_CeO2/Kinetics_All_dsc"),
    os.path.expanduser("~/calc/Kinetics/Cluster/Pd10_CeO2/Kinetics_All_dsc"),
    os.path.expanduser("~/calc/Kinetics/Cluster/Pt10_CeO2/Kinetics_All_dsc"),
]

T_MIN, T_MAX = 350, 1000
EXPORT_TEMPS = list(range(400, 1101, 50))  # Rohwerte; NaN, wenn T nicht exakt in Datei

# Spalten
BASE_COLS   = ["CO2_LH1", "CO2_LH2", "CO2_MvK1"]
MVK2_SINGLE = "CO2_MvK2"
MVK2_INT    = "CO2_MvK2_int"
MVK2_NN     = "CO2_MvK2_nn"
CANDIDATE_COLS = BASE_COLS + [MVK2_SINGLE, MVK2_INT, MVK2_NN]

# derivatives.dat – robuste Suchpfade
DERIVATIVE_CANDIDATES = [
    "range/deriv/derivatives.dat",
    "range/derivatives/derivatives.dat",
    "range/dsc/derivatives.dat",
    "range/drc/derivatives.dat",
    "range/derivative.dat",
    "range/derivatives.dat",
    "deriv/derivatives.dat",
    "derivatives/derivatives.dat",
    "dsc/derivatives.dat",
    "drc/derivatives.dat",
    "derivative.dat",
    "derivatives.dat",
]

# gewünschte Anordnung (2x3): oben Pd1, Pd4, Pd10; unten Pt1, Pt4, Pt10
ORDER_TOP    = ["$Pd_{1}$", "$Pd_{4}$", "$Pd_{10}$"]
ORDER_BOTTOM = ["$Pt_{1}$", "$Pt_{4}$", "$Pt_{10}$"]

# =========================
# (optionale) Styles aus mech.txt
# =========================
farbe_rgb_map = {}
reaction_colors = {}
reaction_labels = {}
species_colors = {}
species_labels = {}
species_linestyles = {}
reaction_linestyles = {}

# Reihenfolge (global) – letzter Eintrag liegt "oben"
GLOBAL_PLOT_ORDER = ["CO2_LH1","CO2_LH2","CO2_MvK1","CO2_MvK2_nn","CO2_MvK2_int","CO2_MvK2"]
PLOT_ORDER_PER_SOURCE = {}  # z.B. {"$Pd_{1}$": ["CO2_MvK1", ...]}

if os.path.exists(MECHANISM_FILE):
    try:
        with open(MECHANISM_FILE, "r", encoding="utf-8") as f:
            code = f.read()
        local_vars = {}
        exec(code, {}, local_vars)
        farbe_rgb_map           = local_vars.get("farbe_rgb_map", {})
        reaction_colors         = local_vars.get("reaction_colors", {})
        reaction_labels         = local_vars.get("reaction_labels", {})
        species_colors          = local_vars.get("species_colors", {})
        species_labels          = local_vars.get("species_labels", {})
        species_linestyles      = local_vars.get("species_linestyles", {})
        reaction_linestyles     = local_vars.get("reaction_linestyles", {})
        GLOBAL_PLOT_ORDER       = local_vars.get("plot_order", GLOBAL_PLOT_ORDER)
        PLOT_ORDER_PER_SOURCE   = local_vars.get("plot_order_per_source", PLOT_ORDER_PER_SOURCE)
    except Exception as e:
        print(f"[warn] Could not parse mech.txt ({e}). Using defaults.")

# --- Palette (hat Vorrang vor mech.txt) ---
USER_COLOR_MAP = {
    "CO2_LH1":      "#66C2FF",
    "CO2_LH2":      "#9ACD32",
    "CO2_MvK1":     "#7F7F7F",
    "CO2_MvK2_int": "#FFA032",
    "CO2_MvK2_nn":  "#C00000",
    "CO2_MvK2":     "#FFA032",  # fallback falls nur kombiniert
}

# --- Linienstile (gültige Matplotlib-Stile) ---
USER_LINESTYLE_MAP = {
    "CO2_MvK1":     "-.",  # .-
    "CO2_MvK2_nn":  ":",   # ..
    "CO2_LH2":      "-.",  # .-
    "CO2_MvK2_int": ":",   # konsistent
    # andere bleiben solid
}

DEFAULT_LABELS = {
    "CO2_LH1":      r"CO$_2$ (LH1)",
    "CO2_LH2":      r"CO$_2$ (LH2)",
    "CO2_MvK1":     r"CO$_2$ (MvK1)",
    "CO2_MvK2":     r"CO$_2$ (MvK2)",
    "CO2_MvK2_int": r"CO$_2$ (MvK2-int)",
    "CO2_MvK2_nn":  r"CO$_2$ (MvK2-nn)",
}
COLOR_KEY_MAP = {
    "CO2_LH1":       "LH1(CO2des)",
    "CO2_LH2":       "LH2(CO2des)",
    "CO2_MvK1":      "MvK1(CO2des)",
    "CO2_MvK2":      "MvK2(CO2des)",
    "CO2_MvK2_int":  "MvK2_int(CO2des)",
    "CO2_MvK2_nn":   "MvK2_nn(CO2des)",
}

# =========================
# Helpers
# =========================
def infer_label_from_path(path: str) -> str:
    m = re.search(r"P([dt])(\d+)_CeO2", path, re.IGNORECASE)
    if m:
        letter = m.group(1).lower()
        num = m.group(2)
        metal = "Pd" if letter == "d" else "Pt"
        return f"${metal}_{{{num}}}$"
    return os.path.basename(path.rstrip("/"))

def find_derivatives_path(folder: str) -> str | None:
    for rel in DERIVATIVE_CANDIDATES:
        p = os.path.join(folder, rel)
        if os.path.exists(p):
            return p
    return None

def determine_cols_for_source(columns: list[str]) -> list[str]:
    cols = [c for c in BASE_COLS if c in columns]
    has_int = MVK2_INT in columns
    has_nn  = MVK2_NN  in columns
    has_single = MVK2_SINGLE in columns
    if has_int or has_nn:
        if has_int: cols.append(MVK2_INT)
        if has_nn:  cols.append(MVK2_NN)
    elif has_single:
        cols.append(MVK2_SINGLE)
    return cols

def get_color_for(colname: str) -> str | None:
    if colname in USER_COLOR_MAP:
        return USER_COLOR_MAP[colname]
    colorname = species_colors.get(colname)
    if not colorname:
        rxn_key = COLOR_KEY_MAP.get(colname)
        if rxn_key:
            colorname = reaction_colors.get(rxn_key)
    if colorname:
        return farbe_rgb_map.get(colorname, None)
    return None

def get_linestyle_for(colname: str) -> str | tuple | None:
    if colname in USER_LINESTYLE_MAP:
        return USER_LINESTYLE_MAP[colname]
    ls = species_linestyles.get(colname)
    if not ls:
        rxn_key = COLOR_KEY_MAP.get(colname)
        if rxn_key:
            ls = reaction_linestyles.get(rxn_key)
    return ls  # None => default '-'

def get_label_for(colname: str) -> str:
    if colname in species_labels:
        return species_labels[colname]
    rxn_key = COLOR_KEY_MAP.get(colname)
    if rxn_key and rxn_key in reaction_labels:
        return reaction_labels[rxn_key]
    return DEFAULT_LABELS.get(colname, colname)

def order_cols_for_source(source: str, cols: list[str]) -> list[str]:
    custom = PLOT_ORDER_PER_SOURCE.get(source, [])
    global_order = GLOBAL_PLOT_ORDER or []
    if not custom and not global_order:
        return cols
    rank = {}
    for i, name in enumerate(custom):
        rank[name] = i
    offset = len(rank)
    for j, name in enumerate(global_order):
        if name not in rank:
            rank[name] = offset + j
    base = {name: k for k, name in enumerate(cols)}
    return sorted(cols, key=lambda n: (rank.get(n, 10_000), base.get(n, 0)))

# ---------- Robust column normalization & aliases ----------
def _norm_colname(c: str) -> str:
    c = c.strip()
    c = c.replace("₂", "2")                       # CO₂ -> CO2
    c = c.replace("-", "_").replace(" ", "_")
    c = re.sub(r"[()]", "_", c)                   # remove brackets
    c = re.sub(r"_+", "_", c)
    return c

ALIASES = [
    (r'(?i)^co2_mvk2(?:_(?:tot|total|comb(?:ined)?)?)?$', "CO2_MvK2"),
    (r'(?i)^co2_mvk2_?int(?:er)?$',                       "CO2_MvK2_int"),
    (r'(?i)^co2_mvk2_?nn$',                               "CO2_MvK2_nn"),
    (r'(?i)^co2_mvk1$',                                   "CO2_MvK1"),
    (r'(?i)^co2_lh1$',                                    "CO2_LH1"),
    (r'(?i)^co2_lh2$',                                    "CO2_LH2"),
]

SAC_SOURCES = {"$Pd_{1}$", "$Pt_{1}$"}  # SACs_edge

# =========================
# Load (Rohdaten, keine Interpolation)
# =========================
frames = []
for folder in FOLDERS:
    path = find_derivatives_path(folder)
    if not path:
        print(f"[warn] No derivatives.dat under: {folder}")
        continue
    try:
        src_label = infer_label_from_path(folder)

        df = pd.read_csv(path, sep=None, engine="python")
        df.columns = [_norm_colname(c) for c in df.columns]  # normalize

        # Aliase anwenden → Standardnamen
        for patt, std in ALIASES:
            for col in list(df.columns):
                if re.match(patt, col):
                    if col != std:
                        df.rename(columns={col: std}, inplace=True)

        # SAC-Sonderfall: falls nur kombiniertes MvK2 vorhanden, als "int" behandeln
        if src_label in SAC_SOURCES:
            if (MVK2_INT not in df.columns) and (MVK2_SINGLE in df.columns):
                df[MVK2_INT] = df[MVK2_SINGLE]  # nicht umbenennen, damit Export/Logik robust bleibt

        # Temperature normalisieren
        if "Temperature" not in df.columns:
            for alt in ["Temp", "T", "temperature", "temp"]:
                if alt in df.columns:
                    df.rename(columns={alt: "Temperature"}, inplace=True)
                    break
        if "Temperature" not in df.columns:
            print(f"[warn] Missing 'Temperature' in {path}")
            continue

        keep = ["Temperature"] + [c for c in CANDIDATE_COLS if c in df.columns]
        mv_missing = [m for m in [MVK2_INT, MVK2_NN, MVK2_SINGLE] if m not in df.columns]
        if mv_missing:
            print(f"[info] {src_label} missing: {', '.join(mv_missing)}")

        df = df[keep].apply(pd.to_numeric, errors="coerce")
        df = df[(df["Temperature"] >= T_MIN) & (df["Temperature"] <= T_MAX)].copy()
        df["Source"] = src_label

        print(f"[ok] {src_label}: {path} | rows={len(df)} | cols={list(df.columns)}")
        frames.append(df)
    except Exception as e:
        print(f"[warn] Failed to read {path}: {e}")

if not frames:
    print("No usable data found.")
    raise SystemExit(0)

all_df = pd.concat(frames, ignore_index=True)

# =========================
# y-Limits auf Basis von log10(Rohdaten)
# =========================
log_vals = []
for col in CANDIDATE_COLS:
    if col in all_df.columns:
        v = all_df[col].to_numpy()
        m = np.isfinite(v) & (v > 0)
        if np.any(m):
            log_vals.append(np.log10(v[m]))

if log_vals:
    all_log = np.concatenate(log_vals)
    y_min = float(np.nanmin(all_log))
    y_max = float(np.nanmax(all_log))
    ypad = max(1.0, 0.05 * (y_max - y_min))
    y_min_global = y_min - ypad
    y_max_global = y_max + ypad
else:
    y_min_global, y_max_global = -25.0, 15.0  # Fallback

# =========================
# Export: Selectivity (Anteile) je Mechanismus, 50-K-Raster
# =========================
os.makedirs("plots", exist_ok=True)

def sel_at_T(row, cols):
    # Numerik robust: nur finite, negatives zu 0
    vals = []
    for c in cols:
        v = float(row[c]) if (c in row and np.isfinite(row[c])) else 0.0
        if v < 0: v = 0.0
        vals.append(v)
    denom = np.sum(vals)
    if denom <= 0:
        return [np.nan] * len(vals)
    return [v / denom for v in vals]

lines = []
header = ["Mechanism", "Series"] + [str(T) for T in EXPORT_TEMPS]
lines.append("\t".join(header))

for source in all_df["Source"].unique():
    sub = all_df[all_df["Source"] == source].copy()
    plot_cols = determine_cols_for_source(list(sub.columns))
    plot_cols = order_cols_for_source(source, plot_cols)
    for col in plot_cols:
        row_out = [source, col]
        for T in EXPORT_TEMPS:
            sel = sub[np.isclose(sub["Temperature"], T)]
            if sel.empty:
                row_out.append("NaN")
                continue
            r = sel.iloc[0]
            cols_here = [c for c in plot_cols if c in sel.columns]
            shares = sel_at_T(r, cols_here)
            # zugehörigen Index der gesuchten Spalte holen
            try:
                idx = cols_here.index(col)
                row_out.append(f"{shares[idx]:.6e}")
            except ValueError:
                row_out.append("NaN")
        lines.append("\t".join(row_out))

with open(os.path.join("plots", "CO2_sel_all.txt"), "w") as f:
    f.write("\n".join(lines))
print("Wrote plots/CO2_sel_all.txt")

# =========================
# Plot (2×3): Selectivity (0..1, keine Log-Skala)
# =========================
plt.rcParams['font.size'] = 12
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['legend.fontsize'] = 12

fig, axs = plt.subplots(2, 3, figsize=(11.5, 7.5), sharex=True)
axs = axs.reshape(2, 3)

sources = list(all_df["Source"].unique())
present_top = [s for s in ORDER_TOP if s in sources]
present_bot = [s for s in ORDER_BOTTOM if s in sources]
rest = [s for s in sources if s not in present_top + present_bot]
for s in rest:
    if len(present_top) < 3: present_top.append(s)
    elif len(present_bot) < 3: present_bot.append(s)
grid_order = [present_top, present_bot]

# XTicks nur alle 100 K
xt0 = int(np.ceil(T_MIN / 100.0) * 100)
xt1 = int(np.floor(T_MAX / 100.0) * 100)
xticks_100 = list(range(xt0, xt1 + 1, 100))

def make_title(src: str) -> str:
    return rf"{src}/CeO$_{{2}}$(111)"

for r, row_sources in enumerate(grid_order):
    for c, source in enumerate(row_sources[:3]):
        ax = axs[r, c]
        sub = all_df[all_df["Source"] == source].sort_values("Temperature").copy()
        if sub.empty:
            ax.axis("off"); continue
        ax.set_title(make_title(source))

        plot_cols = determine_cols_for_source(list(sub.columns))
        plot_cols = order_cols_for_source(source, plot_cols)

        # Denominator pro Temperatur (robust)
        denom = np.zeros(len(sub), dtype=float)
        numers = {}
        for col in plot_cols:
            y = sub[col].to_numpy(dtype=float)
            y[~np.isfinite(y)] = 0.0
            y[y < 0.0] = 0.0
            numers[col] = y
            denom += y
        denom[denom <= 0.0] = np.nan  # vermeidet Division durch 0

        # zeichnen (kein Log)
        T = sub["Temperature"].to_numpy()
        for z, col in enumerate(plot_cols, start=1):
            sel = numers[col] / denom
            color = get_color_for(col)
            label = get_label_for(col)
            linestyle = get_linestyle_for(col) or "-"
            if color is not None:
                ax.plot(T, sel, linewidth=2.0, color=color, linestyle=linestyle, label=label, zorder=z)
            else:
                ax.plot(T, sel, linewidth=2.0, linestyle=linestyle, label=label, zorder=z)

        ax.set_xlim(T_MIN, T_MAX)
        from matplotlib.ticker import MultipleLocator, FixedLocator, FormatStrFormatter, NullFormatter

        ax.set_ylim(0.0, 1.02)
        ax.yaxis.set_major_locator(FixedLocator([0.0, 0.5, 1.0]))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        ax.yaxis.set_minor_locator(MultipleLocator(0.10))
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.grid(which='major', axis='y', color='0.6', linestyle='-',  linewidth=0.6, alpha=0.9)
        ax.grid(which='minor', axis='y', color='0.8', linestyle='--', linewidth=0.4, alpha=0.8)
        ax.grid(linewidth=0.45, alpha=0.6)
        ax.axhline(0.0, color="black", linewidth=0.75, linestyle="--")

        ax.set_xticks(xticks_100)
        ax.set_xticklabels([str(t) for t in xticks_100])

        if c == 0:
            ax.set_ylabel("Selectivity [-]")
        else:
            ax.tick_params(labelleft=False)

        if r == 1:
            ax.set_xlabel("Temperature [K]")
        else:
            ax.tick_params(labelbottom=False)

# gemeinsame Legende unten
handles, labels = [], []
used = set()
for ax in axs.flatten():
    for ln in ax.get_lines():
        lb = ln.get_label()
        if lb and not lb.startswith("_") and lb not in used:
            handles.append(ln); labels.append(lb); used.add(lb)

plt.tight_layout(rect=[0.02, 0.15, 1.00, 1.00])
legend_ax = fig.add_axes([0.00, 0.00, 1.00, 0.22])
legend_ax.axis("off")
legend_ax.legend(
    handles, labels,
    loc="center",
    ncol=min(6, len(labels)),
    frameon=True, fancybox=True,
    fontsize=12,
    handlelength=1.8,
    columnspacing=0.8,
    handletextpad=0.3,
    borderpad=0.4,
    labelspacing=0.3
)

plt.savefig(os.path.join("plots", "CO2_sel_all.png"), dpi=600)
plt.savefig(os.path.join("plots", "CO2_sel_all.pdf"))
plt.close()
print("Saved plots/CO2_sel_all.(png|pdf)")

