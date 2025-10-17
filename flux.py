import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import FancyArrowPatch, RegularPolygon
from PIL import Image
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

# ----------------------------------------
# Farbdefinitionen für Mechanismen
pathway_colors = {
    "MvK": (1.0, 0.5, 0.0),        # Orange
    "sLH": (0.5, 0.8, 1.0),        # Hellblau
    "cLH": (0.2, 0.2, 1.0),        # Dunkelblau
    "sMix": (0.6, 0.6, 0.6),       # Grau
    "cMix": (0.2, 0.8, 0.4)        # Grün
}

# ----------------------------------------
# Aus mech.txt laden
mech_path = "flux/mech.txt"
local_vars = {}
with open(mech_path, "r") as f:
    exec(f.read(), {}, local_vars)

Zuordnung_mit_grid = local_vars["Zuordnung_mit_grid"]
mechanism_assignment = local_vars["mechanism_assignment"]

# ----------------------------------------
def read_rates(file_path, temperature):
    with open(file_path, "r") as f:
        lines = f.readlines()
    reaction_names = lines[0].strip().split('\t')[1:]
    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.strip().split('\t')
        temp = float(parts[0])
        if abs(temp - temperature) < 1e-2:
            rates = parts[1:]
            return dict(zip(reaction_names, rates))
    return None

def map_rates_to_alpha(rates_dict, min_alpha=0.2, max_alpha=1.0, min_exp=-15):
    if not rates_dict:
        return {}
    log_rates = [np.log10(max(float(r), 1e-300)) for r in rates_dict.values() if float(r) > 0]
    if not log_rates:
        return {}
    min_log = max(min(log_rates), min_exp)
    max_log = max(log_rates)
    styles = {}
    for reaction, rate in rates_dict.items():
        try:
            rate_val = float(rate)
            if rate_val <= 0:
                continue
            log_r = np.log10(max(rate_val, 1e-300))
            log_r = max(log_r, min_exp)
            norm = (log_r - min_log) / (max_log - min_log + 1e-12)
            alpha = min_alpha + norm * (max_alpha - min_alpha)
            styles[reaction] = alpha
        except:
            continue
    return styles

def draw_arrows(ax, rates_dict, grid_mapping, shift_distance=0.1):
    alpha_values = map_rates_to_alpha(rates_dict)
    for reaction, rate_str in rates_dict.items():
        if '->' in reaction:
            forward = True
            lhs, rhs = reaction.split('->')
        elif '<-' in reaction:
            forward = False
            rhs, lhs = reaction.split('<-')
        else:
            continue
        lhs_species = [s.strip() for s in lhs.strip().split('+')]
        rhs_species = [s.strip() for s in rhs.strip().split('+')]
        lhs_mapped = next((s for s in lhs_species if s in grid_mapping), None)
        rhs_mapped = next((s for s in rhs_species if s in grid_mapping), None)
        start = lhs_mapped if forward else rhs_mapped
        end = rhs_mapped if forward else lhs_mapped
        if not start or not end:
            continue

        y0, x0 = grid_mapping[start]["grid"]
        y1, x1 = grid_mapping[end]["grid"]
        dx, dy = x1 - x0, y1 - y0
        length = np.sqrt(dx**2 + dy**2)
        if length == 0:
            continue
        ux, uy = dx / length, dy / length
        perp_dx, perp_dy = -uy, ux
        shift = shift_distance if forward else -shift_distance
        x0a = x0 + perp_dx * shift
        y0a = y0 + perp_dy * shift
        x1a = x1 + perp_dx * shift
        y1a = y1 + perp_dy * shift

        alpha = alpha_values.get(reaction, 0.6)
        mechanism = mechanism_assignment.get(reaction, None)
        base_color = pathway_colors.get(mechanism, (1.0, 0.4, 0.7))
        color_with_alpha = (*base_color, alpha)

        # Pfeil zeichnen
        ax.add_patch(FancyArrowPatch((x0a, y0a), (x1a, y1a), arrowstyle='-', linewidth=15, color=color_with_alpha, zorder=2))
        print(f"{lhs.strip():<40} {'→' if forward else '←'} {rhs.strip():<40}  Rate: {float(rate_str):.3e}")

        # Text-Position + Rotation
        tx = (x0a + x1a) / 2
        ty = (y0a + y1a) / 2
        dx_vis = x1a - x0a
        dy_vis = -(y1a - y0a)  # invertierte y-Achse
        angle = np.degrees(np.arctan2(dy_vis, dx_vis))
        if angle > 90 or angle < -90:
            angle += 180

        ax.text(tx, ty, f"{float(rate_str):.1e}", fontsize=10, ha='center', va='center', rotation=angle,
                rotation_mode='anchor', color="black", alpha=1.0, zorder=100)

        # kleine Pferdchen – angepasst für CO2-Fälle
        co2_case = any(s.strip() == "CO2" for s in lhs_species + rhs_species)
        show_pferdchen = (forward and not co2_case) or (not forward and co2_case)
        if show_pferdchen:
            dxn, dyn = dx / length, dy / length
            perp_dx, perp_dy = -dyn, dxn
            offset_dir = -1 if (not forward and co2_case) else 1  # CO2-Fälle: Pferdchen unten
            center_x = tx + perp_dx * 0.12 * offset_dir
            center_y = ty + perp_dy * 0.12 * offset_dir
            for i in [-1, 0, 1]:
                offset_x = center_x + i * 0.06 * dxn
                offset_y = center_y + i * 0.06 * dyn
                # Orientierung drehen, wenn CO2-Edukt → Rückreaktion
                if not forward and co2_case:
                    orientation = np.arctan2(dy, dx) + np.pi / 2
                else:
                    orientation = np.arctan2(dy, dx) - np.pi / 2
                triangle = RegularPolygon((offset_x, offset_y), numVertices=3, radius=0.05,
                                          orientation=orientation,
                                          color=color_with_alpha, zorder=101)
                ax.add_patch(triangle)




# ----------------------------------------
fig, ax = plt.subplots(figsize=(9, 10))
ax.set_xlim(0, 5)
ax.set_ylim(0, 6)
ax.invert_yaxis()
ax.grid(False)

image_folder = "flux"
image_size = 0.55
gas_phase_species = {"O2", "CO", "CO2"}

for species, props in Zuordnung_mit_grid.items():
    if species in gas_phase_species:
        continue
    row, col = props["grid"]
    image_path = os.path.join(image_folder, props["image"])
    if os.path.exists(image_path):
        img = mpimg.imread(image_path)
        left = col - image_size / 2
        right = col + image_size / 2
        bottom = row - image_size / 2
        top = row + image_size / 2
        ax.imshow(img, extent=(left, right, top, bottom), zorder=10)
        image_height = abs(top - bottom)
        quarter_height = image_height / 4
        balken_y = bottom
        ax.add_patch(plt.Rectangle((left, balken_y), right - left, quarter_height,
                                   facecolor=(0, 0, 0, 0.5), edgecolor='none', zorder=11))
        label_text = props.get("label", species)
        ax.text((left + right) / 2, balken_y + quarter_height / 2, label_text,
                color='white', alpha=1.0, fontweight='bold', fontsize=9,
                ha='center', va='center', zorder=100)
    else:
        ax.text(col, row, props["image"], ha='center', va='center', fontsize=9,
                bbox=dict(boxstyle="round", fc="salmon"), zorder=12)
T = 750
rates_dat_path = "Kinetics_All/range/rates.dat"
output_folder = "flux"
output_file = f"rates_{T}.txt"
output_path = os.path.join(output_folder, output_file)

rates_dict = read_rates(rates_dat_path, T)
if rates_dict:
    os.makedirs(output_folder, exist_ok=True)
    with open(output_path, "w") as f:
        for reaction, rate in rates_dict.items():
            f.write(f"{reaction}: {rate} (0.0,0.0)\n")
    draw_arrows(ax, rates_dict, Zuordnung_mit_grid)
else:
    print(f"⚠️ Keine Raten für T={T} K gefunden in {rates_dat_path}")
ax.set_aspect('equal')
plt.tight_layout()
plt.show()
