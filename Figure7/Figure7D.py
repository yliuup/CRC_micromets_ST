import os
import gzip
import cv2
import imageio
import numpy as np
import pandas as pd
from collections import defaultdict
import skimage.draw as ski
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap, Normalize

# ===================== Paths and configuration =====================

BASE_CRC = "./plot"
BASE_SEG = "./visium_HD"
SAMPLES = ["ML4"]

BLUE_WHITE_RED = LinearSegmentedColormap.from_list(
    "blue_white_red", ["#4DBBD5", "white", "#E64B35"]
)
NORM = Normalize(vmin=-0.3, vmax=0.7, clip=True)

MY_COLS = [
    "#984EA3", "#E41A1C", "#6A3D9A", "#FF7F00", "#33A02C", "#A6CEE3", "#1F78B4",
    "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6", "#FFFF99", "#B15928",
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4",
    "#91D1C2", "#7E6148", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
    "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
    "#CCEBC5", "#FFED6F", "#FFB3B3", "#99CCCC", "#FFCC99"
]

os.chdir(BASE_CRC)

# ===================== Utility functions =====================

def hex_to_rgb(hex_color):
    hex_color = hex_color.lstrip("#")
    return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))


def fast_find_cell_boundaries(mask):
    coord_map = defaultdict(list)
    it = np.nditer(mask, flags=['multi_index'])

    for val in it:
        cid = int(val)
        if cid != 0:
            coord_map[cid].append(it.multi_index)

    boundaries = {}
    for cid, coords in coord_map.items():
        coords = np.array(coords)
        y_min, x_min = coords.min(axis=0)
        y_max, x_max = coords.max(axis=0)
        roi = np.zeros((y_max - y_min + 1, x_max - x_min + 1), dtype=np.uint8)
        roi_coords = coords - [y_min, x_min]
        roi[roi_coords[:, 0], roi_coords[:, 1]] = 1
        contours, _ = cv2.findContours(roi, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        if not contours:
            continue
        boundary = max(contours, key=cv2.contourArea)[:, 0, :]
        boundary[:, 0] += x_min
        boundary[:, 1] += y_min
        boundaries[cid] = boundary
    return boundaries


def plot_fill_and_boundary_on_image(
    img,
    cell_boundary,
    border_color_for_cid,
    fill_value_for_cid,
    cmap,
    norm,
    thick=2,
    alpha=0.8
):
    draw = img.copy().astype(np.float32)
    H, W = draw.shape[:2]

    for cid, poly in cell_boundary.items():
        if poly.shape[0] < 3:
            continue
        x = np.clip(poly[:, 0], 0, W - 1)
        y = np.clip(poly[:, 1], 0, H - 1)
        rr, cc = ski.polygon(y, x, shape=(H, W))
        val = fill_value_for_cid.get(cid, np.nan)
        if np.isnan(val):
            continue
        rgba = cmap(norm(val))
        fill_rgb = np.array([255 * c for c in rgba[:3]], dtype=np.float32)
        draw[rr, cc] = (1 - alpha) * draw[rr, cc] + alpha * fill_rgb

    for cid, poly in cell_boundary.items():
        if poly.shape[0] < 3:
            continue
        x = np.clip(poly[:, 0], 0, W - 1)
        y = np.clip(poly[:, 1], 0, H - 1)
        rr, cc = ski.polygon_perimeter(y, x, shape=(H, W), clip=True)
        border_hex = border_color_for_cid.get(cid, "#B0B0B0")
        border_rgb = np.array(hex_to_rgb(border_hex), dtype=np.float32)
        for yy, xx in zip(rr, cc):
            rr_disk, cc_disk = ski.disk((yy, xx), radius=thick, shape=(H, W))
            draw[rr_disk, cc_disk] = border_rgb

    return np.clip(draw, 0, 255).astype(np.uint8)


def draw_legend_discrete(name_to_hex, out_path, title="tum"):
    fig, ax = plt.subplots(figsize=(4, max(1.5, len(name_to_hex) * 0.25)))
    handles = [mpatches.Patch(color=v, label=k) for k, v in name_to_hex.items()]
    ax.legend(handles=handles, loc='center left', frameon=False, title=title)
    ax.axis('off')
    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)


def save_with_colorbar(img_array, cmap, norm, out_path, label="signature2_1"):
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(img_array)
    ax.axis("off")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(label, fontsize=12)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def process_sample(sample):
    print(f"Processing {sample}...")

    mask = np.load(os.path.join(BASE_SEG, sample, "segmentation_stardist.npy"))
    with gzip.open(os.path.join(BASE_SEG, sample, "image_to_segmentation.npy.gz"), "rb") as f:
        he_image = np.load(f, allow_pickle=True)
    if he_image.ndim == 2:
        he_image = cv2.cvtColor(he_image, cv2.COLOR_GRAY2RGB)

    meta = pd.read_csv(os.path.join(BASE_CRC, f"{sample}_meta.csv"))
    meta["bar"] = meta["bar"].astype(str)

    mask_ids = np.unique(mask)
    mask_ids = mask_ids[mask_ids != 0].astype(str)
    meta = meta[meta["bar"].isin(mask_ids)].copy()

    tum_levels = list(pd.unique(meta["tum"].astype(str)))
    tum_to_color = {t: MY_COLS[i % len(MY_COLS)] for i, t in enumerate(tum_levels)}

    border_color_for_cid = {int(row["bar"]): tum_to_color[str(row["tum"])] for _, row in meta.iterrows()}
    fill_value_for_cid = {int(row["bar"]): float(row["signature2_1"]) for _, row in meta.iterrows()}

    cell_boundary = fast_find_cell_boundaries(mask)
    keep_cids = set(int(x) for x in meta["bar"])
    cell_boundary = {cid: b for cid, b in cell_boundary.items() if cid in keep_cids}

    overlay_img = plot_fill_and_boundary_on_image(
        he_image, cell_boundary,
        border_color_for_cid,
        fill_value_for_cid,
        BLUE_WHITE_RED, NORM
    )

    imageio.imwrite(f"{BASE_CRC}/{sample}_overlay.jpg", overlay_img, quality=95)
    save_with_colorbar(overlay_img, BLUE_WHITE_RED, NORM, f"{BASE_CRC}/{sample}_overlay_cbar.jpg")
    draw_legend_discrete(tum_to_color, f"{BASE_CRC}/{sample}_legend.jpg")


if __name__ == "__main__":
    for sample in SAMPLES:
        process_sample(sample)
