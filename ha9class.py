#!/usr/bin/env python3
"""
ha9class.py
Simple H/Alpha 9-zone classifier for polarimetric SAR (compact, practical).

Usage:
  # GDAL-readable inputs -> writes a single paletted GeoTIFF 'ha9.tif' (default)
  python ha9class.py --h_file H.tif --alpha_file alpha.tif --out ha9.tif

  # raw binaries (need to define width and height)
  python ha9class.py --h_file H.bin --alpha_file alpha.bin --width 2048 --height 2048 --dtype float32 --out ha9.tif

Outputs:
  - single-band paletted GeoTIFF (byte): values 0..9 (0 = nodata/unclassified, 1..9 classes)
  - quicklook PNG

Thresholds:
  - H thresholds: (default: 0.5, 0.9)
    --h_thresh H_LOW H_HIGH

  - Alpha thresholds (degrees):
    These are H-dependent, allowing for the curved zone boundaries.

    For Low H (H <= H_LOW):
    --alpha_thresh_low_h A_L1 A_H1  (default: 42.5, 52.5)

    For Med H (H_LOW < H <= H_HIGH):
    --alpha_thresh_med_h A_L2 A_H2  (default: 40.0, 50.0)

    For High H (H > H_HIGH):
    --alpha_thresh_high_h A_L3 A_H3 (default: 45.0, 55.0)

S. Cloude and E. Pottier: "An entropy based classification scheme for land applications of polarimetric SAR". 
Geoscience and Remote Sensing, IEEE Transactions on, vol. 35, no. 1, Jan. 1997, pp. 68 - 78.
"""
import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
from osgeo import gdal

# ------------------ IO helpers ------------------
def read_with_gdal(path):
    """Read a raster file using GDAL and return array with metadata."""
    ds = gdal.Open(path, gdal.GA_ReadOnly)
    if ds is None:
        return None
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray().astype(np.float32)
    return {
        'array': arr,
        'geo': ds.GetGeoTransform(),
        'proj': ds.GetProjection(),
        'nodata': band.GetNoDataValue(),
    }

def read_raw_binary(path, width, height, dtype='float32', byteorder='little'):
    """Read a raw binary file and reshape to 2D array."""
    dt = np.dtype(dtype)
    dt = dt.newbyteorder('<' if byteorder == 'little' else '>')
    count = width * height
    a = np.fromfile(path, dtype=dt, count=count)
    if a.size != count:
        msg = f"Raw file size mismatch: expected {count} elements, got {a.size}"
        raise ValueError(msg)
    return a.reshape((height, width))

def ensure_loaded(path, width, height, dtype, byteorder):
    """Load data from GDAL-readable file or raw binary."""
    g = read_with_gdal(path)
    if g is not None:
        return g['array'], g.get('geo'), g.get('proj'), g.get('nodata')
    if width is None or height is None:
        msg = "width & height required for raw binary inputs"
        raise ValueError(msg)
    arr = read_raw_binary(path, width, height, dtype=dtype, byteorder=byteorder)
    return arr, None, None, None

# ------------------ H/alpha helpers ------------------
def autodetect_alpha_units(alpha_arr):
    """Detect if alpha is in radians or degrees based on max value."""
    amax = np.nanmax(alpha_arr)
    # radians if values <= ~1.7 (pi/2)
    return 'radians' if amax <= 1.7 else 'degrees'


def convert_alpha_to_degrees(alpha_arr):
    """Convert alpha to degrees if it's in radians."""
    units = autodetect_alpha_units(alpha_arr)
    return np.degrees(alpha_arr) if units == 'radians' else alpha_arr

# ------------------ explicit H/alpha -> class mapping ------------------
def classify_9zones(
    H,
    alpha_deg,
    h_thresh=(0.5, 0.9),
    alpha_thresh_low_h=(42.5, 52.5),
    alpha_thresh_med_h=(40.0, 50.0),
    alpha_thresh_high_h=(45.0, 55.0),
    nodata_mask=None,
):
    """
    Classify using H-dependent alpha thresholds, applying the
    user-specified (scrambled) class mapping.
    This is fully vectorized and efficient.
    Returns uint8 array with 0..9 (0 = nodata/unclassified).
    """
    h_lo, h_hi = h_thresh
    a_l1, a_h1 = alpha_thresh_low_h  # Thresholds for Low H
    a_l2, a_h2 = alpha_thresh_med_h  # Thresholds for Med H
    a_l3, a_h3 = alpha_thresh_high_h  # Thresholds for High H

    cls = np.zeros(H.shape, dtype=np.uint8)  # 0 = nodata/unclassified

    # Create masks for the 3 Entropy bins
    mask_h_low = (H <= h_lo)
    mask_h_med = (H > h_lo) & (H <= h_hi)
    mask_h_high = (H > h_hi)

    # --- Apply mapping for Low H bin ---
    # (0,0) -> Class 7
    cls[mask_h_low & (alpha_deg <= a_l1)] = 7
    # (0,1) -> Class 4
    cls[mask_h_low & (alpha_deg > a_l1) & (alpha_deg <= a_h1)] = 4
    # (0,2) -> Class 1
    cls[mask_h_low & (alpha_deg > a_h1)] = 1

    # --- Apply mapping for Med H bin ---
    # (1,0) -> Class 8
    cls[mask_h_med & (alpha_deg <= a_l2)] = 8
    # (1,1) -> Class 5
    cls[mask_h_med & (alpha_deg > a_l2) & (alpha_deg <= a_h2)] = 5
    # (1,2) -> Class 2
    cls[mask_h_med & (alpha_deg > a_h2)] = 2

    # --- Apply mapping for High H bin ---
    # (2,0) -> Class 9
    cls[mask_h_high & (alpha_deg <= a_l3)] = 9
    # (2,1) -> Class 6
    cls[mask_h_high & (alpha_deg > a_l3) & (alpha_deg <= a_h3)] = 6
    # (2,2) -> Class 3
    cls[mask_h_high & (alpha_deg > a_h3)] = 3

    # Apply nodata mask
    if nodata_mask is None:
        nodata_mask = np.isnan(H) | np.isnan(alpha_deg)

    cls[nodata_mask] = 0

    return cls

# ------------------ Color table ------------------
def get_esa_colortable():
    """
    ESA mapping, modified for standard Z1-Z9 numbering (indices 0..9).
    0 = nodata (transparent).
    """
    return {
        0: (0, 0, 0, 0),        # 0 = nodata (transparent)
        1: (255,165,  0,255),   # 1 dihedral -> orange
        2: (204,  0,  0,255),   # 2 forestry double-bounce -> red
        3: (128,  0, 64,255),   # 3 branch/crown -> maroon
        4: (255,204,  0,255),   # 4 dipole -> yellow
        5: (0, 255, 0,255),     # 5 vegetation -> lime green
        6: (0,100,  0,255),     # 6 anisotropic needles -> dark green
        7: (0,204,204,255),     # 7 Bragg surface -> cyan
        8: (0, 51,153,255),     # 8 surface roughness -> dark blue
        9: (153, 51,153,255),   # 9 no feasible region -> violet
    }

CLASS_DESC = {
    0: "0: nodata / unclassified",
    1: "1: Dihedral scatterer (orange) [Z1]",
    2: "2: Forestry / double-bounce (red) [Z2]",
    3: "3: Branch / crown structure (maroon) [Z3]",
    4: "4: Dipole (yellow) [Z4]",
    5: "5: Vegetation (lime green) [Z5]",
    6: "6: Cloud of anisotropic needles (dark green) [Z6]",
    7: "7: Bragg surface (cyan) [Z7]",
    8: "8: Surface roughness / propagation (dark blue) [Z8]",
    9: "9: No feasible region / ambiguous (violet) [Z9]",
}

# ------------------ Write paletted GeoTIFF ------------------
def write_paletted_geotiff(
    path, arr_byte, geo=None, proj=None, nodata=0, class_desc=None
):
    """
    Write single-band paletted GeoTIFF with explicit palette entries 0..9.
    All remaining palette indices (10..255) are assigned transparent (0 alpha).
    Embeds class descriptions in dataset metadata as CLASS_0..CLASS_9.
    """
    if class_desc is None:
        class_desc = CLASS_DESC

    driver = gdal.GetDriverByName('GTiff')
    h, w = arr_byte.shape
    ds = driver.Create(
        path, w, h, 1, gdal.GDT_Byte, options=["COMPRESS=LZW", "TILED=YES"]
    )
    if geo:
        ds.SetGeoTransform(geo)
    if proj:
        ds.SetProjection(proj)

    band = ds.GetRasterBand(1)
    band.WriteArray(arr_byte.astype(np.uint8))
    band.SetNoDataValue(int(nodata))

    # Build color table: set 0..9 to our values; 10..255 -> transparent (0,0,0,0)
    ct = gdal.ColorTable()
    cmap = get_esa_colortable()
    # set defined entries 0..9
    for idx in sorted(cmap.keys()):
        r, g, b, a = cmap[idx]
        ct.SetColorEntry(int(idx), (int(r), int(g), int(b), int(a)))
    # set remaining entries to transparent (alpha=0) to avoid opaque black for unused indices
    for i in range(10, 256):
        ct.SetColorEntry(i, (0, 0, 0, 0))

    band.SetRasterColorTable(ct)
    try:
        band.SetRasterColorInterpretation(gdal.GCI_PaletteIndex)
    except Exception:
        pass

    # metadata describing each class
    meta = {}
    for k in sorted(class_desc.keys()):
        meta[f"CLASS_{int(k)}"] = class_desc[k]
    ds.SetMetadata(meta)

    try:
        band.SetDescription("H/alpha 9-zone class map (0=nodata, 1..9 classes)")
    except Exception:
        pass

    ds.FlushCache()
    ds = None

# ------------------ main ------------------
def main():
    """Main entry point for H/Alpha 9-zone classification."""
    p = argparse.ArgumentParser(
        description="H/Alpha 9-Zone Classifier",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    g_io = p.add_argument_group('IO')
    g_io.add_argument(
        '--h_file', required=True, help="Input Entropy (H) file (GeoTIFF or raw)"
    )
    g_io.add_argument(
        '--alpha_file',
        required=True,
        help="Input Alpha (alpha) file (GeoTIFF or raw)",
    )
    g_io.add_argument(
        '--out', default='ha9.tif', help="Output paletted GeoTIFF class map"
    )

    g_raw = p.add_argument_group('Raw Binary IO')
    g_raw.add_argument(
        '--width', type=int, help="Required if inputs are raw binary"
    )
    g_raw.add_argument(
        '--height', type=int, help="Required if inputs are raw binary"
    )
    g_raw.add_argument(
        '--dtype', default='float32', help="Raw binary data type (default: float32)"
    )
    g_raw.add_argument(
        '--byteorder', choices=('little', 'big'), default='little'
    )

    g_thresh = p.add_argument_group(
        'Classification Thresholds',
        description="""
    H/Alpha thresholds. Defaults follow 'best practice' H-dependent curves.
    Example: ... --h_thresh 0.5 0.9 --alpha_thresh_low_h 42.5 52.5 ...
    """,
    )
    g_thresh.add_argument(
        '--h_thresh',
        nargs=2,
        type=float,
        default=(0.5, 0.9),
        metavar=('H_LOW', 'H_HIGH'),
        help="Entropy thresholds (default: 0.5 0.9)",
    )
    g_thresh.add_argument(
        '--alpha_thresh_low_h',
        nargs=2,
        type=float,
        default=(42.5, 52.5),
        metavar=('A_L1', 'A_H1'),
        help="Alpha thresholds for Low H bin (default: 42.5 52.5)",
    )
    g_thresh.add_argument(
        '--alpha_thresh_med_h',
        nargs=2,
        type=float,
        default=(40.0, 50.0),
        metavar=('A_L2', 'A_H2'),
        help="Alpha thresholds for Med H bin (default: 40.0 50.0)",
    )
    g_thresh.add_argument(
        '--alpha_thresh_high_h',
        nargs=2,
        type=float,
        default=(45.0, 55.0),
        metavar=('A_L3', 'A_H3'),
        help="Alpha thresholds for High H bin (default: 45.0 55.0)",
    )
    args = p.parse_args()

    # Load inputs
    print("Loading H file:", args.h_file)
    H, geo, proj, nodata = ensure_loaded(
        args.h_file, args.width, args.height, args.dtype, args.byteorder
    )
    print("Loading Alpha file:", args.alpha_file)
    Alpha, geo2, proj2, nodata2 = ensure_loaded(
        args.alpha_file, args.width, args.height, args.dtype, args.byteorder
    )

    if geo is None and geo2 is not None:
        geo = geo2
    if proj is None and proj2 is not None:
        proj = proj2

    alpha_deg = convert_alpha_to_degrees(Alpha)

    # nodata mask
    if nodata is not None:
        mask1 = H == nodata
    else:
        mask1 = np.isnan(H)

    if nodata2 is not None:
        mask2 = Alpha == nodata2
    else:
        mask2 = np.isnan(Alpha)

    nodata_mask = mask1 | mask2

    # Clamp H values (can be slightly > 1.0 due to noise)
    H[H > 1.0] = 1.0
    H[H < 0.0] = 0.0

    print("Classifying zones...")
    cls = classify_9zones(
        H,
        alpha_deg,
        h_thresh=tuple(args.h_thresh),
        alpha_thresh_low_h=tuple(args.alpha_thresh_low_h),
        alpha_thresh_med_h=tuple(args.alpha_thresh_med_h),
        alpha_thresh_high_h=tuple(args.alpha_thresh_high_h),
        nodata_mask=nodata_mask,
    )

    # Write single paletted GeoTIFF
    print(f"Wrote paletted class GeoTIFF: {args.out}")
    write_paletted_geotiff(
        args.out, cls.astype(np.uint8), geo=geo, proj=proj, nodata=0
    )

    print("Thresholds used:")
    print(f"  H: {args.h_thresh}")
    print(f"  Alpha (Low H): {args.alpha_thresh_low_h}")
    print(f"  Alpha (Med H): {args.alpha_thresh_med_h}")
    print(f"  Alpha (High H): {args.alpha_thresh_high_h}")

    # Create quicklook RGB using the same palette
    # Build a 10-element lookup table (for classes 0-9)
    # 10 rows (one for each class 0-9), 3 columns (R,G,B)
    cmap = get_esa_colortable()
    h, w = cls.shape
    color_lookup = np.zeros((10, 3), dtype=np.uint8)
    for idx, (r, g, b, a) in cmap.items():
        if 0 <= idx <= 9:  # ensure valid class index
            color_lookup[idx] = [r, g, b]
    # Apply lookup table in one vectorized operation
    # This maps each class index in to its corresponding [R,G,B] triplet
    rgb = color_lookup[cls]
    quicklook = os.path.splitext(args.out)[0] + '_ql.png'
    plt.imsave(quicklook, rgb)
    print(f"Wrote quicklook PNG: {quicklook}")


if __name__ == '__main__':
    main()
