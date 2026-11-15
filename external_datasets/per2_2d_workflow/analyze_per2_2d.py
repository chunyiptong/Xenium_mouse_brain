#!/usr/bin/env python3
"""2D Per2 spot detection, raw-image photometry, and mask-based assignment.

See README.md. No channel selection, Z projection, or cell segmentation is
performed here. Inputs must already share the same unscaled 2D pixel grid.
"""
import argparse
import csv
import json
from importlib.metadata import version, PackageNotFoundError
from pathlib import Path

import numpy as np
import tifffile
from PIL import Image
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from roifile import roiread
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path as PolygonPath


def read_image(path):
    if Path(path).suffix.lower() in ('.tif', '.tiff'):
        a = tifffile.imread(path)
    else:
        a = np.asarray(Image.open(path))
    if a.ndim != 2:
        raise ValueError(f'{path}: expected a single grayscale 2D plane, got {a.shape}')
    return a


def write_table(path, rows, fields):
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def roi_membership(xy, rois):
    names = [[] for _ in range(len(xy))]
    for name, vertices in rois:
        # Tiny positive/negative buffer includes the boundary for either winding.
        poly = PolygonPath(np.vstack((vertices, vertices[0])), closed=True)
        hit = poly.contains_points(xy, radius=1e-7) | poly.contains_points(xy, radius=-1e-7)
        for i in np.flatnonzero(hit):
            names[i].append(name)
    return [';'.join(n) for n in names]


def analyze(raw, labels, points, scores, rois, args):
    sx, sy = args.pixel_size_x, args.pixel_size_y
    if args.expand_um:
        distance, index = ndi.distance_transform_edt(labels == 0, sampling=(sy, sx), return_indices=True)
        assigned = labels.copy()
        expand = (labels == 0) & (distance <= args.expand_um)
        assigned[expand] = labels[tuple(index[:, expand])]
    else:
        assigned = labels.copy()
    # A boundary flag is advisory; ambiguous assignments are retained for review.
    boundary = ndi.maximum_filter(assigned, size=3) != ndi.minimum_filter(assigned, size=3)
    xy = points[:, ::-1]
    roi_ids = roi_membership(xy, rois)
    tree = cKDTree(points) if len(points) else None
    rows = []
    h, w = raw.shape
    mode = ('cell_mask' if args.mask_kind == 'cells' else
            'estimated_cell_from_nucleus' if args.expand_um else 'nucleus_only')
    for i, (y, x) in enumerate(points):
        iy, ix = np.floor([y + .5, x + .5]).astype(int)
        iy, ix = np.clip(iy, 0, h-1), np.clip(ix, 0, w-1)
        radius = args.bg_outer
        y0, y1 = max(0, int(np.floor(y-radius))), min(h, int(np.ceil(y+radius))+1)
        x0, x1 = max(0, int(np.floor(x-radius))), min(w, int(np.ceil(x+radius))+1)
        yy, xx = np.mgrid[y0:y1, x0:x1]
        d2 = (yy-y)**2 + (xx-x)**2
        aperture = d2 <= args.spot_radius**2
        annulus = (d2 >= args.bg_inner**2) & (d2 <= radius**2)
        neighbours = [j for j in tree.query_ball_point([y, x], radius + args.spot_radius) if j != i]
        for j in neighbours:
            ny, nx = points[j]
            annulus &= (yy-ny)**2 + (xx-nx)**2 > args.spot_radius**2
        tile = raw[y0:y1, x0:x1].astype(float)
        n_bg = int(annulus.sum())
        background = float(np.median(tile[annulus])) if n_bg >= 10 else np.nan
        n_signal = int(aperture.sum())
        integrated = float(tile[aperture].sum() - background*n_signal) if n_signal else np.nan
        overlap = any(np.linalg.norm(points[j]-points[i]) < 2*args.spot_radius for j in neighbours)
        label = int(assigned[iy, ix])
        rows.append(dict(
            image_id=args.image_id, spot_id=f'p{i+1:06d}', roi_id=roi_ids[i],
            roi_overlap=';' in roi_ids[i], x_px=float(x), y_px=float(y),
            x_um=float(x*sx), y_um=float(y*sy), score=float(scores[i]),
            peak_raw=float(tile[aperture].max()) if n_signal else np.nan,
            background_local=background, integrated_bgsub=integrated,
            aperture_n_pixels=n_signal, background_n_pixels=n_bg,
            cell_id=f'{args.image_id}_c{label:05d}' if label else '',
            mask_label=label, assignment=mode if label else 'unassigned',
            near_mask_boundary=bool(boundary[iy, ix]),
            aperture_overlaps_other_spot=overlap,
            edge_window=bool(x-radius < 0 or y-radius < 0 or x+radius > w-1 or y+radius > h-1),
        ))
    ids = np.unique(assigned)
    ids = ids[ids > 0]
    centers = ndi.center_of_mass(np.ones_like(assigned, dtype=float), assigned, ids)
    cell_roi = roi_membership(np.array(centers)[:, ::-1], rois)
    counts = np.bincount(assigned.ravel())
    total_counts, scn_counts = {}, {}
    for row in rows:
        k = row['mask_label']
        total_counts[k] = total_counts.get(k, 0) + 1
        if row['roi_id']:
            scn_counts[k] = scn_counts.get(k, 0) + 1
    cells = []
    for k, (y, x), roi in zip(ids, centers, cell_roi):
        cells.append(dict(image_id=args.image_id, cell_id=f'{args.image_id}_c{k:05d}',
                          mask_label=int(k), assignment_method=mode,
                          roi_at_centroid=roi, centroid_x_px=x, centroid_y_px=y,
                          area_um2=float(counts[k]*sx*sy),
                          n_spots_all=total_counts.get(k, 0), n_spots_in_scn=scn_counts.get(k, 0)))
    return rows, cells, assigned


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--per2', required=True)
    p.add_argument('--mask', required=True, help='Cellpose integer label PNG/TIFF; background=0')
    p.add_argument('--mask-kind', required=True, choices=['cells', 'nuclei'])
    p.add_argument('--rois', required=True, help='Fiji Polygon/Freehand .roi or ROI .zip')
    p.add_argument('--image-id', required=True)
    p.add_argument('--image-mode', required=True, choices=['single_plane', 'max_projection'])
    p.add_argument('--source-description', required=True, help='ND2 filename, channel, Z plane/range, timepoint')
    p.add_argument('--pixel-size-x', required=True, type=float, help='Micrometers per pixel')
    p.add_argument('--pixel-size-y', required=True, type=float)
    p.add_argument('--expand-um', default=0., type=float, help='Optional nucleus expansion; default 0')
    p.add_argument('--spot-radius', default=2., type=float, help='Aperture radius in pixels; validate on images')
    p.add_argument('--bg-inner', default=3., type=float)
    p.add_argument('--bg-outer', default=5., type=float)
    p.add_argument('--prob-thresh', default=None, type=float, help='None uses pretrained model default')
    p.add_argument('--points-tsv', help='Optional precomputed table with y_px,x_px,score; bypass inference')
    p.add_argument('--out', required=True)
    args = p.parse_args()
    if min(args.pixel_size_x, args.pixel_size_y) <= 0 or args.expand_um < 0:
        p.error('Pixel sizes must be positive and expansion nonnegative')
    if not (1 <= args.spot_radius < args.bg_inner < args.bg_outer):
        p.error('Require 1 <= spot-radius < bg-inner < bg-outer')
    if args.mask_kind == 'cells' and args.expand_um:
        p.error('Expansion is only supported for nuclei masks')
    if args.prob_thresh is not None and not 0 <= args.prob_thresh <= 1:
        p.error('prob-thresh must be between 0 and 1')
    out = Path(args.out)
    if out.exists() and any(out.iterdir()):
        p.error('Output directory is nonempty; use a new directory to preserve previous results')
    raw, labels = read_image(args.per2), read_image(args.mask)
    if raw.shape != labels.shape or not np.isfinite(raw).all():
        p.error('Image/mask dimensions must match; image must contain finite values')
    if not np.issubdtype(labels.dtype, np.integer) or labels.min() < 0 or labels.max() == 0:
        p.error('Mask must have positive integer object IDs and background 0')
    if labels.max() > 1000000:
        p.error('Unexpectedly large label IDs; export a compact integer label image')
    if len(np.unique(labels[labels > 0])) == 1:
        print('WARNING: mask contains only one object. Verify it is a label image, not a binary mask.')
    records = roiread(args.rois)
    if not isinstance(records, list):
        records = [records]
    rois = []
    for i, roi in enumerate(records):
        if roi.roitype.name not in ('POLYGON', 'FREEHAND', 'TRACED') or roi.shape_roi_size:
            raise ValueError('Use simple Polygon/Freehand ROIs; composite ROIs are not supported')
        v = np.asarray(roi.coordinates(), dtype=float)
        if v.ndim != 2 or v.shape[1] != 2 or len(v) < 3:
            raise ValueError('Invalid polygon ROI')
        rois.append((roi.name or f'SCN_{i+1}', v))
    if not rois or len({n for n, _ in rois}) != len(rois):
        p.error('At least one ROI is required, with unique names')
    if args.points_tsv:
        with open(args.points_tsv) as f:
            supplied = list(csv.DictReader(f, delimiter='\t'))
        points = np.array([[float(r['y_px']), float(r['x_px'])] for r in supplied]).reshape(-1, 2)
        scores = np.array([float(r['score']) for r in supplied])
    else:
        from spotiflow.model import Spotiflow
        model = Spotiflow.from_pretrained('general')
        points, details = model.predict(raw, prob_thresh=args.prob_thresh)
        points = np.asarray(points).reshape(-1, 2)
        scores = np.asarray(details.prob).reshape(-1)
    if len(scores) != len(points) or not np.isfinite(points).all():
        raise ValueError('Invalid detection output')
    if len(points) and ((points < 0).any() or (points[:, 0] >= raw.shape[0]).any() or (points[:, 1] >= raw.shape[1]).any()):
        raise ValueError('Spot coordinates outside image')
    rows, cells, assigned = analyze(raw, labels, points, scores, rois, args)
    out.mkdir(parents=True, exist_ok=True)
    spot_fields = ('image_id spot_id roi_id roi_overlap x_px y_px x_um y_um score peak_raw background_local '
                   'integrated_bgsub aperture_n_pixels background_n_pixels cell_id mask_label assignment '
                   'near_mask_boundary aperture_overlaps_other_spot edge_window').split()
    write_table(out/'per2_spots.txt', rows, spot_fields)
    write_table(out/'cells.txt', cells, list(cells[0]))
    vertices = [dict(roi_id=name, vertex_order=i, x_px=float(x), y_px=float(y))
                for name, vv in rois for i, (x, y) in enumerate(vv)]
    write_table(out/'contours.txt', vertices, ['roi_id', 'vertex_order', 'x_px', 'y_px'])
    tifffile.imwrite(out/'assignment_labels.tif', assigned.astype(np.uint32))
    fig, ax = plt.subplots(figsize=(12, 8))
    low, high = np.percentile(raw, [1, 99.8])
    ax.imshow(raw, cmap='gray', vmin=low, vmax=high, origin='upper')
    edge = ndi.maximum_filter(assigned, 3) != ndi.minimum_filter(assigned, 3)
    overlay = np.zeros((*raw.shape, 4))
    overlay[edge] = (0, 1, 1, .6)
    ax.imshow(overlay, origin='upper')
    for name, v in rois:
        vv = np.vstack((v, v[0]))
        ax.plot(vv[:, 0], vv[:, 1], color='yellow', linewidth=1)
        ax.text(*v[0], name, color='yellow', fontsize=8)
    if len(points):
        colors = ['lime' if r['cell_id'] else 'magenta' for r in rows]
        ax.scatter(points[:, 1], points[:, 0], s=12, facecolors='none', edgecolors=colors, linewidths=.5)
    ax.set_title('2D Per2: assigned=green, unassigned=magenta, labels=cyan, SCN=yellow')
    ax.set_xlabel('x (0-based pixels)'); ax.set_ylabel('y (0-based pixels)')
    fig.tight_layout(); fig.savefig(out/'qc_overlay.png', dpi=200); plt.close(fig)
    packages = {}
    for name in ['spotiflow', 'numpy', 'scipy', 'roifile', 'tifffile', 'matplotlib']:
        try:
            packages[name] = version(name)
        except PackageNotFoundError:
            packages[name] = None
    metadata = dict(parameters=vars(args), shape_yx=list(raw.shape), versions=packages,
                    model='general' if not args.points_tsv else 'external_coordinates',
                    coordinate_convention='0-based; x rightward, y downward; um relative to image origin',
                    intensity_definition='Raw 2D input; aperture sum minus median annulus background times aperture pixel count. No deblending.',
                    notes='Default aperture sizes are initial values requiring validation. Projection measurements are 2D only. Expanded nuclear territories are estimated. Inspect boundary/overlap flags.')
    (out/'run_info.json').write_text(json.dumps(metadata, indent=2, ensure_ascii=False))
    print(f'Wrote {len(rows)} spots and {len(cells)} objects to {out.resolve()}')


if __name__ == '__main__':
    main()
