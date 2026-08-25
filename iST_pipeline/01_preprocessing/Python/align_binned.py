"""
Purpose:  Rasterise fresh source (Xenium/MERSCOPE) and target (VisiumHD H&E)
          images, then fit an affine + LDDMM diffeomorphic registration
          (STalign) from a pre-selected set of source/target landmark points
          and apply the fitted transform to the source bin coordinates. The
          rasterised images are recomputed from raw inputs on every run —
          nothing is loaded from a pre-existing static image file.
Inputs:   --coords_csv : CSV with columns "cell_id","x","y" giving the
              platform's bin centroids — used both to rasterise the source
              image (uniform weighting, matching the validated STalign_
              reticulate.Rmd workflow) and as the points to transform.
          --target_image : PNG of the VisiumHD H&E image (values in [0, 1]),
              written by align_binned.R from the VisiumHD reference object's
              image slot.
          --source_landmarks / --target_landmarks : CSV files with columns
              "y","x" (one row per landmark, matching STalign's point
              convention) — converted from the original landmark-picker RDS
              output by align_binned.R.
Outputs:  --output_csv : CSV with columns "cell_id","x","y" — the input bin
              centroids mapped into VisiumHD (target) coordinate space.
Usage:    python align_binned.py \
              --target_image <path> \
              --source_landmarks <path> --target_landmarks <path> \
              --coords_csv <path> --output_csv <path>

NOTE ON POINT ORDER: STalign represents every point as (y, x) — row before
column — not (x, y). The landmark CSVs and the intermediate arrays built here
all follow that convention; only the final output is written as (x, y) to
match Seurat's expected column order.
"""

import argparse

import matplotlib
matplotlib.use("Agg")  # SLURM compute nodes are headless; STalign.rasterize() opens a figure internally
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from STalign import STalign


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--target_image", required=True,
                    help="PNG of the VisiumHD H&E image, values in [0, 1]")
    p.add_argument("--source_landmarks", required=True,
                    help="CSV of source landmarks, columns y,x")
    p.add_argument("--target_landmarks", required=True,
                    help="CSV of target landmarks, columns y,x")
    p.add_argument("--coords_csv", required=True,
                    help="CSV of bin coordinates to rasterise and transform, columns cell_id,x,y")
    p.add_argument("--output_csv", required=True,
                    help="Output path for transformed coordinates, columns cell_id,x,y")
    p.add_argument("--dx", type=float, default=30,
                    help="Source rasterisation grid resolution in microns [default: %(default)s]")
    p.add_argument("--device", default="auto",
                    help="'auto' (cuda if available, else cpu), or an explicit torch device")
    p.add_argument("--niter", type=int, default=500)
    p.add_argument("--sigmaP", type=float, default=2e-1)
    p.add_argument("--sigmaM", type=float, default=0.2)
    p.add_argument("--sigmaB", type=float, default=0.3)
    p.add_argument("--sigmaA", type=float, default=0.3)
    p.add_argument("--diffeo_start", type=int, default=100)
    p.add_argument("--epL", type=float, default=5e-11)
    p.add_argument("--epT", type=float, default=5e-4)
    p.add_argument("--epV", type=float, default=5e1)
    return p.parse_args()


def main():
    args = parse_args()

    device = args.device
    if device == "auto":
        device = "cuda:0" if torch.cuda.is_available() else "cpu"
    print(f"Using device: {device}")
    torch.set_default_device(device)

    # Landmarks: (y, x) order, one row per named landmark
    pointsI = pd.read_csv(args.source_landmarks)[["y", "x"]].to_numpy()
    pointsJ = pd.read_csv(args.target_landmarks)[["y", "x"]].to_numpy()
    print(f"Loaded {pointsI.shape[0]} source landmarks, {pointsJ.shape[0]} target landmarks")

    # Bin coordinates: rasterised fresh into the source image below, and also
    # the points transformed at the end of this script
    coords_df = pd.read_csv(args.coords_csv)

    # --- Rasterise source (platform) image from raw bin centroids -------------
    # Uniform weighting (not nCount-weighted) — matches the ##fix chunk in
    # STalign_reticulate.Rmd, which is what the validated samples (e.g. wt709)
    # were registered with.
    xI = coords_df["x"].to_numpy()
    yI = coords_df["y"].to_numpy()
    gI = np.ones(len(coords_df))
    XI, YI, I, _ = STalign.rasterize(xI, yI, g=gI, dx=args.dx)
    I = np.vstack((I, I, I))  # make into 3xNxM
    I = STalign.normalize(I)

    # --- Build target (VisiumHD H&E) image from the raw image PNG -------------
    V = plt.imread(args.target_image)
    print(f"Target image shape: {V.shape}")
    Vnorm = STalign.normalize(V)
    J = Vnorm.transpose(2, 0, 1)
    YJ = np.arange(J.shape[1]) * 1.  # needs to be float, not int, for STalign.transform later
    XJ = np.arange(J.shape[2]) * 1.

    # Initial affine transform from landmark point pairs
    L, T = STalign.L_T_from_points(pointsI, pointsJ)
    A = STalign.to_A(torch.tensor(L), torch.tensor(T))

    # Full diffeomorphic (LDDMM) registration, guided by landmarks + images.
    # Parameters below match STalign_reticulate.Rmd exactly.
    print(f"Running LDDMM for {args.niter} iterations...")
    out = STalign.LDDMM(
        [YI, XI], I, [YJ, XJ], J,
        L=L, T=T,
        niter=args.niter,
        pointsI=pointsI, pointsJ=pointsJ,
        device=device,
        sigmaP=args.sigmaP, sigmaM=args.sigmaM, sigmaB=args.sigmaB, sigmaA=args.sigmaA,
        diffeo_start=args.diffeo_start,
        epL=args.epL, epT=args.epT, epV=args.epV,
    )

    A_fit = out["A"]
    v = out["v"]
    xv = out["xv"]

    # Apply the fitted transform to the requested bin coordinates. STalign
    # points are (y, x); the input CSV is (x, y), so stack in swapped order.
    points_yx = np.stack([coords_df["y"].to_numpy(), coords_df["x"].to_numpy()], axis=1)

    # transform_points_source_to_target requires float64 inputs on GPU;
    # otherwise it fails with "result type Double can't be cast to the
    # desired output type Long"
    xv = [x.double() for x in xv]
    v = v.double()
    points_yx = points_yx.double() if hasattr(points_yx, "double") else torch.tensor(points_yx, dtype=torch.float64)

    tpoints = STalign.transform_points_source_to_target(xv, v, A_fit, points_yx)
    tpoints = tpoints.detach().cpu().numpy()

    # Switch back from (y, x) to (x, y) for the output, matching Seurat's
    # expected centroid coordinate column order
    aligned_y = tpoints[:, 0]
    aligned_x = tpoints[:, 1]

    out_df = pd.DataFrame({
        "cell_id": coords_df["cell_id"],
        "x": aligned_x,
        "y": aligned_y,
    })
    out_df.to_csv(args.output_csv, index=False)
    print(f"Saved {len(out_df)} transformed coordinates to {args.output_csv}")


if __name__ == "__main__":
    main()
