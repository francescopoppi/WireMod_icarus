import argparse
import os
import numpy as np
import pandas as pd
import uproot
from tqdm import tqdm
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

parser = argparse.ArgumentParser()

parser.add_argument("--input",  required=True, help="Slim ROOT file")
parser.add_argument("--output", required=True, help="Output PKL")

parser.add_argument("--chunk-size", type=int, default=200)
parser.add_argument("--plane", type=int, required=True, choices=[0, 1, 2])

parser.add_argument("--y-min",    type=float, default=-180)
parser.add_argument("--y-max",    type=float, default=130)
parser.add_argument("--n-y-bins", type=int,   default=15)

parser.add_argument("--z-min",    type=float, default=-900)
parser.add_argument("--z-max",    type=float, default=900)
parser.add_argument("--n-z-bins", type=int,   default=20)

parser.add_argument(
    "--bootstrap-iters",
    type=int,
    default=200,
    help="Bootstrap samples for ITM error estimate (0 = classical SEM)"
)
parser.add_argument(
    "--bootstrap-max-n",
    type=int,
    default=10000,
    help="Bins with more than N entries use SEM instead of bootstrap"
)
parser.add_argument(
    "--bootstrap-seed",
    type=int,
    default=42,
    help="RNG seed for bootstrap"
)
parser.add_argument(
    "--save-root",
    action="store_true",
    help="Save ROOT histograms after PKL creation"
)

args = parser.parse_args()

RNG = np.random.default_rng(args.bootstrap_seed)


# =========================================================
# CUTS
# =========================================================

def cut_integral(x):
    return 1.35696 * np.exp(0.0786976 * x) - 24.1874


cut_width_params = [
    (0.108951, 0.0574775, 0.850119),
    (0.0444325, 0.0709458, 0.455733),
    (0.0833033, 0.0544162, 0.732181)
]


def cut_width(x, plane_index):
    p = cut_width_params[plane_index]
    return p[0] * np.exp(p[1] * x) + p[2]


# =========================================================
# ITM
# =========================================================

def _itm_core(values, sig_down=-2, sig_up=1.75, tol=1e-4, max_iter=100):
    if len(values) == 0:
        return 0., np.array([], dtype=np.float64)
    if len(values) == 1:
        return values[0], values

    working = values
    result  = np.mean(working)

    for _ in range(max_iter):
        mean   = np.mean(working)
        rms    = np.std(working)
        if rms == 0:
            return mean, working
        median = np.median(working)
        mask   = (
            (working >= median + sig_down * rms) &
            (working <= median + sig_up   * rms)
        )
        new_values = working[mask]
        if len(new_values) == 0:
            break
        new_mean = np.mean(new_values)
        if np.abs(new_mean - mean) < tol:
            return new_mean, new_values
        working = new_values
        result  = new_mean

    return result, working


def iterative_truncated_mean(
    values,
    sig_down=-2,
    sig_up=1.75,
    tol=1e-4,
    max_iter=100,
    n_bootstrap=0,
    bootstrap_max_n=10000,
    rng=None
):
    values = np.asarray(values, dtype=np.float64)
    if len(values) == 0:
        return 0., 0.
    if len(values) == 1:
        return values[0], 0.

    final_mean, final_values = _itm_core(values, sig_down, sig_up, tol, max_iter)

    use_bootstrap = (
        n_bootstrap > 1
        and len(values) > 1
        and len(values) <= bootstrap_max_n
    )

    if use_bootstrap:
        if rng is None:
            rng = np.random.default_rng()
        n = len(values)
        boot_means = np.empty(n_bootstrap, dtype=np.float64)
        for i in range(n_bootstrap):
            sample        = rng.choice(values, size=n, replace=True)
            boot_m, _     = _itm_core(sample, sig_down, sig_up, tol, max_iter)
            boot_means[i] = boot_m
        err = float(np.std(boot_means, ddof=1))
    else:
        err = (
            float(np.std(final_values) / np.sqrt(len(final_values)))
            if len(final_values) > 1 else 0.
        )

    return float(final_mean), err


# =========================================================
# ROOT HISTOGRAMS
# =========================================================

def save_yz_histograms_root(bins_pre, bins_post, y_edges, z_edges, output_root):

    histograms = {}
    n_y_bins = len(y_edges) - 1
    n_z_bins = len(z_edges) - 1

    for y_idx in range(n_y_bins):
        for z_idx in range(n_z_bins):
            for tag, bdict in (("pre", bins_pre), ("post", bins_post)):
                b = bdict[(y_idx, z_idx)]

                arr_i = (
                    np.concatenate(b["integral"]).astype(np.float64)
                    if b["integral"] else np.array([], dtype=np.float64)
                )
                arr_w = (
                    np.concatenate(b["width"]).astype(np.float64)
                    if b["width"] else np.array([], dtype=np.float64)
                )

                if len(arr_i) > 0:
                    counts, edges = np.histogram(arr_i, bins=400, range=(0., 2000.))
                    histograms[f"hIntegral_{tag}_y{y_idx}_z{z_idx}"] = (counts, edges)

                if len(arr_w) > 0:
                    counts, edges = np.histogram(arr_w, bins=300, range=(0., 30.))
                    histograms[f"hWidth_{tag}_y{y_idx}_z{z_idx}"] = (counts, edges)

    tmp_root = output_root + ".tmp"
    with uproot.recreate(tmp_root) as fout:
        for name, hist in histograms.items():
            fout[name] = hist
    os.replace(tmp_root, output_root)
    print(f"Saved ROOT histograms -> {output_root}")


# =========================================================
# MAIN
# =========================================================

def process_yz_slim(
    slice_file,
    save_to=None,
    chunk_size_mb=200,
    plane_index=0,
    y_min=-180,
    y_max=130,
    n_y_bins=15,
    z_min=-900,
    z_max=900,
    n_z_bins=20
):
    y_edges  = np.linspace(y_min, y_max, n_y_bins + 1)
    z_edges  = np.linspace(z_min, z_max, n_z_bins + 1)

    bins_pre  = {}
    bins_post = {}
    for y_idx in range(n_y_bins):
        for z_idx in range(n_z_bins):
            bins_pre [(y_idx, z_idx)] = {"integral": [], "width": []}
            bins_post[(y_idx, z_idx)] = {"integral": [], "width": []}

    # =====================================================
    # LOOP ROOT
    # =====================================================

    with uproot.open(slice_file) as f:
        tree     = f["nominal"]
        branches = ["dirX", "pitch", "integral", "width", "y", "z"]

        for chunk in tqdm(
            tree.iterate(branches, step_size=f"{chunk_size_mb} MB", library="np"),
            desc="Processing tree nominal",
            unit="chunk"
        ):
            thetaX   = np.degrees(np.arctan(chunk["dirX"] * chunk["pitch"] / 0.3))
            y        = chunk["y"]
            z        = chunk["z"]
            integral = chunk["integral"]
            width    = chunk["width"]

            # =============================================
            # PRE-CUT STORAGE
            # =============================================

            yb_pre = np.digitize(y, y_edges) - 1
            zb_pre = np.digitize(z, z_edges) - 1

            valid_pre = (
                (yb_pre >= 0) & (yb_pre < n_y_bins) &
                (zb_pre >= 0) & (zb_pre < n_z_bins)
            )

            for y_idx in range(n_y_bins):
                for z_idx in range(n_z_bins):
                    mask_bin_pre = (
                        valid_pre &
                        (yb_pre == y_idx) &
                        (zb_pre == z_idx)
                    )
                    if mask_bin_pre.any():
                        bins_pre[(y_idx, z_idx)]["integral"].append(
                            integral[mask_bin_pre].astype(np.float32)
                        )
                        bins_pre[(y_idx, z_idx)]["width"].append(
                            width[mask_bin_pre].astype(np.float32)
                        )

            # =============================================
            # APPLY CUTS
            # =============================================

            mask = (
                (integral > cut_integral(thetaX)) &
                (width    > cut_width(thetaX, plane_index))
            )

            y        = y[mask]
            z        = z[mask]
            integral = integral[mask]
            width    = width[mask]

            # =============================================
            # POST-CUT STORAGE
            # =============================================

            yb = np.digitize(y, y_edges) - 1
            zb = np.digitize(z, z_edges) - 1

            valid = (
                (yb >= 0) & (yb < n_y_bins) &
                (zb >= 0) & (zb < n_z_bins)
            )

            for y_idx in range(n_y_bins):
                for z_idx in range(n_z_bins):
                    mask_bin = (
                        valid &
                        (yb == y_idx) &
                        (zb == z_idx)
                    )
                    if mask_bin.any():
                        bins_post[(y_idx, z_idx)]["integral"].append(
                            integral[mask_bin].astype(np.float32)
                        )
                        bins_post[(y_idx, z_idx)]["width"].append(
                            width[mask_bin].astype(np.float32)
                        )

    # =====================================================
    # COMPUTE ITM
    # =====================================================

    results  = []
    n_boot   = args.bootstrap_iters
    max_boot = args.bootstrap_max_n

    for (y_idx, z_idx), b in tqdm(
        bins_post.items(),
        desc="Computing ITM per (Y, Z) bin"
    ):
        y_min_bin = y_edges[y_idx]
        y_max_bin = y_edges[y_idx + 1]
        z_min_bin = z_edges[z_idx]
        z_max_bin = z_edges[z_idx + 1]

        arr_i = (
            np.concatenate(b["integral"])
            if b["integral"] else np.array([], dtype=np.float32)
        )
        arr_w = (
            np.concatenate(b["width"])
            if b["width"] else np.array([], dtype=np.float32)
        )

        n_entries = len(arr_i)

        mean_i, err_i = iterative_truncated_mean(
            arr_i, n_bootstrap=n_boot, bootstrap_max_n=max_boot, rng=RNG
        )
        mean_w, err_w = iterative_truncated_mean(
            arr_w, n_bootstrap=n_boot, bootstrap_max_n=max_boot, rng=RNG
        )

        results.append({
            "y_min":         y_min_bin,
            "y_max":         y_max_bin,
            "z_min":         z_min_bin,
            "z_max":         z_max_bin,
            "n_entries":     int(n_entries),
            "mean_integral": float(mean_i),
            "err_integral":  float(err_i),
            "mean_width":    float(mean_w),
            "err_width":     float(err_w),
            "err_method": (
                "bootstrap"
                if n_boot > 1 and n_entries <= max_boot
                else "sem"
            ),
        })

    # =====================================================
    # SAVE PKL
    # =====================================================

    df = pd.DataFrame(results)

    if save_to:
        tmp_pkl = save_to + ".tmp"
        df.to_pickle(tmp_pkl)
        os.replace(tmp_pkl, save_to)
        print(f"Saved PKL -> {save_to}")

    # =====================================================
    # SAVE ROOT
    # =====================================================

    if args.save_root and save_to:
        root_output = save_to.replace(".pkl", "_histo.root")
        save_yz_histograms_root(
            bins_pre=bins_pre,
            bins_post=bins_post,
            y_edges=y_edges,
            z_edges=z_edges,
            output_root=root_output
        )

    return df


# =========================================================
# RUN
# =========================================================

process_yz_slim(
    slice_file=args.input,
    save_to=args.output,
    chunk_size_mb=args.chunk_size,
    plane_index=args.plane,
    y_min=args.y_min,
    y_max=args.y_max,
    n_y_bins=args.n_y_bins,
    z_min=args.z_min,
    z_max=args.z_max,
    n_z_bins=args.n_z_bins,
)
