import argparse
import numpy as np
import pandas as pd
import uproot
from tqdm import tqdm
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True, help="Slim ROOT file")
parser.add_argument("--output", required=True, help="Output PKL")
parser.add_argument("--chunk-size", type=int, default=200)
parser.add_argument("--plane", type=int, required=True, choices=[0,1,2])
parser.add_argument("--z-min", type=float, default=-150, help="Min Z")
parser.add_argument("--z-max", type=float, default=150, help="Max Z")
args = parser.parse_args()


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

def iterative_truncated_mean(values, sig_down=-2, sig_up=1.75, tol=1e-4, max_iter=100):
    values = np.array(values)
    if len(values) == 0:
        return 0., 0.
    result = np.mean(values)
    for i in range(max_iter):
        mean = np.mean(values)
        rms  = np.std(values)
        median = np.median(values)

        mask = (values >= median + sig_down*rms) & (values <= median + sig_up*rms)
        new_values = values[mask]

        if len(new_values) == 0:
            break

        new_mean = np.mean(new_values)

        if np.abs(new_mean - mean) < tol:
            return new_mean, np.std(new_values) / np.sqrt(len(new_values))

        values = new_values
        result = new_mean

    return result, np.std(values)/np.sqrt(len(values))

def compute_3x3_subbins(values, theta_vals, x_vals, theta_min, theta_max, x_min, x_max):
    results = {}
    theta_edges_sub = np.linspace(theta_min, theta_max, 4)
    x_edges_sub     = np.linspace(x_min, x_max, 4)
    sub_id = 1
    for j in range(3):        # theta direction
        for i in range(3):    # X direction
            t_lo, t_hi = theta_edges_sub[j], theta_edges_sub[j+1]
            x_lo, x_hi = x_edges_sub[i], x_edges_sub[i+1]

            mask = (
                (theta_vals >= t_lo) & (theta_vals < t_hi) &
                (x_vals >= x_lo) & (x_vals < x_hi)
            )

            sub_vals = np.array(values)[mask]
            if len(sub_vals) > 0:
                m, e = iterative_truncated_mean(sub_vals)
            else:
                m, e = np.nan, np.nan

            results[sub_id] = (m, e)
            sub_id += 1

    return results

def process_z_slice_slim(slice_file, save_to=None, chunk_size_mb=200, plane_index=0,
                         z_min=-150, z_max=150):
    y_edges = np.linspace(-180, 120, 15+1)   # 31 bin Y
    n_y_bins = len(y_edges) - 1

    z_edges = np.array([z_min, z_max])       # slice Z già fissata
    n_z_bins = len(z_edges) - 1

    bins = {}
    for y_idx in range(n_y_bins):
        for z_idx in range(n_z_bins):
            bins[(y_idx, z_idx)] = {
                "integral": [],
                "width": [],
                "y": [],
                "z": []
            }

    with uproot.open(slice_file) as f:
        tree = f["nominal"]
        branches = ["dirX", "pitch", "integral", "width", "y", "z"]

        for chunk in tqdm(
            tree.iterate(branches, step_size=f"{chunk_size_mb} MB", library="np"),
            desc="Processing tree nominal",
            unit="chunk"
        ):
            thetaX = np.degrees(np.arctan(chunk["dirX"] * chunk["pitch"] / 0.3))

            y        = chunk["y"]
            z        = chunk["z"]
            integral = chunk["integral"]
            width    = chunk["width"]

            mask_integral = integral > cut_integral(thetaX)
            mask_width    = width > cut_width(thetaX, plane_index)
            mask = mask_integral & mask_width

            y        = y[mask]
            z        = z[mask]
            integral = integral[mask]
            width    = width[mask]

            yb = np.digitize(y, y_edges) - 1
            zb = np.digitize(z, z_edges) - 1

            valid = (yb >= 0) & (yb < n_y_bins) & (zb >= 0) & (zb < n_z_bins)

            for y_idx, z_idx, i, w, yy, zz in zip(yb[valid], zb[valid],
                                                  integral[valid], width[valid],
                                                  y[valid], z[valid]):
                b = bins[(y_idx, z_idx)]
                b["integral"].append(i)
                b["width"].append(w)
                b["y"].append(yy)
                b["z"].append(zz)

    results = []
    for (y_idx, z_idx), b in tqdm(bins.items(), desc="Computing ITM per (Y,Z) bin"):
        y_min, y_max = y_edges[y_idx], y_edges[y_idx+1]
        z_min_bin, z_max_bin = z_edges[z_idx], z_edges[z_idx+1]

        n_entries = len(b["integral"])
        mean_i, err_i = iterative_truncated_mean(b["integral"])
        mean_w, err_w = iterative_truncated_mean(b["width"])

        row = {
            "y_bin": pd.Interval(y_min, y_max),
            "z_bin": pd.Interval(z_min_bin, z_max_bin),
            "n_entries": n_entries,
            "mean_integral": mean_i,
            "err_integral": err_i,
            "mean_width": mean_w,
            "err_width": err_w
        }

        results.append(row)

    df = pd.DataFrame(results)
    if save_to:
        df.to_pickle(save_to)
        print(f"Saved Z slice with Y/Z bins to {save_to}")

    return df

process_z_slice_slim(
    slice_file=args.input,
    save_to=args.output,
    chunk_size_mb=args.chunk_size,
    plane_index=args.plane,
    z_min=args.z_min,
    z_max=args.z_max
)
