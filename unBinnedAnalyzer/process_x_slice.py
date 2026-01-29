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
parser.add_argument("--x-min", type=float, default=-150, help="Min X")
parser.add_argument("--x-max", type=float, default=150, help="Max X")
parser.add_argument("--x-nbins", type=int, default=15, help="Number of X bins")
args = parser.parse_args()

# theta bin edges
theta_edges_bin = np.array([0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 54, 65, 90])

# cut functions
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

# compute 3x3 subbins
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

# main processing function
def process_x_slice_slim(slice_file, save_to=None, chunk_size_mb=200, plane_index=0, x_min=-150, x_max=150, x_nbins=15):
    theta_edges = theta_edges_bin
    theta_centers = 0.25 * (theta_edges[1:] - theta_edges[:-1])
    n_theta_bins = len(theta_edges) - 1

    x_edges = np.linspace(x_min, x_max, x_nbins+1)
    n_x_bins = len(x_edges) - 1

    # initialize bins
    bins = {}
    for t_idx in range(n_theta_bins):
        for x_idx in range(n_x_bins):
            bins[(t_idx, x_idx)] = {
                "integral": [],
                "width": [],
                "theta": [],
                "x": []
            }

    # open ROOT file
    with uproot.open(slice_file) as f:
        tree = f["nominal"]
        branches = ["dirX", "pitch", "integral", "width", "x"]

        for chunk in tqdm(
            tree.iterate(branches, step_size=f"{chunk_size_mb} MB", library="np"),
            desc="Processing tree nominal",
            unit="chunk"
        ):
            thetaX = np.degrees(np.arctan(chunk["dirX"] * chunk["pitch"] / 0.3))
            x   = chunk["x"]
            integral = chunk["integral"]
            width    = chunk["width"]

            mask_integral = integral > cut_integral(thetaX)
            mask_width = width > cut_width(thetaX, plane_index)
            mask = mask_integral & mask_width

            thetaX = thetaX[mask]
            x   = x[mask]
            integral = integral[mask]
            width = width[mask]

            tb = np.digitize(thetaX, theta_edges) - 1
            xb = np.digitize(x, x_edges) - 1

            valid = (
                (tb >= 0) & (tb < n_theta_bins) &
                (xb >= 0) & (xb < n_x_bins)
            )

            for t, x, i, w, th, xx in zip(tb[valid], xb[valid], integral[valid], width[valid], thetaX[valid], x[valid]):
                b = bins[(t, x)]
                b["integral"].append(i)
                b["width"].append(w)
                b["theta"].append(th)
                b["x"].append(xx)

    # compute results
    results = []
    for (t_idx, x_idx), b in tqdm(bins.items(), desc="Computing ITM per (theta,X) bin"):
        theta_min, theta_max = theta_edges[t_idx], theta_edges[t_idx+1]
        x_min_bin, x_max_bin = x_edges[x_idx], x_edges[x_idx+1]

        n_entries = len(b["integral"])

        mean_i, err_i = iterative_truncated_mean(b["integral"])
        mean_w, err_w = iterative_truncated_mean(b["width"])

        sub_i = compute_3x3_subbins(
            b["integral"], np.array(b["theta"]), np.array(b["x"]),
            theta_min, theta_max, x_min_bin, x_max_bin
        )

        sub_w = compute_3x3_subbins(
            b["width"], np.array(b["theta"]), np.array(b["x"]),
            theta_min, theta_max, x_min_bin, x_max_bin
        )

        row = {
            "theta_bin": pd.Interval(theta_min, theta_max),
            "x_bin": pd.Interval(x_min_bin, x_max_bin),
            "n_entries": n_entries,
            "mean_integral": mean_i,
            "err_integral": err_i,
            "mean_width": mean_w,
            "err_width": err_w,
        }

        for k in range(1, 10):
            row[f"sub{k}_mean_integral"] = sub_i[k][0]
            row[f"sub{k}_err_integral"]  = sub_i[k][1]
            row[f"sub{k}_mean_width"]    = sub_w[k][0]
            row[f"sub{k}_err_width"]     = sub_w[k][1]

        results.append(row)

    df = pd.DataFrame(results)
    if save_to:
        df.to_pickle(save_to)
        print(f"Saved X slice with theta/X bins to {save_to}")

    return df

process_x_slice_slim(
    slice_file=args.input,
    save_to=args.output,
    chunk_size_mb=args.chunk_size,
    plane_index=args.plane,
    x_min=args.x_min,
    x_max=args.x_max,
    x_nbins=args.x_nbins
)
