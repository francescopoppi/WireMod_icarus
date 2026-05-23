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
args = parser.parse_args()

theta_edges_bin = np.array([0,6,12,18,24,30,38,46,56,90])

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

def save_theta_histograms_root(
    bins_pre,
    bins_post,
    theta_edges,
    output_root="theta_histograms.root"
):
    """
    Salva histogrammi ROOT compatibili via uproot.

    Per ogni bin di theta salva:
      - integral pre-cut
      - integral post-cut
      - width pre-cut
      - width post-cut
    """

    histograms = {}

    n_theta_bins = len(theta_edges) - 1

    for t_idx in range(n_theta_bins):

        theta_min = theta_edges[t_idx]
        theta_max = theta_edges[t_idx + 1]

        bpre = bins_pre[(t_idx, 0)]

        integral_pre = np.array(bpre["integral"], dtype=np.float64)
        width_pre    = np.array(bpre["width"], dtype=np.float64)



        if len(integral_pre) > 0:

            #xmin = np.percentile(integral_pre, 0.5)
            xmin = 0.
            #xmax = np.percentile(width_pre, 99.5)
            xmax = 2000.
            if xmin == xmax:
                xmax += 1.

            counts, edges = np.histogram(
                integral_pre,
                bins=400,
                range=(xmin, xmax)
            )

            histograms[
                f"hIntegral_pre_theta{t_idx}"
            ] = (counts, edges)

        # -------------------------
        # width pre
        # -------------------------

        if len(width_pre) > 0:

            #xmin = np.percentile(integral_pre, 0.5)
            xmin = 0.
            #xmax = np.percentile(width_pre, 99.5)
            xmax = 30.
            
            if xmin == xmax:
                xmax += 1.

            counts, edges = np.histogram(
                width_pre,
                bins=300,
                range=(xmin, xmax)
            )

            histograms[
                f"hWidth_pre_theta{t_idx}"
            ] = (counts, edges)

        # =========================================================
        # POST CUT
        # =========================================================

        bpost = bins_post[(t_idx, 0)]

        integral_post = np.array(bpost["integral"], dtype=np.float64)
        width_post    = np.array(bpost["width"], dtype=np.float64)

        # -------------------------
        # integral post
        # -------------------------

        if len(integral_post) > 0:

            #xmin = np.percentile(integral_pre, 0.5)
            xmin = 0.
            #xmax = np.percentile(width_pre, 99.5)
            xmax = 2000.

            if xmin == xmax:
                xmax += 1.

            counts, edges = np.histogram(
                integral_post,
                bins=400,
                range=(xmin, xmax)
            )

            histograms[
                f"hIntegral_post_theta{t_idx}"
            ] = (counts, edges)

        # -------------------------
        # width post
        # -------------------------

        if len(width_post) > 0:

            #xmin = np.percentile(integral_pre, 0.5)
            xmin = 0.
            #xmax = np.percentile(width_pre, 99.5)
            xmax = 30.
            
            if xmin == xmax:
                xmax += 1.

            counts, edges = np.histogram(
                width_post,
                bins=300,
                range=(xmin, xmax)
            )

            histograms[
                f"hWidth_post_theta{t_idx}"
            ] = (counts, edges)

    with uproot.recreate(output_root) as fout:
        for name, hist in histograms.items():
            fout[name] = hist

    print(f"Saved histograms to {output_root}")

def process_x_slice_slim(slice_file, save_to=None, chunk_size_mb=200, plane_index=0, x_min=-150, x_max=150, x_nbins=15):
    theta_edges = theta_edges_bin
    theta_centers = 0.25 * (theta_edges[1:] - theta_edges[:-1])
    n_theta_bins = len(theta_edges) - 1

    x_edges = np.array([x_min, x_max])
    n_x_bins = 1

    bins_pre = {}
    bins = {}
    for t_idx in range(n_theta_bins):
        for x_idx in range(n_x_bins):
            bins_pre[(t_idx, x_idx)] = {
                "integral": [],
                "width": []
            }
            bins[(t_idx, x_idx)] = {
                "integral": [],
                "width": []
            }

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

            for t_idx in range(n_theta_bins):
                for x_idx in range(n_x_bins):
                    mask_bin = valid & (tb == t_idx) & (xb == x_idx)
                    if mask_bin.any():
                        bins[(t_idx, x_idx)]["integral"].append(integral[mask_bin].astype(np.float32))
                        bins[(t_idx, x_idx)]["width"].append(width[mask_bin].astype(np.float32))

    results = []
    for (t_idx, x_idx), b in tqdm(bins.items(), desc="Computing ITM per (theta,X) bin"):
        theta_min, theta_max = theta_edges[t_idx], theta_edges[t_idx+1]
        x_min_bin, x_max_bin = x_edges[x_idx], x_edges[x_idx+1]

        arr_i = np.concatenate(b["integral"]) if b["integral"] else np.array([], dtype=np.float32)
        arr_w = np.concatenate(b["width"])    if b["width"]    else np.array([], dtype=np.float32)
        n_entries = len(arr_i)

        mean_i, err_i = iterative_truncated_mean(arr_i)
        mean_w, err_w = iterative_truncated_mean(arr_w)

        #sub_i = compute_3x3_subbins(
        #    b["integral"], np.array(b["theta"]), np.array(b["x"]),
        #    theta_min, theta_max, x_min_bin, x_max_bin
        #)
#
        #sub_w = compute_3x3_subbins(
        #    b["width"], np.array(b["theta"]), np.array(b["x"]),
        #    theta_min, theta_max, x_min_bin, x_max_bin
        #)

        row = {
            "theta_bin": pd.Interval(theta_min, theta_max),
            "x_bin": pd.Interval(x_min_bin, x_max_bin),
            "n_entries": n_entries,
            "mean_integral": mean_i,
            "err_integral": err_i,
            "mean_width": mean_w,
            "err_width": err_w,
        }

        #for k in range(1, 10):
        #    row[f"sub{k}_mean_integral"] = sub_i[k][0]
        #    row[f"sub{k}_err_integral"]  = sub_i[k][1]
        #    row[f"sub{k}_mean_width"]    = sub_w[k][0]
        #    row[f"sub{k}_err_width"]     = sub_w[k][1]

        results.append(row)
        save_theta_histograms_root(
            bins_pre=bins_pre,
            bins_post=bins,
            theta_edges=theta_edges,
            output_root=save_to.replace(".pkl", "_histo.root")
        )
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
    x_max=args.x_max
)
