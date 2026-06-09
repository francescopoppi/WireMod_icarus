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

parser.add_argument("--input", required=True, help="Slim ROOT file")
parser.add_argument("--output", required=True, help="Output PKL")

parser.add_argument("--chunk-size", type=int, default=200)
parser.add_argument("--plane", type=int, required=True, choices=[0,1,2])

parser.add_argument("--x-min", type=float, default=-150)
parser.add_argument("--x-max", type=float, default=150)

# Bootstrap per la stima dell'errore sull'ITM.
# --bootstrap-iters 0 disabilita il bootstrap e ricade sulla SEM classica.
# --bootstrap-max-n: per bin con più eventi di questa soglia si usa la SEM
#   (per grandi N il bootstrap converge alla SEM, quindi non serve calcolarlo).
parser.add_argument(
    "--bootstrap-iters",
    type=int,
    default=200,
    help="Numero di campioni bootstrap per la stima dell'errore ITM (0 = SEM classica)"
)
parser.add_argument(
    "--bootstrap-max-n",
    type=int,
    default=10000,
    help="Soglia: se il bin ha più di N eventi si usa la SEM invece del bootstrap"
)
parser.add_argument(
    "--bootstrap-seed",
    type=int,
    default=42,
    help="Seed per il generatore di numeri casuali del bootstrap"
)

# opzionale: crea anche il ROOT histo DOPO il pickle
parser.add_argument(
    "--save-root",
    action="store_true",
    help="Save ROOT histograms after PKL creation"
)

args = parser.parse_args()

theta_edges_bin = np.array([0,6,12,18,24,30,38,46,56,90])

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
    """
    Calcola solo il valore medio dell'ITM (senza errore).
    Usato internamente dal bootstrap.
    """
    if len(values) == 0:
        return 0., np.array([], dtype=np.float64)
    if len(values) == 1:
        return values[0], values

    working = values
    result = np.mean(working)

    for _ in range(max_iter):
        mean = np.mean(working)
        rms  = np.std(working)
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
    """
    Restituisce (mean, err) della media troncata iterativa.

    Stima dell'errore:
      - Se n_bootstrap > 0  e  len(values) <= bootstrap_max_n:
          bootstrap resampling sul dataset originale (pre-taglio).
          Cattura l'incertezza dovuta al confine di taglio adattivo,
          più accurato della SEM sul campione troncato per N moderati.
      - Altrimenti:
          SEM classica sul campione post-taglio: std/sqrt(N_trunc).
          Per grandi N il bootstrap converge a questo valore,
          quindi non aggiunge informazione e costa meno.
    """
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
            sample = rng.choice(values, size=n, replace=True)
            boot_m, _ = _itm_core(sample, sig_down, sig_up, tol, max_iter)
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

def save_theta_histograms_root(
    bins_pre,
    bins_post,
    theta_edges,
    output_root
):

    histograms = {}

    n_theta_bins = len(theta_edges) - 1

    for t_idx in range(n_theta_bins):

        # =====================================================
        # PRE
        # =====================================================

        bpre = bins_pre[(t_idx, 0)]

        integral_pre = (
            np.concatenate(bpre["integral"]).astype(np.float64)
            if bpre["integral"]
            else np.array([], dtype=np.float64)
        )

        width_pre = (
            np.concatenate(bpre["width"]).astype(np.float64)
            if bpre["width"]
            else np.array([], dtype=np.float64)
        )

        # integral pre

        if len(integral_pre) > 0:

            counts, edges = np.histogram(
                integral_pre,
                bins=400,
                range=(0., 2000.)
            )

            histograms[f"hIntegral_pre_theta{t_idx}"] = (
                counts,
                edges
            )

        # width pre

        if len(width_pre) > 0:

            counts, edges = np.histogram(
                width_pre,
                bins=300,
                range=(0., 30.)
            )

            histograms[f"hWidth_pre_theta{t_idx}"] = (
                counts,
                edges
            )

        # =====================================================
        # POST
        # =====================================================

        bpost = bins_post[(t_idx, 0)]

        integral_post = (
            np.concatenate(bpost["integral"]).astype(np.float64)
            if bpost["integral"]
            else np.array([], dtype=np.float64)
        )

        width_post = (
            np.concatenate(bpost["width"]).astype(np.float64)
            if bpost["width"]
            else np.array([], dtype=np.float64)
        )

        # integral post

        if len(integral_post) > 0:

            counts, edges = np.histogram(
                integral_post,
                bins=400,
                range=(0., 2000.)
            )

            histograms[f"hIntegral_post_theta{t_idx}"] = (
                counts,
                edges
            )

        # width post

        if len(width_post) > 0:

            counts, edges = np.histogram(
                width_post,
                bins=300,
                range=(0., 30.)
            )

            histograms[f"hWidth_post_theta{t_idx}"] = (
                counts,
                edges
            )

    tmp_root = output_root + ".tmp"

    with uproot.recreate(tmp_root) as fout:

        for name, hist in histograms.items():
            fout[name] = hist

    os.replace(tmp_root, output_root)

    print(f"Saved ROOT histograms -> {output_root}")


# =========================================================
# MAIN
# =========================================================

def process_x_slice_slim(
    slice_file,
    save_to=None,
    chunk_size_mb=200,
    plane_index=0,
    x_min=-150,
    x_max=150
):

    theta_edges = theta_edges_bin

    n_theta_bins = len(theta_edges) - 1

    x_edges = np.array([x_min, x_max])
    n_x_bins = 1

    bins_pre = {}
    bins_post = {}

    for t_idx in range(n_theta_bins):
        for x_idx in range(n_x_bins):

            bins_pre[(t_idx, x_idx)] = {
                "integral": [],
                "width": []
            }

            bins_post[(t_idx, x_idx)] = {
                "integral": [],
                "width": []
            }

    # =====================================================
    # LOOP ROOT
    # =====================================================

    with uproot.open(slice_file) as f:

        tree = f["nominal"]

        branches = [
            "dirX",
            "pitch",
            "integral",
            "width",
            "x"
        ]

        for chunk in tqdm(
            tree.iterate(
                branches,
                step_size=f"{chunk_size_mb} MB",
                library="np"
            ),
            desc="Processing tree nominal",
            unit="chunk"
        ):

            thetaX = np.degrees(
                np.arctan(
                    chunk["dirX"] * chunk["pitch"] / 0.3
                )
            )

            x = chunk["x"]

            integral = chunk["integral"]
            width = chunk["width"]

            # =============================================
            # PRE-CUT STORAGE
            # =============================================

            tb_pre = np.digitize(thetaX, theta_edges) - 1
            xb_pre = np.digitize(x, x_edges) - 1

            valid_pre = (
                (tb_pre >= 0) &
                (tb_pre < n_theta_bins) &
                (xb_pre >= 0) &
                (xb_pre < n_x_bins)
            )

            for t_idx in range(n_theta_bins):
                for x_idx in range(n_x_bins):

                    mask_bin_pre = (
                        valid_pre &
                        (tb_pre == t_idx) &
                        (xb_pre == x_idx)
                    )

                    if mask_bin_pre.any():

                        bins_pre[(t_idx, x_idx)]["integral"].append(
                            integral[mask_bin_pre].astype(np.float32)
                        )

                        bins_pre[(t_idx, x_idx)]["width"].append(
                            width[mask_bin_pre].astype(np.float32)
                        )

            # =============================================
            # APPLY CUTS
            # =============================================

            mask_integral = (
                integral > cut_integral(thetaX)
            )

            mask_width = (
                width > cut_width(thetaX, plane_index)
            )

            mask = mask_integral & mask_width

            thetaX = thetaX[mask]
            x = x[mask]
            integral = integral[mask]
            width = width[mask]

            # =============================================
            # POST-CUT STORAGE
            # =============================================

            tb = np.digitize(thetaX, theta_edges) - 1
            xb = np.digitize(x, x_edges) - 1

            valid = (
                (tb >= 0) &
                (tb < n_theta_bins) &
                (xb >= 0) &
                (xb < n_x_bins)
            )

            for t_idx in range(n_theta_bins):
                for x_idx in range(n_x_bins):

                    mask_bin = (
                        valid &
                        (tb == t_idx) &
                        (xb == x_idx)
                    )

                    if mask_bin.any():

                        bins_post[(t_idx, x_idx)]["integral"].append(
                            integral[mask_bin].astype(np.float32)
                        )

                        bins_post[(t_idx, x_idx)]["width"].append(
                            width[mask_bin].astype(np.float32)
                        )

    # =====================================================
    # COMPUTE ITM
    # =====================================================

    results = []

    n_boot   = args.bootstrap_iters
    max_boot = args.bootstrap_max_n

    for (t_idx, x_idx), b in tqdm(
        bins_post.items(),
        desc="Computing ITM per (theta,X) bin"
    ):

        theta_min = theta_edges[t_idx]
        theta_max = theta_edges[t_idx + 1]

        x_min_bin = x_edges[x_idx]
        x_max_bin = x_edges[x_idx + 1]

        arr_i = (
            np.concatenate(b["integral"])
            if b["integral"]
            else np.array([], dtype=np.float32)
        )

        arr_w = (
            np.concatenate(b["width"])
            if b["width"]
            else np.array([], dtype=np.float32)
        )

        n_entries = len(arr_i)

        mean_i, err_i = iterative_truncated_mean(
            arr_i,
            n_bootstrap=n_boot,
            bootstrap_max_n=max_boot,
            rng=RNG
        )
        mean_w, err_w = iterative_truncated_mean(
            arr_w,
            n_bootstrap=n_boot,
            bootstrap_max_n=max_boot,
            rng=RNG
        )

        row = {

            "theta_min": theta_min,
            "theta_max": theta_max,

            "x_min": x_min_bin,
            "x_max": x_max_bin,

            "n_entries": int(n_entries),

            "mean_integral": float(mean_i),
            "err_integral":  float(err_i),

            "mean_width": float(mean_w),
            "err_width":  float(err_w),

            # flag per sapere se l'errore e' da bootstrap o SEM
            "err_method": (
                "bootstrap"
                if n_boot > 1 and n_entries <= max_boot
                else "sem"
            ),
        }

        results.append(row)

    # =====================================================
    # SAVE PKL FIRST
    # =====================================================

    df = pd.DataFrame(results)

    if save_to:

        tmp_pkl = save_to + ".tmp"

        df.to_pickle(tmp_pkl)

        os.replace(tmp_pkl, save_to)

        print(f"Saved PKL -> {save_to}")

    # =====================================================
    # SAVE ROOT AFTER PKL
    # =====================================================

    if args.save_root and save_to:

        root_output = save_to.replace(
            ".pkl",
            "_histo.root"
        )

        save_theta_histograms_root(
            bins_pre=bins_pre,
            bins_post=bins_post,
            theta_edges=theta_edges,
            output_root=root_output
        )

    return df


# =========================================================
# RUN
# =========================================================

process_x_slice_slim(
    slice_file=args.input,
    save_to=args.output,
    chunk_size_mb=args.chunk_size,
    plane_index=args.plane,
    x_min=args.x_min,
    x_max=args.x_max
)
