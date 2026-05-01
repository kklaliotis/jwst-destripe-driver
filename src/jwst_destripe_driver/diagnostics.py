#!/usr/bin/env python3
"""Generate diagnostic figures from JWST destriping outputs.

This script is a standalone version of ``Diagnostic_Figs.ipynb``. It can
process either a single ``*_crf.fits`` file or every matching file in a
directory. In directory mode it also writes ``med_scatter_stats.pkl`` and the
summary histogram ``RRM_sigmas.png``.
"""

from __future__ import annotations

import argparse
import pickle
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import scipy
from astropy.io import fits
from scipy.fft import fft


NROWS = 2048
AMPCOLS = 512

cmap = plt.get_cmap("inferno")
dark_color = cmap(0.15)
mid_color = cmap(0.55)
light_color = cmap(0.9)


def open_fits(name, i="SCI"):
    with fits.open(name) as hdul:
        return np.copy(hdul[i].data)


def image_fromparams(allparams, i, n_rows=NROWS):
    """Make an image of stripes from the best-fit parameter vector.
    
    Parameters
    ----------
    allparams : 2D array
        Parameters table, shape (n_scas, n_params_per_sca)
    i : int
        Row index (SCA index)
    n_rows : int
        Number of rows in the image (default 2048 for JWST)
    Returns
    -------
    2D array
        Stripe image from parameters
    """
    allparams = np.asarray(allparams)
    row_image = allparams[i, :].reshape(-1, 1)

    if AMPCOLS is None or AMPCOLS <= 0:
        # Row-only model
        return np.tile(row_image, (1, NROWS))
    else:
        # Row + column-block model
        n_cols = NROWS
        n_col_blocks = n_cols // AMPCOLS
        params_per_row = 1  # Assuming constant model; adjust if linear
        n_row_params = n_rows * params_per_row
        
        # Extract row and column parts
        row_params = allparams[i, :n_row_params].reshape(-1, 1)
        col_params = allparams[i, n_row_params:n_row_params + n_col_blocks]
        
        # Build row image
        row_image = np.tile(row_params, (1, NROWS))
        
        # Build column-block image
        col_image = np.zeros((NROWS, NROWS))
        for b in range(n_col_blocks):
            j_start = b * AMPCOLS
            j_end = j_start + AMPCOLS
            col_image[:, j_start:j_end] = col_params[b]
        
        return row_image + col_image


def get_row_medians(image):
    """Get the median value in each row."""

    x_min = 0
    x_max = image.shape[1]

    medians = []
    for row in range(image.shape[0]):
        medians.append(np.median(image[row, x_min:x_max]))

    return np.array(medians)


def power_spectrum_1d(medians):
    """Calculate the 1D power spectrum of the row medians."""

    fft_median = fft(medians - np.mean(medians))
    return np.abs(fft_median) ** 2


def power_spectrum_2d(image, scale=1):
    im_ft = np.abs(scipy.fft.fft2(image * scale))
    im_ft = scipy.fft.fftshift(im_ft)
    return im_ft**2 / (NROWS**2)


def plot_three_images(
    img1,
    img2,
    img3,
    titles=("Image 1", "Image 2", "Image 3"),
    suptitle="Suptitle",
    cmap="inferno",
    cbar_label="Flux [DN/fr]",
    savefig=None,
    pctile=(20, 80),
):
    fig, axes = plt.subplots(1, 3, figsize=(12, 6), constrained_layout=True)
    ims = []
    for ax, img, title in zip(axes, [img1, img2, img3], titles):
        im = ax.imshow(
            img,
            cmap=cmap,
            origin="lower",
            vmin=np.percentile(img, pctile[0]),
            vmax=np.percentile(img, pctile[1]),
        )
        ax.set_title(title, fontsize=12)
        ax.axis("off")
        ims.append(im)

    for ax, im in zip(axes, ims):
        cbar = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.048, pad=0.07)
        cbar.set_label(cbar_label, fontsize=12)

    plt.suptitle(suptitle, y=0.9, x=0.33)
    if savefig is not None:
        fig.savefig(savefig, bbox_inches="tight")
    plt.close(fig)


def row_median_resid_scatterhist(med_pre, med_post, savefig, extrareturn=False):
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.005

    avg_pre = np.mean(med_pre)
    avg_post = np.mean(med_post)

    rect_scatter = [left, bottom, width, height]
    rect_histy = [left + width + spacing, bottom, 0.2, height]
    fig = plt.figure(figsize=(10, 5))

    ax = fig.add_axes(rect_scatter)
    ax_histy = fig.add_axes(rect_histy, sharey=ax)

    ax_histy.tick_params(axis="y", labelbottom=False)

    x = np.arange(0, 2048)
    y = avg_pre - med_pre[0:2048]
    x2 = np.arange(0, 2048)
    y2 = avg_post - med_post[0:2048]
    ax.scatter(x=x, y=y, marker="+", color=dark_color, label="Pre", alpha=1)
    ax.scatter(x=x2, y=y2, marker="+", color=mid_color, label="Post", alpha=0.8)
    ax.legend()
    ax.set_xlabel("Row Index", fontsize=12)
    ax.set_ylabel(r"Residual Row Median [e]", fontsize=12)
    ax.set_ylim((np.min(y), np.max(y)))
    ax.axhline(0, color="black", linestyle="dashed", alpha=0.5)

    binwidth = 0.25
    xymax = max(np.max(np.abs(x)), np.max(np.abs(y)))
    lim = (int(xymax / binwidth) + 1) * binwidth

    bins = np.arange(-lim, lim + binwidth, binwidth)
    ax_histy.hist(
        y,
        bins=bins,
        orientation="horizontal",
        histtype="step",
        color=dark_color,
        label=f"Pre\n$\sigma$={np.std(y):.3}",
    )
    ax_histy.hist(
        y2,
        bins=bins,
        orientation="horizontal",
        histtype="step",
        color=mid_color,
        label=f"Post\n$\sigma$={np.std(y2):.3}",
    )
    ax_histy.legend(fontsize="small")
    if savefig is not None:
        fig.savefig(savefig, bbox_inches="tight")
    if extrareturn:
        plt.close(fig)
        return y, y2
    plt.close(fig)


def adaptive_bin(x, y, bin_edges=(10, 100, 1000, 10000), bin_factors=(4, 40, 400, 400)):
    """Adaptive binning of (x, y)."""

    x_binned, y_binned = [], []
    n = len(y)
    idx = 0

    first_chunk = min(bin_edges[0], n)
    x_binned.extend(x[:first_chunk])
    y_binned.extend(y[:first_chunk])
    idx = first_chunk

    factor_idx = 0
    while idx < n:
        if factor_idx < len(bin_factors):
            factor = bin_factors[factor_idx]
        else:
            factor *= 10

        chunk = slice(idx, min(idx + factor, n))
        if chunk.stop <= chunk.start:
            break
        x_binned.append(np.mean(x[chunk]))
        y_binned.append(np.mean(y[chunk]))
        idx += factor

        if factor_idx < len(bin_edges) - 1 and idx >= bin_edges[factor_idx + 1]:
            factor_idx += 1

    return np.array(x_binned), np.array(y_binned)


def stem_from_crf_path(path):
    name = Path(path).name
    if name.startswith("jw") and name.endswith("_crf.fits"):
        return name[2:-9]
    if name.endswith("_crf.fits"):
        return name[:-9]
    return Path(path).stem


def load_params_table(params_path):
    with fits.open(params_path) as hdul:
        return np.asarray(hdul[0].data)


def get_param_row(params_table, index):
    params = np.asarray(params_table)
    if params.ndim == 1:
        return params
    if index is None:
        index = 0
    if index < 0 or index >= params.shape[0]:
        raise IndexError(
            f"Parameter index {index} is out of range for table with {params.shape[0]} rows."
        )
    return params[index]


def load_mask(mask_path, shape):
    if mask_path.exists():
        mask = open_fits(mask_path, i=0)
        if mask.shape != shape:
            raise ValueError(f"Mask shape {mask.shape} does not match image shape {shape} for {mask_path}")
        return mask
    return np.ones(shape, dtype=float)


def save_noise_distribution_plot(pre_medians, post_medians, save_path):
    y1 = pre_medians - np.median(pre_medians)
    y2 = post_medians - np.median(post_medians)

    fig = plt.figure(figsize=(7, 5))
    plt.hist(y1, bins=100, histtype="step", density=True, label="Before stripe removal", color=dark_color)
    plt.hist(y2, bins=100, histtype="step", density=True, label="After stripe removal", color=mid_color)
    plt.xlabel("Noise value (median-centered)")
    plt.ylabel("Normalized count")
    plt.title("Noise distributions (clipped to 0.5-99.5%)")
    plt.legend()
    fig.savefig(save_path, bbox_inches="tight")
    plt.close(fig)


def save_power_spectrum_change_plot(noisy_image, ds_image, mask, file_stem, save_path):
    noise_initial = noisy_image * mask
    noise_final = ds_image * mask

    fft_noise = scipy.fft.fft2(noise_initial - np.median(noise_initial))
    pspec_noise = np.abs(fft_noise) ** 2

    fft_noise_post = scipy.fft.fft2((noise_final - np.median(noise_final)) * mask)
    pspec_noise_post = np.abs(fft_noise_post) ** 2

    pspec_noise = np.ravel(pspec_noise.T)
    pspec_noise_post = np.ravel(pspec_noise_post.T)

    freqs = np.arange(len(pspec_noise))
    freqs_binned, pspec_noise_binned = adaptive_bin(freqs, pspec_noise / len(pspec_noise))
    freqs_post_binned, pspec_noise_post_binned = adaptive_bin(freqs, pspec_noise_post / len(pspec_noise_post))

    left, width = 0.1, 0.7
    bottom, height = 0.2, 0.8

    rect_ps = [left, bottom, width, height]
    rect_frac = [left, bottom - 0.4, 0.7, 0.3]
    fig = plt.figure(figsize=(10, 5))

    ax = fig.add_axes(rect_ps)
    ax_frac = fig.add_axes(rect_frac, sharex=ax)
    ax_frac.tick_params(axis="y", labelbottom=False)

    ax.plot(freqs_binned, pspec_noise_binned, label="Pre-Destripe", color=dark_color, marker=".")
    ax.plot(freqs_post_binned, pspec_noise_post_binned, label="Post-Destripe", color=mid_color, alpha=0.8, marker=".")
    ax.legend(loc="center right")
    ax.set_ylabel(r"Power [$e^2$ px]", fontsize=12)

    ax_frac.plot(
        freqs_post_binned,
        pspec_noise_post_binned / pspec_noise_binned,
        label="Fractional Change",
        color=dark_color,
        marker=".",
    )
    ax_frac.axhline(y=1, c="k", ls="--", label="No change (y=1)")

    plt.suptitle(f"Power Spectrum of JW{file_stem}", fontsize=14)
    plt.xlabel("Frequency [cycles/SCA]", fontsize=12)
    plt.ylabel(r"$P_{post}/P_{pre}$", fontsize=12)
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    ax_frac.grid(True, which="both", linestyle="--", linewidth=0.5)
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax_frac.set_yscale("log")
    ax_frac.set_xscale("log")
    ax.set_xlim(0.9, 2048**2 / 2)
    plt.legend(loc="best")
    fig.savefig(save_path, bbox_inches="tight")
    plt.close(fig)


def save_rrm_sigma_summary(summary_data, save_path):
    med_scatter_state = np.array(summary_data)

    fig = plt.figure(figsize=(7, 5))
    plt.hist(
        med_scatter_state[0, :],
        color=dark_color,
        label="Pre-Destripe",
        bins=30,
        alpha=0.6,
    )
    plt.axvline(
        sig_pre := np.mean(med_scatter_state[0, :]),
        linestyle="--",
        color=dark_color,
        label=r"$\mu=$" + f"{sig_pre:.3f}",
    )
    plt.hist(
        med_scatter_state[1, :],
        color=mid_color,
        label="Post-Destripe",
        bins=30,
        alpha=0.6,
    )
    plt.axvline(
        sig_post := np.mean(med_scatter_state[1, :]),
        linestyle="--",
        color=mid_color,
        label=r"$\mu=$" + f"{sig_post:.3f}",
    )
    plt.title(r"Distribution of Row Median Residual $\sigma$")
    plt.xlabel(r"$\sigma_{RRM}$ [e]")
    plt.ylabel("Counts")
    plt.xlim(0, 3)
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.legend()
    fig.savefig(save_path, bbox_inches="tight")
    plt.close(fig)


def make_figures_for_file(crf_path, params_table, fig_dir, mask_dir, param_index=None):
    crf_path = Path(crf_path)
    file_stem = stem_from_crf_path(crf_path)
    noisy_image = open_fits(crf_path)
    noisy_image = np.where(np.isnan(noisy_image), 0, noisy_image)

    params_table = np.asarray(params_table)
    row_index = 0 if params_table.ndim == 1 else (param_index if param_index is not None else 0)
    stripe_image = image_fromparams(params_table, i=row_index)
    ds_image = noisy_image - stripe_image

    mask_path = Path(mask_dir) / f"{file_stem}_mask.fits"
    mask = load_mask(mask_path, noisy_image.shape)

    pre_centered = noisy_image - np.median(noisy_image)
    post_centered = ds_image - np.median(ds_image)
    pre_medians = get_row_medians(pre_centered)
    post_medians = get_row_medians(post_centered)

    fig_dir = Path(fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)

    plot_three_images(
        pre_centered,
        post_centered,
        stripe_image,
        titles=("Pre-Destripe", "Post-Destripe", "Difference"),
        suptitle=f"Image: {file_stem}",
        cbar_label="MJy/sr",
        savefig=fig_dir / f"triptych_unmasked_{file_stem}.png",
    )

    plot_three_images(
        pre_centered * mask,
        post_centered * mask,
        stripe_image,
        titles=("Pre-Destripe", "Post-Destripe", "Difference"),
        suptitle=f"Masked Image: {file_stem}",
        cbar_label="MJy/sr",
        savefig=fig_dir / f"triptych_masked_{file_stem}.png",
    )

    row_median_resid_scatterhist(
        pre_medians,
        post_medians,
        savefig=fig_dir / f"row_median_residuals_{file_stem}.png",
    )

    save_noise_distribution_plot(
        pre_medians,
        post_medians,
        fig_dir / f"row_median_distribution_{file_stem}.png",
    )

    save_power_spectrum_change_plot(
        noisy_image,
        ds_image,
        mask,
        file_stem,
        fig_dir / f"power_spectrum_change_{file_stem}.png",
    )

    return np.std(pre_medians), np.std(post_medians)


def build_med_scatter_stats(stats, save_path):
    save_path = Path(save_path)
    with save_path.open("wb") as handle:
        pickle.dump([list(stats[0]), list(stats[1])], handle)


def parse_sca_list_from_logfile(logfile_path):
    """Parse SCA list from destripe .out logfile.
    
    Expected format:
        SCA 0: jw03707032001_04101_00001_nrcb1
        SCA 1: jw03707032001_04101_00001_nrcb2
        ...
    
    Returns a dict: {file_stem -> sca_index}
    """
    sca_map = {}
    try:
        with open(logfile_path, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith("SCA "):
                    parts = line.split(": ", 1)
                    if len(parts) == 2:
                        sca_num_str = parts[0].replace("SCA ", "").strip()
                        sca_num = int(sca_num_str)
                        file_path = parts[1].strip()
                        file_stem = stem_from_crf_path(file_path)
                        sca_map[file_stem] = sca_num
    except (OSError, ValueError) as e:
        print(f"Warning: Could not parse SCA list from {logfile_path}: {e}")
    print(sca_map)
    return sca_map


def iter_crf_files(target):
    target = Path(target)
    if target.is_file():
        return [target]
    return sorted(target.glob("*_crf.fits"))


def main(argv=None):
    parser = argparse.ArgumentParser(description="Generate JWST destriping diagnostic figures.")
    parser.add_argument("target", help="A single *_crf.fits file or a directory containing them.")
    parser.add_argument(
        "--output-dir",
        default=".",
        help="Directory containing final_params.fits and masks/ (default: current directory).",
    )
    parser.add_argument(
        "--params",
        default=None,
        help="Path to final_params.fits (default: <output-dir>/final_params.fits).",
    )
    parser.add_argument(
        "--mask-dir",
        default=None,
        help="Directory containing <stem>_mask.fits files (default: <output-dir>/masks).",
    )
    parser.add_argument(
        "--fig-dir",
        default=None,
        help="Directory for saved figures (default: <output-dir>/figures).",
    )
    parser.add_argument(
        "--param-index",
        type=int,
        default=None,
        help="Row index in final_params.fits for single-file processing or to override directory order.",
    )
    parser.add_argument(
        "--logfile",
        default=None,
        help="Name of .out logfile to parse SCA list for automatic index lookup. Full path or just filename if in output-dir.",
    )
    args = parser.parse_args(argv)

    target = Path(args.target)
    output_dir = Path(args.output_dir)
    params_path = Path(args.params) if args.params is not None else output_dir / "final_params.fits"
    mask_dir = Path(args.mask_dir) if args.mask_dir is not None else output_dir / "masks"
    fig_dir = Path(args.fig_dir) if args.fig_dir is not None else output_dir / "figures"
    log_file = output_dir / str(args.logfile) 
    print(f"Using parameters from: {params_path}")
    print(f"Using logfile: {log_file}")

    params_table = load_params_table(params_path)
    files = iter_crf_files(target)
    if not files:
        raise FileNotFoundError(f"No *_crf.fits files found under {target}")

    fig_dir.mkdir(parents=True, exist_ok=True)

    # Parse SCA list if logfile provided
    sca_map = {}
    if log_file:
        sca_map = parse_sca_list_from_logfile(log_file)

    pre_sigmas = []
    post_sigmas = []

    if target.is_file():
        # Determine index: explicit param-index > logfile lookup > default 0
        if args.param_index is not None:
            single_index = args.param_index
        else:
            file_stem = stem_from_crf_path(target)
            single_index = sca_map.get(file_stem, 0)
        pre_sigma, post_sigma = make_figures_for_file(
            target,
            params_table,
            fig_dir=fig_dir,
            mask_dir=mask_dir,
            param_index=single_index,
        )
        pre_sigmas.append(pre_sigma)
        post_sigmas.append(post_sigma)
    else:
        for file_idx, crf_path in enumerate(files):
            # Determine index: explicit param-index > logfile lookup > enumerate order
            if args.param_index is not None:
                row_index = args.param_index
            else:
                file_stem = stem_from_crf_path(crf_path)
                row_index = sca_map.get(file_stem, file_idx)  # Default to file order
            pre_sigma, post_sigma = make_figures_for_file(
                crf_path,
                params_table,
                fig_dir=fig_dir,
                mask_dir=mask_dir,
                param_index=row_index,
            )
            pre_sigmas.append(pre_sigma)
            post_sigmas.append(post_sigma)

    if target.is_dir():
        summary_path = fig_dir / "med_scatter_stats.pkl"
        build_med_scatter_stats((pre_sigmas, post_sigmas), summary_path)
        save_rrm_sigma_summary((pre_sigmas, post_sigmas), fig_dir / "RRM_sigmas.png")


if __name__ == "__main__":
    main()