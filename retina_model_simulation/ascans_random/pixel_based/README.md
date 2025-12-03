# Pixel-based Random Retinal A-scan Generator

This folder contains the **pixel-based** model for generating synthetic retinal A-scans.

Unlike the cluster-based model, here each layer is represented directly on a **fixed depth grid** with a constant `pixel_size`, and each pixel in a layer is assigned a refractive index with some random deviation.

The model is simpler and more direct, and corresponds to an **earlier version** of the retinal A-scan simulator.

---

## Overview of the method

The main script in this folder:

- Defines a set of layers with baseline refractive indices `N_1` and thickness ranges `Depth_ranges`.
- Converts physical thickness ranges into pixel ranges using `pixel_size`.
- For each sample:
  - Randomly selects the number of layers (e.g. 9–11).
  - For each layer:
    - Samples a random thickness in pixels.
    - Draws a refractive index for each pixel:
      - Around a base value with maximum deviation `dn_layer`.
      - Optionally applies a **smooth transition** region between neighbouring layers to avoid sharp discontinuities.
    - Stores thickness and refractive index profiles in `thick_cells` and `refidx_cells`.
  - Concatenates all layers into:
    - `LD` – pixel-wise thickness array
    - `LN` – pixel-wise refractive index array
  - Appends the final background index at the end.

This produces a **fine-grained refractive-index profile** along depth, purely in a pixel-wise fashion.

---

## OCT simulation and analysis

The pixel-based model uses an older OCT forward configuration:

- Source parameters (e.g. lambda_c, Bandwidth, FWHM) are defined.
- Spectrometer sampling and SBW are modelled via:
  - `d_sbw`, `d_convert`, `sbw_Window`, and `sbw_convertor`.
- The reference arm is generated using `General_Multilayer_Fresnel_V10` with a `[N0, N0]` stack (older mirror model).
- The sample arm is generated via `General_Multilayer_Fresnel_V10(lambda, [N0, LN], [z0 LD], ...)`.
- The interferogram is optionally smoothed and passed through a spectrometer model.
- `OCT_Analyse` is used to convert the spectral interferogram to:
  - `Xaxis` – depth axis
  - `Depth` – recovered A-scan

An example plot is included in the script that overlays:

- The reconstructed A-scan
- Vertical lines at cumulative `OPL_2` to mark layer boundaries

Saving to disk is currently commented out but can be enabled (e.g. saving `OPL_2` and `Depth` to `.mat` files).

---

## Recommended usage

This pixel-based model is useful for:

- Reproducing older experiments or baselines  
- Analysing per-pixel refractive-index noise models  
- Comparing the effect of cluster-based vs pixel-based structural assumptions

For large-scale dataset generation and richer structural realism, the **cluster-based model** is usually preferred.
