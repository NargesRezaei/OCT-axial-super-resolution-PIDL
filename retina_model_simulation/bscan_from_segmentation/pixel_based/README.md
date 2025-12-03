# Pixel-based Retinal B-scan Simulation

This folder contains the **pixel-based** retinal B-scan simulator.  
In this model, each A-scan inside the B-scan is generated directly from the **segmented retinal geometry** and a simple **per-layer random refractive-index model**.

It provides a clean baseline for generating synthetic B-scans that follow real retinal layer thicknesses but without internal cluster structure.

---

## Overview

The main script (e.g., `bscan_pixel_based.m`) performs the following steps:

### 1. Load segmented retina
A `.mat` file containing a segmentation-like map is loaded:

- `bscn.nI` contains non-zero values where retinal tissue exists.
- Each column corresponds to an A-scan.
- For each column, consecutive non-zero segments define the thickness of each anatomical layer.

These thickness profiles are extracted into `peak{jj}`.

---

### 2. Build per-layer refractive-index profile
For each A-scan:

- `pxl_ranges` defines the number of pixels per layer (from segmentation).
- A subset of predefined refractive indices `N_1` is selected.
- Each layer is assigned a refractive index using:

  - A base refractive index for the layer  
  - A random deviation controlled by `dn_layer(i)`  
  - An optional transition zone between neighboring layers for smoothness

The results are stored in:

- `thick_cells{i}` – thickness values for each pixel  
- `refidx_cells{i}` – refractive index per pixel

These are concatenated to form:

- `LD` – full thickness array  
- `LN` – full refractive-index profile

---

### 3. OCT forward simulation
The OCT signal is generated using:

- `General_Multilayer_Fresnel_V10` for both sample and reference arms  
- An interferogram computed as:
 `E_sum  = Er + ErR_2`
 `I_OCT2 = smooth(0.5 * (E_sum .* conj(E_sum)))`

A spectrometer sampling model and SBW (`sbw_convertor`) are then applied to obtain the spectrally filtered signal `I_OCT3`.

Finally, reconstruction is done using:

- ` [Xaxis, Depth] = OCT_Analyse(I_OCT3, ref3, lambda_sbw , LD);`

Each reconstructed A-scan (`Depth`) is stacked into the B-scan matrix:

- `B_Depth(:, jj)`

Ground-truth cumulative optical-path-length boundaries for that column are stored as:

- `GT{jj}`

---

## Output

The model produces:

- **`B_Depth`** – the simulated B-scan (linear intensity)  
- **`bscan_db`** – log-compressed B-scan for visualization  
- **`GT{jj}`** – layer boundary ground-truth for each A-scan  

A B-scan image may also be saved (for example: `b-Scan_out.png`).

---

## Use cases

The pixel-based B-scan is suitable for:

- Baseline dataset generation  
- Reproducing older simulation styles  
- Debugging OCT forward models  
- Comparing against more advanced cluster-based B-scan models  

It provides a simple and stable way to convert segmentation-based thickness profiles into synthetic OCT B-scans.

