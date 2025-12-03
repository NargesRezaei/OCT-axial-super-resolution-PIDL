# Cluster-based Retinal B-scan Simulation

This folder contains the **cluster-based** retinal B-scan simulator.

It combines:

- Real retinal geometry extracted from a segmented OCT image.
- An anatomically inspired **cluster-based refractive-index model** inside each layer.
- A modern OCT forward model using `General_Multilayer_V11` and `Reference_Mirror`.

The goal is to generate synthetic B-scans that follow the **true layer thicknesses** from segmentation, while adding realistic **clustered scattering structure** within each retinal layer.

---

## Overview

The main script (e.g. `bscan_cluster_based.m`) performs the following steps:

### 1. Load segmented retina

A `.mat` file containing segmentation data is loaded, for example:

- `bscn.nI` is a 2D map where non-zero values indicate retinal tissue.
- Each column corresponds to one A-scan.

For each column `jj`:

- Non-zero indices `q` are found along depth.
- `diff(q)` is used to derive contiguous layer thicknesses.
- An `offset` is applied to start slightly above the first layer.
- Thickness values are converted from the original image pixel size (`pixel_size`) to a finer modelling grid with step `n_pixel_size`.

The result is stored in `peak{jj}` as a vector of layer thicknesses (in fine pixels).

---

### 2. Cluster-based refractive-index model per layer

For each A-scan:

- `pxl_ranges = peak{jj}` gives the number of fine pixels per layer.
- A subset of baseline refractive indices `N_1` is chosen according to how many layers are present.
- Layer-dependent parameters are defined, for example:

  - `DN_cluster` – extra index contrast inside clusters.  
  - `dn_all` – background index variation within the layer.  
  - `p_cluster_full` – probability of forming a cluster in each layer.  
  - `cluster_size_range_full` – allowed cluster size range.  
  - `env_after_cluster_range_full` – background length between clusters.  
  - `bright_spot_layers`, `bright_spot_prob_layers` – occasionally add stronger scattering “bright spots” in certain layers (e.g. GCL, IPL, OPL).

For each layer `i`:

- A transition region at the top of the layer can be smoothed by interpolating between the previous layer index and the current layer index.
- Then the rest of the layer is filled using a **cluster process**:

  - With probability `p_cluster`, a cluster is created:
    - A random cluster length is chosen within `cluster_range`.
    - A cluster refractive index is set around `base_value` plus `DN_cluster`, with some random jitter.
    - If `base_value` belongs to one of the `bright_spot_layers`, a bright cluster may be created with probability `bright_spot_prob_layers`.
  - After each cluster, one or more environment pixels are inserted with index near `base_value`.
  - If no cluster is formed, a run of environment pixels of random length within `env_range` is created.

Thickness and index for each layer are stored in:

- `thick_cells{i}` – all equal to `n_pixel_size`.
- `refidx_cells{i}` – the cluster-based refractive-index profile.

All layers are concatenated into:

- `LD` – full thickness array (fine grid).
- `LN` – full refractive-index array, plus a final background value.

---

### 3. OCT forward simulation

The OCT simulation uses:

- A source centered at `lambda_c = 840 nm` with bandwidth `Bandwidth = 60 nm`.
- A spectrometer model defined by `d_sbw`, `lambda_sbw` and a spectral window matrix `sbw_convertor`.
- A physical reference arm field computed with `Reference_Mirror` using parameters `N0`, `z0_R`, `theta0`, `pol`, and `AmpSpectrum`.

For the sample arm:

- `General_Multilayer_Fresnel_V11` is used to compute the reflected field `Er`.

The interferogram is formed as:

- `E_sum = Er + ErR`  
- `I_OCT2 = 0.5 * (E_sum .* conj(E_sum))`

Then the spectrometer / SBW model is applied:

- `I_OCT3 = sbw_convertor * I_OCT2'` (spectral-domain measurement)  
- A reference intensity `ref2` is computed and passed through the same SBW model to obtain `ref3`.

---

### 4. Reconstruction and B-scan formation

The A-scan is reconstructed using:

- `OCT_Analyse(I_OCT3, ref3, lambda_sbw, LD)`

The output depth axis `Xaxis` and intensity `Depth` are:

- Rounded or truncated as needed.
- Stored as one column of the B-scan:

  - `B_Depth(:, jj) = Depth`.

After iterating over all columns, `B_Depth` becomes a 2D B-scan image.

A log-compressed representation is computed for visualization:

- `bscan_db = 10 * log10(B_Depth + 1e-5)`

which can be displayed.

---

## Output

The cluster-based B-scan simulator provides:

- **`B_Depth`** – linear-scale B-scan intensity.  
- **`bscan_db`** – log-compressed B-scan for display.  
- Optionally, **layer-wise OPL information** (e.g. `OPL_2`) can be stored per A-scan for ground-truth boundary positions.

---

## Use cases

This model is useful for:

- Generating anatomically realistic synthetic retinal B-scans.  
- Creating paired data: **B-scan + ground-truth layer boundaries**.  
- Training and evaluating deep-learning models for segmentation, denoising, or super-resolution.  
- Comparing against simpler pixel-based models to study the effect of intra-layer clustering on OCT appearance.
