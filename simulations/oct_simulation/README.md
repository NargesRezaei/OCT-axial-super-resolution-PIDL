# OCT A-scan Simulation Module

This folder contains the OCT forward model used to convert the multilayer Fresnel reflection spectrum into a physically consistent OCT A-scan.  
It includes two OCT simulators with different levels of spectral sampling, a clean reference-arm model, and an FFT-based reconstruction module.

---

## Folder structure

Files in this folder:
```
oct_simulation/
│
├── oct_forward_simulator.m          # Main OCT simulator (current version)
├── Reference_Mirror.m               # Clean physical model of the reference arm
├── OCT_Analyse.m                    # FFT-based reconstruction to A-scan
├── README.md                        # This file
└── versions/
    ├── oct_forward_simulator_v02.m  # Older / baseline simulator
```
---

## 1. OCT A-scan simulation scripts

Two OCT forward simulators are provided, differing in spectral detail and computational cost.

### Script 1 – OCT A-scan simulator  
*(oct_forward_simulator.m – main script)*

- Uses a moderate spectral bandwidth and fewer wavelength samples.  
- Uses a simpler source and spectrometer model.  
- Designed as a **clean, easy-to-follow implementation** for testing:
  - The multilayer Fresnel model
  - The OCT data pipeline end-to-end  
- Faster and lighter to run:
  - Good for debugging
  - Good for generating quick example A-scans  

> **Summary:** Baseline, simple, and lightweight.

---

### Script 2 – High-resolution wide-band OCT simulator  
*(oct_forward_simulator_v02.m, in `versions/`)*

- Uses a **wider spectral bandwidth**.  
- Uses **more wavelength samples**, leading to:
  - Denser sampling in the spectrum
  - **Higher axial resolution** in the A-scan  
- Employs an updated OCT and SBW modelling approach:
  - More realistic source and spectrometer behaviour
  - Better suited for **fine-scale, realistic simulations**  
- Computationally heavier, but more accurate and detailed.

> **Summary:** More wavelengths + wider spectrum → **higher resolution and more realism**, at higher computational cost.

---

## 2. Reference arm model – `Reference_Mirror.m`

`Reference_Mirror` implements the OCT reference arm as a **single physical mirror** at distance `z0` in a medium with refractive index `N0`.

### Physical behaviour

- Adds a well-defined **round-trip phase term**:  
  exp(−i · 2k · N₀ · z₀ · cos(θ₀))  
- Includes a **controllable reflectance** via `Rcoeff` (0…1).  
- Uses a physically meaningful sign for the reflection coefficient:
  - `mode = 'ideal'` → always `r = −√Rcoeff` (π phase shift)  
  - `mode = 'pec'`   → `r = −√Rcoeff` for s-polarization, `+√Rcoeff` for p-polarization  

### Coupler consistency

- The input spectrum is split by **1/√2** (`Ein = AmpSpectrum / √2`),  
  consistent with a 50/50 coupler model.

### Comparison to old approach

Previously, the reference arm was approximated as a fake two-layer multilayer stack `[N0, N0]` passed through the general Fresnel solver.  
`Reference_Mirror` improves upon this by being:

- **Simpler and more explicit** (no artificial multilayer stack).  
- Free from unnecessary numerical overhead.  
- Easier to control and interpret (only `z0`, `Rcoeff`, and `mode`).  

> **In short:** `Reference_Mirror` provides a clean, physically consistent reference-arm model that replaces the previous “fake [N0, N0] multilayer” construction.

---

## 3. Spectral processing and FFT – `OCT_Analyse.m`

`OCT_Analyse` takes the spectral interferogram and reference signal and reconstructs the depth-resolved OCT A-scan.

It performs:

1. **Background removal**  
   - Removes DC / background components from the interferogram.  

2. **Zero-padding**  
   - Extends the spectrum to improve FFT sampling in depth.  

3. **Resampling from λ to k**  
   - Converts the spectrum from **wavelength domain** to a **linear k-grid** (wavenumber domain), which is required for standard OCT reconstruction.  

4. **FFT-based reconstruction**  
   - Applies a 1D FFT on the linear-k interferogram.  
   - Produces the complex depth profile and the corresponding magnitude A-scan.  

5. **Depth axis generation**  
   - Computes the physical depth axis `Zaxis` associated with the A-scan, based on:
     - Sampling in k-space  
     - Refractive index assumptions  
     - FFT scaling  

> **Output:** A calibrated A-scan and `Zaxis` that can be directly plotted or used as input for dataset generation and deep learning.

---

## Recommended usage

- Use `oct_forward_simulator.m` for:
  - High-resolution simulations  
  - Generating high-quality ground-truth A-scans  
  - Realistic OCT signal modelling  

- Use `oct_forward_simulator_v02.m` (in `versions/`) for:
  - Quick tests  
  - Debugging the Fresnel module  
  - Lightweight experiments when speed matters more than ultimate resolution  

`Reference_Mirror.m` and `OCT_Analyse.m` are shared components used by both simulators to ensure a consistent, physically meaningful OCT forward model.

