# OCT A-scan Simulation Module Based on Ray Tracing Matrix 

This folder contains the OCT forward model used to convert the multilayer Fresnel reflection spectrum into a physically consistent OCT A-scan.  
The two simulators use different spectral bandwidths and sampling densities, which makes one of them **more accurate** and **higher-resolution** than the other.

---

## Folder structure

```
oct_simulation_raytracing/
│
├── oct_forward_simulator.m          # Main OCT simulator (efficient, dataset-friendly)
├── Reference_Mirror.m               # Clean physical model of the reference arm
├── OCT_Analyse.m                    # FFT-based reconstruction
├── README.md                        # This file
└── versions/
    ├── oct_forward_simulator_v02.m  # Higher-detail, heavier version
```

---

## 1. OCT A-scan simulation scripts

### Script 1 – OCT A-scan simulator  
*(oct_forward_simulator.m – main script)*

- Efficient and clean implementation  
- Moderate number of wavelength samples  
- Stable and suitable for:
  - Dataset generation  
  - Debugging  
  - Most standard experiments  

> **Balanced choice:** fast, reliable, and accurate for large-scale simulation.

---

### Script 2 – High-detail OCT simulator  
*(oct_forward_simulator_v02.m — in `versions/`)*

- Uses **denser spectral sampling** over the same bandwidth  
- Produces smoother interferograms and slightly more accurate spectral structure  
- Computation time is noticeably higher  

> **Best for:** validation and precision-focused experiments.

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

- **Use `oct_forward_simulator.m`** for most applications:  
  - Dataset generation  
  - Training pipelines  
  - Large-scale simulations  
  - General OCT modelling  

  It is **computationally efficient** while preserving the same axial resolution as the high-detail version.

---

- **Use `oct_forward_simulator_v02.m`** only when:
  - You need maximum spectral fidelity  
  - You are performing validation or numerical accuracy checks  

  It samples the same bandwidth more densely, giving slightly more detailed interferograms at the cost of speed.


