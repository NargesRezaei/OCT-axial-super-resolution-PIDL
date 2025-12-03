# Random Retinal A-scan Simulation

This folder contains scripts for generating **synthetic retinal A-scans** using two different modelling strategies:

1. A **cluster-based model** that builds each layer from refractive-index clusters and background pixels.
2. A **pixel-based model** that assigns a refractive index to each depth pixel directly.

Both models use the same OCT forward-physics implemented in the `simulations/` module (Fresnel + OCT forward + analysis) and only differ in **how the retinal refractive-index profile is generated**.

---

## Folder structure

- `cluster_based/`  
  Cluster-based random A-scan generation (newer and richer structural model).

- `pixel_based/`  
  Pixel-based random A-scan generation (older, more direct per-pixel model).

Each subfolder contains its own main MATLAB script and a README describing the details of that approach.

---

## Dependencies

These scripts **do not re-implement** the OCT physics.  
They call the core functions from the main simulation modules, for example:

- `General_Multilayer_Fresnel_V11` or `General_Multilayer_V10`
- `Reference_Mirror`
- `OCT_Analyse`

Make sure the `simulations/` folder is on the MATLAB path before running:

```
addpath('../../simulations');
addpath('../../simulations/multilayer_fresnel');
addpath('../../simulations/oct_simulation');
```

---

## Recommended usage

**Use the *cluster-based model* when you want:**
- Layer-specific structure
- Random clusters and bright spots
- More realistic retinal texture along depth

**Use the *pixel-based model* for:**
- Reproducing older experiments
- Comparing with the original, simpler random model
- Cases where per-pixel randomness is preferred

