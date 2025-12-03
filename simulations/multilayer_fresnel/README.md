# Multilayer Fresnel Solver (V7–V11)

This folder contains several generations of a multilayer Fresnel solver for simulating reflection and transmission in stratified media.  
The focus is on evolving from a simple loop-based implementation to a fully vectorized, multi-wavelength solver suitable for OCT simulations.

The **current final implementation** is V11.  
Earlier versions (V7, V71, V10) are preserved in the `versions/` folder to document the development history.

---

## Folder structure

```
multilayer_fresnel/
│
├── General_Multilayer_Fresnel_V11.m     # Final multilayer Fresnel implementation
├── transfer_V7_vectorized.m             # Vectorized interface matrices (helper)
├── prop_V7_vectorized_fullspec.m        # Vectorized propagation matrices (helper)
├── README.md                            
└── versions/
    ├── General_Multilayer_Fresnel_V7.m  # Baseline loop-based version
    ├── General_Multilayer_Fresnel_V71.m # Partially vectorized version
    ├── General_Multilayer_Fresnel_V10.m # First multi-wavelength cell-array version
```

---
## Final implementation – V11

**Main function**
- `General_Multilayer_Fresnel_V11.m`

**New helpers (fully numeric, array-based)**

- `transfer_V7_vectorized.m`  
  Returns interface matrices as a **2×2×N numeric array** (`M_ifc`),  
  where **N** is the number of interfaces.  
  Also returns `theta2_all` as the transmitted angle after the last interface.  
  No cell arrays, no replication per wavelength.

- `prop_V7_vectorized_fullspec.m`  
  Returns propagation matrices as a **2×2×N×L numeric array** (`M_prop`),  
  where **N** is the number of layers and **L** is the number of wavelengths.

**Characteristics**

- Fully numeric 3D/4D arrays designed to work with `pagemtimes` (batched matrix multiplication).  
- Multi-wavelength handling is clean and vectorized:
  - Layers are cascaded in a simple loop over layer index.
  - All wavelengths are processed in parallel.
- No `cellfun`, no `reshape` / `mod` hacks.
- This is the **cleanest, fastest, and most scalable** multi-wavelength implementation in the series and is considered the **final version**.

---

## Legacy versions (V7–V10)

The `versions/` folder contains older implementations preserved for reproducibility and documentation.

---

### V7 — baseline loop-based implementation

**Main function**
- `General_Multilayer_Fresnel_V7.m`

**Helpers (conceptual)**
- `transref_V7`
- `prop_V7`

**Characteristics**
- Fully loop-based:
  - Outer loop over wavelengths
  - Inner loop over interfaces/layers
- `transref_V7`:
  - Computes each interface sequentially
  - Updates angle layer-by-layer (physically correct)
  - Returns a **1×(N−1) cell array** of 2×2 interface matrices
- `prop_V7`:
  - Builds **2×2 diagonal propagation matrices** in a 1×N cell array
- Simple and clear, but not optimized for multi-wavelength performance.

---

### V71 — partially vectorized interface handling

**Main function**
- `General_Multilayer_Fresnel_V71.m`

**Helpers**
- `transref_V7_cell`
- `prop_V7` (same as in V7)

**Characteristics**
- Still loops over wavelengths, but interface calculations are partially vectorized.
- `transref_V7_cell`:
  - Computes Fresnel coefficients vectorized
  - Returns a **1×(N−1) cell array** of 2×2 interface matrices
- Propagation remains cell-based.
- Slightly faster and cleaner than V7, but still based on loops and cells.

---

### V10 — first multi-wavelength cell-array version

**Main function**
- `General_Multilayer_Fresnel_V10.m`

**Helpers (older vectorized style)**
- `transref_V7_vectorized(n, theta_in, p, q)`  
  Returns a **q×(N−1) cell array** of 2×2 interface matrices.
- `prop_V7_vectorized`  
  Returns propagation matrices as a **cell array**, one per layer/wavelength.

**Characteristics**
- Handles **all wavelengths in one call**.
- Uses:
  - `repmat` on interface matrices across wavelengths
  - `cellfun`, `reshape`, and `mod(...)` tricks to cascade and regroup matrices
- Faster than loop-based versions but structurally complex due to heavy use of cell arrays.

---

## Recommended usage

Use **V11** for all new experiments and production code:
- `General_Multilayer_Fresnel_V11.m`
- `transfer_V7_vectorized.m`
- `prop_V7_vectorized_fullspec.m`

Use **V7, V71, V10** only if:
- You need to reproduce older results, or
- You want to study the evolution toward a fully vectorized implementation.

---
