# Retina Model Simulation

This folder contains retina-specific OCT simulation code used to **generate synthetic datasets for deep learning**.

It provides:

- **Random retinal A-scans** (cluster-based & pixel-based)
- **Retinal B-scan simulation** using real segmentation-derived geometry

These simulations are used to create paired data (signal + layer geometry) for training and evaluating neural networks.

---

## Structure

```
retina_model_simulation/
│
├── ascans_random/
│   ├── cluster_based/   # Anatomically inspired cluster model for A-scans
│   ├── pixel_based/     # Simpler per-pixel random A-scan model
│   └── README.md
│
├── bscan_simulation/
│   ├── cluster_based/   # B-scans with intra-layer clusters and bright spots
│   ├── pixel_based/     # B-scans with per-pixel random refractive index
│   └── README.md
│
├── sample_outputs/        
│   ├── pixel_based_example.png
│   └── cluster_based_example.png
│
└── README.md
```

All models use the core OCT physics implemented in the main ‍‍‍‍`simulations/` module
(Fresnel multilayer + OCT forward simulation + FFT reconstruction).
