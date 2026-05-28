# OCT Signal Simulator Using Angular Spectrum

This repository contains a simple MATLAB-based simulator for generating and analyzing Optical Coherence Tomography (OCT) signals from a multilayer sample.

The code models the propagation, reflection, and transmission of a broadband optical wave through multiple layers with different refractive indices and thicknesses. The reflected waves from the sample arm are combined with a reference arm, and the resulting interference spectrum is analyzed to reconstruct depth information.

## Files

| File | Description |
|---|---|
| `OCTSignal.m` | Main script for defining OCT parameters, generating the broadband source spectrum, simulating the sample and reference arms, and plotting the OCT signal and depth profile. |
| `OCTsimulator.m` | Simulates wave propagation inside the multilayer sample and sums all reflected waves returning to the initial medium. |
| `tracer.m` | Traces each transmitted and reflected wave at layer interfaces. It updates the wave direction, current layer, reflection count, and propagation state. |
| `Propagate_ASM.m` | Applies phase propagation to the optical field through a layer using a simplified angular spectrum method. |
| `rforward.m` | Calculates the reflection coefficient at an interface between two media. |
| `tforward.m` | Calculates the transmission coefficient at an interface between two media. |
| `OCT_Analyse.m` | Processes the simulated OCT interference signal, removes the background, resamples the spectrum from wavelength space to wavenumber space, and reconstructs the depth profile. |

## Simulation Overview

The simulation follows these main steps:

1. Define OCT source parameters such as center wavelength, bandwidth, and spectral resolution.
2. Define the multilayer sample using refractive indices and layer thicknesses.
3. Generate a broadband Gaussian source spectrum.
4. Propagate the incident wave through the sample.
5. Calculate reflected and transmitted waves at each interface.
6. Sum all reflected waves returning from the sample arm.
7. Generate the reference arm signal.
8. Combine sample and reference fields to obtain the OCT interference signal.
9. Analyze the interference spectrum to reconstruct the depth profile.

## Main Parameters

Some important parameters in `OCTSignal.m` are:

```matlab
lambda_c  = 850e-9;     % Center wavelength
Bandwidth = 240e-9;     % Source bandwidth
d_lambda  = 0.1e-9;     % Wavelength step

N0 = 1;                 % Refractive index of the surrounding medium
Ni = [1 1.5 1.2 1];     % Refractive indices of the sample layers
Di = [300e-6 250e-6 150e-6]; % Thicknesses of the sample layers

Num_reflections = 4;    % Maximum number of reflections to calculate
