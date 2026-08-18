# Retarded van der Waals Deposition Model

MATLAB code from my PhD (Chemical Engineering, University of Melbourne, 2023) for predicting how colloidal particles, homogeneous, coated, and air-filled capsules, deposit onto a surface while flowing through a microfluidic channel.

Full thesis: *Deposition Behaviour of Coated Particles and Capsules* (Tanweepriya Das, 2023). Related publication: Das, T.; Smith, J. D.; Uddin, M. H.; Dagastine, R. R. "Anisotropic Particle Fabrication Using Thermal Scanning Probe Lithography." *ACS Appl. Mater. Interfaces* 2022, 14(17), 19878–19888. https://pubs.acs.org/doi/10.1021/acsami.2c02885

## The problem this solves

A particle moving through a microfluidic channel gets pushed toward the wall by diffusion and shear, and whether it actually sticks depends on the balance of two surface forces acting over the last few nanometres: van der Waals attraction and electrostatic (double-layer) repulsion. The existing Brownian dynamics model this work builds on (Cejas et al.) simulated that balance for homogeneous polystyrene particles, but assumed the van der Waals interaction could be treated as a **non-retarded, constant-Hamaker** attraction. That assumption breaks down for two reasons this thesis addresses:

1. **Retardation.** At the separation distances that actually govern deposition (a few to tens of nanometres), the electromagnetic propagation delay between the two interacting surfaces measurably weakens the van der Waals attraction. A constant Hamaker function ignores this and systematically overestimates attraction — the thesis shows this leads to ~67% overestimation of deposition counts for 1 μm silica particles in 1 M solution (Chapter 5, Fig. 5.2).
2. **Composite geometry.** A coated particle (e.g. titania coated with silica, or a polymer-shelled air capsule) has two materials with different optical properties contributing to the interaction, and the Hamaker function stops being a single constant — it depends on coating thickness and varies non-monotonically with separation distance (Chapter 5, Fig. 5.4).

This repo is the code that: (a) computes the retarded, composite-aware Hamaker/van der Waals interaction from Lifshitz theory, and (b) runs the resulting force through a Langevin/Brownian-dynamics particle-tracking simulation to predict deposition kinetics — validated in the thesis against a purpose-built experimental technique (Resonance Imaging Microscopy).

## Repo structure

```
hamaker-function/       Retarded Hamaker function calculation for a homogeneous
                         particle-surface pair (the building block everything else extends)

homogeneous-particles/  BD simulation for uncoated particles (silica, titania) —
                         the baseline case, one material only

composite-particles/    BD simulation for coated/composite particles: titania-alumina,
                         titania-silica, silica-polydopamine, silica-polystyrene.
                         This is the core contribution — the van der Waals force here
                         accounts for both core and coating material properties

capsules/                BD simulation for air-filled polystyrene capsules, where the
                         van der Waals interaction is long-range *repulsive* rather
                         than attractive (Chapter 7)

legacy-baseline-model/  The earlier, simpler model this thesis started from and modified
                         — kept here for contrast, not as current code
```

Each material-system folder is a representative sample (one or two concentrations) pulled from a much larger parameter sweep — the full study ran the same simulation across ~10 electrolyte concentrations and multiple particle sizes per system to build the deposition kinetics and phase diagrams reported in the thesis. This repo shows the model itself, not the full sweep.

## How the code fits together

Each material system generally has up to three kinds of file:

- **`fun_*.m`** — the simulation kernel. For a single particle radius, runs `nstatmax` batches of `N` particles through the channel using a Langevin equation: Poiseuille flow profile (solved as a Fourier series across the channel), hindered diffusion coefficients near the wall (`beta_x/y/z`, from lubrication theory), the van der Waals force (piecewise: an analytical near-contact fit, a table lookup interpolated from the theoretical retarded Hamaker curve at intermediate range, and a far-field fit — see Chapter 5, Fig. 5.5 for why a single functional form doesn't work), and a linearised (Debye–Hückel) electrostatic double-layer force. Returns which particles adsorbed, when, and where.
- **`main_*.m`** — the driver. Sets the physical parameters for one system + one electrolyte concentration (channel geometry, zeta potentials, the array of particle radii to sweep, particle concentration φ) and calls the matching `fun_*.m` for each radius, then saves the raw trajectories to a `.mat` file.
- **`force_check_*.m` / `DLVO_*.m` / `DLVOE_*.m` / `PE_*.m`** — diagnostic scripts. These reload the precomputed van der Waals lookup table for a system+concentration and plot the van der Waals, electrostatic, and total DLVO force or energy profile vs. separation distance — used to sanity-check the interaction before trusting a `main_*` run, and to generate figures like Chapter 5's Fig. 5.1.

`hamaker-function/silica-silica/fun_SiO2_1000mM.m` is the odd one out — it's a kernel that uses a precomputed global Hamaker "constant" `A` directly rather than the piecewise-interpolated retarded profile, used for the retarded-vs-non-retarded comparison in Chapter 5.

`composite-particles/titania-silica/` also has two analysis-side helpers: `dep_kinetics_func.m` and `diff_DLVO.m`, used to build cross-system plots like deposition rate vs. Debye length (Chapter 6).

## Legacy baseline model

`legacy-baseline-model/` holds the pre-existing approach this thesis modified: `PS_call_fun.m` calculates van der Waals force from purely empirical piecewise power-law fits (no Lifshitz theory, no explicit retardation treatment), and `analysis_deposition.m` is the general post-processing pipeline (collector efficiency vs. time/particle size from raw `.mat` trajectory output) that both the old and new models share. `polyvalfitting.m` is an early curve-fitting experiment from that phase. Per the thesis preface, the original MATLAB script for calculating retarded van der Waals force (prior to my modification) and the image-analysis code were written by Dr Avinash Ashok; the per-material-system models in `composite-particles/`, `capsules/`, and `homogeneous-particles/` are the modification and extension of that starting point developed during this thesis.

## Running this code

Requires MATLAB (developed with the Curve Fitting Toolbox for the `.sfit` fits referenced in some diagnostic scripts, not included here). Each `main_*.m` is self-contained once its matching `fun_*.m` is on the path — set the `save(...)` path near the bottom of the driver to a local results folder before running. The van der Waals lookup tables (`NvdW*.xlsx`) that several kernels call via `interp1` are generated separately from the Hamaker function calculation and aren't all included in this trimmed repo; ask if you want the generation scripts too.

## Status

This is portfolio code from a completed PhD project, cleaned up and reorganised for readability — not an actively maintained package. Original per-concentration duplication (the real project has near-identical files for every electrolyte concentration tested) has been trimmed to representative examples per system.

## License

Not yet decided — TODO before making this public.
