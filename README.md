# Optimal Reaction Fluxes and Ensemble Fluxes in various Catalyst Arrangements

(c) 2021 Florian Hinzpeter

This repository comprises the software for simulations that were used to generate the results presented in the manuscript "Trade-offs and design principles in the spatial organization of catalytic particles" by Hinzpeter F., Tostevin F., Buchner A. & Gerland U. (Nat. Phys. 2021).

## Description

The code is structured into a folder with the simulations of the 2 dimensional system (2_Dim) and a folder with the simulations of the 3 dimensional system (3_Dim).
Each folder contains the following subfolders with the corresponding functionalities
- **arrangement_optimization:**
  This folder contains the code that computes the flux optimizing arrangements
  for a single first catalyst at the system center and several second catalysts around it.
  The optimization procedure is a random search procedure.
  To facilitate the search process mirror and spherical symmetries were imposed.
- **sample_cat_arrangements**:
  This folder contains the code that computes random samples of delocalized, clustered
  and pair arrangements. For each of these classes arrangement ensembles are generated.
- **fluxes_comsol_livelink**:
  This folder contains the code that computes individual fluxes given a catalysts
  arrangement and all fluxes for an catalyst arrangement ensemble.
  The flux computation is done using Comsol together with the Matlab interface LiveLink.
- **statistics_and_visualizations**:
  This folder contains functionality for visualizations and computation of statistics on the
  flux ensembles for the different arrangement classes.

## Dependencies

* Matlab 9.0
* Comsol 5.2 with the LiveLink interface for Matlab.

## License

This project is licensed under the MIT License - see the LICENSE.md file for details
