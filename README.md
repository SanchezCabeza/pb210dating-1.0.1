# pb210dating

R package for classical Pb-210 sediment geochronology, uncertainty propagation and age-model validation.

## Overview

pb210dating provides a reproducible workflow for dating sediment cores from Pb-210 activity profiles measured by alpha or gamma spectrometry.

The package implements the Constant Flux (CF) and Constant Flux Constant Sedimentation (CFCS) dating models, following the unified formulation and nomenclature of Sanchez-Cabeza and Ruiz-Fernandez (2012).

It also provides functions for:

- data input and preparation;
- radioactive decay correction;
- supported Ra-226 estimation;
- equilibrium assessment;
- missing Pb-210 inventory estimation;
- calculation of sediment accumulation rate (SAR) and mass accumulation rate (MAR);
- Monte Carlo propagation of dating uncertainties;
- age-model validation;
- visualization of radionuclide profiles and age models;
- standardized output of dating results.

## Version

This repository contains pb210dating version 1.0.1, corresponding to the software version used in the associated manuscript.

## Requirements

- R >= 4.3
- lubridate

## Installation

The released package is available from CRAN and can be installed with:

    install.packages("pb210dating")

Then load the package:

    library(pb210dating)

## Quick test

A quick test is provided in:

    examples/quick_test.R

After installing the package, run the test from the root directory of this repository with:

    source("examples/quick_test.R")

The script runs the complete Pb-210 dating workflow using the package-provided TEHUA-II example dataset.

The quick test uses 10 Monte Carlo simulations for both the CF and CFCS models so that it can be executed rapidly. For scientific applications, substantially larger numbers of Monte Carlo simulations should be used.

## Example data

Two example datasets are included in the package:

- TehuaII.csv — example data for a complete sediment core.
- TehuaIIMI.csv — example data for a sediment core with incomplete Pb-210 inventory.

The files are located in the package inst/extdata directory and can be accessed from R with:

    system.file("extdata", "TehuaII.csv", package = "pb210dating")

and:

    system.file("extdata", "TehuaIIMI.csv", package = "pb210dating")

## Main functions

The package provides functions for the main stages of the Pb-210 dating workflow, including:

- ReadData() — reads sediment-core data.
- Prepare4Dating() — prepares data before dating.
- DecayCorrection() — applies radioactive decay corrections.
- CompleteProfile() — completes the Pb-210 profile when required.
- ConstantRa() — estimates supported Ra-226 activity.
- Equilibrium() — determines the equilibrium layer.
- EquilibriumSet() — sets the dating interval at a specified layer.
- Pb210CF() — calculates the Constant Flux model.
- Pb210CFCS() — calculates the Constant Flux Constant Sedimentation model.
- MissingInventory() — estimates missing Pb-210 inventory.
- AgeModel() — calculates and plots age models.
- PlotProfiles() — plots radionuclide profiles.
- PlotMARSAR() — plots mass and sediment accumulation rates.
- Output() — exports dating results.
- Pb210Dating() — runs the complete Pb-210 dating workflow.

Function documentation can be accessed from R, for example:

    ?Pb210CF
    ?Pb210CFCS
    ?AgeModel
    ?Pb210Dating

## Documentation

Documentation for pb210dating is available from CRAN:

https://CRAN.R-project.org/package=pb210dating

## Reproducibility

This repository contains the complete source code of pb210dating version 1.0.1, the example datasets included in the package, and a quick-test script.

The repository provides the source code and example material required to inspect and run the software version used in the associated manuscript.

## License

pb210dating is distributed under the GNU General Public License (GPL >= 2).

## References

Sanchez-Cabeza, J.A. and Ruiz-Fernandez, A.C. (2012). 210Pb sediment radiochronology: an integrated formulation and classification of dating models. Geochimica et Cosmochimica Acta, 82, 183-200. https://doi.org/10.1016/j.gca.2010.12.024

Sanchez-Cabeza, J.A., Ruiz-Fernandez, A.C., Ontiveros-Cuadras, J.F., Bernal, L.H.P. and Olid, C. (2014). Monte Carlo uncertainty calculation of 210Pb chronologies and accumulation rates of sediments and peat bogs. Quaternary Geochronology, 23, 80-93. https://doi.org/10.1016/j.quageo.2014.06.002
