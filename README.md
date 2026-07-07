[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21238413.svg)](https://doi.org/10.5281/zenodo.21238413)


# spatsoc
Behavioural plasticity at the spatial-social interface.

## Project description

This repository provides code and processed data supporting the manuscript entitled: **Behavioural plasticity at the spatial-social interface: Predation risk modulates density-dependent breeding dispersion**.

The study was supported by the **National Science Centre (NCN) in Poland (grant no. 2018/29/B/NZ8/00066)**. The computational resources supporting this work were provided by the **Poznań Supercomputing and Networking Centre (PCSS) (grant no. pl0090-01)**.

## Repository structure
- R
  - 01_index.R - demonstration of Normalised Spatial Dispersion Index (NSDI) properties, simulation of the dummy data and NSDI calculation
  - 02_GAMM.R - fitting GAMM models to the empirical, processed data
  - 03_visualisation.R - plotting figures for publication
  - utils.R - utility functions for plotting and NSDI calculation

- data
  - data.RDS - main dataset containing processed data used in GAMM modelling (R binary format)
  - data.csv - equivalent version of the dataset in plain text format

- Figures
  - Figure_1 - Frequency distribution of the NSDI across all surveyed northern lapwing breeding plots.
  - Figure_2 - Marginal main effects of socio-ecological predictors on the spatial breeding dispersion of lapwings, quantified by the NSDI.
  - Figure_3 - Predicted interaction effects between local conspecific density and predation risk on the spatial breeding dispersion (quantified by the NSDI) of lapwings.

## Data

The raw bird counts are provided by BirdLife Poland and the Chief Inspectorate for Environmental Protection (GIOS). The GIOS data policy does not permit the release of raw data. For this reason, we are only making the processed data available (in the 'data/data.RDS' and 'data/data.csv' files), which includes densities estimated using the distance sampling method and the NSDI calculated. The original observers' and plots' IDs were pseudonymised. However, to enable the replication of the entire data processing chain, we provide sample data, which is an exact replica of the original data. We then use this data to demonstrate how to calculate the NSDI. 

- The processed dataset includes:
  - hooded crow (*Corvus cornix*) and northern lapwing (*Vanellus vanellus*) densities [pairs km⁻²],
  - NSDI calculated across all surveyed northern lapwing breeding plots,
  - Land cover and climatic variables obtained from MODIS (MCD12Q1 Version 6.1) instrument and TerraClimate database.

## Requirements

- **R 4.5.2**
- **Required R packages:** `mgcv`, `pbapply`, `dplyr`

## Workflow & Usage

1. Open as R project in RStudio or set the working directory to the repository root on your machine.
2. Run the scripts sequentially:
   - 1. 01_index.R
   - 2. 02_GAMM.R 
   - 3. 03_visualisation.R
