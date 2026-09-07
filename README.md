# Multilayer-dynamic-network-analysis

## Description

Analysis of the structure and temporal evolution of a multilayer dynamic network of international interactions using ICEWS data. The analysis considers four interaction layers: Material-, Material+, Verbal-, and Verbal+.

The study includes network centrality analysis, structural network features, Matrix Autoregressive (MAR) modeling, impulse response analysis, residual diagnostics, and simulation experiments.

## Project Structure

- `analysis.R`: Data analysis, MAR model estimation, diagnostics, impulse response analysis, and simulation experiments
- `visualizations.Rmd`: Static and dynamic network visualizations
- `icews_tensor_processed.RData`: Processed ICEWS data
- `README.md`: Project documentation

## Methodology

1. Selection of the 25 most central countries using PageRank centrality
2. Computation of structural network features across layers and time
3. Trend and stationarity analysis joint with standardization
4. Estimation of a Matrix Autoregressive (MAR) model
5. Impulse response and residual diagnostics
6. Monte Carlo simulations to evaluate model performance

## Data

The analysis uses ICEWS event data described in Hoff (2015). The data consist of weekly interactions between countries classified into four types: negative material, positive material, negative verbal, and positive verbal interactions.

The processed data are stored as a tensor containing the interactions by country, layer, and time period.

Reference:

Hoff, P. D. (2015). *Multilinear tensor regression for longitudinal relational data*. Annals of Applied Statistics, 9(3), 1169–1193. https://doi.org/10.1214/15-AOAS839

## Reproducibility

- **R:** 4.4.2
- **RStudio:** 2023.09.1+494
- **Operating system:** macOS Big Sur 11.7.11
- **Processor:** Dual-Core Intel Core i5, 2.6 GHz
- **Memory:** 8 GB RAM
- **Architecture:** x86_64

### R Packages

- `dplyr` 1.1.4
- `expm` 1.0-0
- `fable` 0.5.0
- `feasts` 0.5.0
- `ggplot2` 4.0.3
- `ggraph` 2.2.2
- `igraph` 2.2.1
- `ndtv` 0.13.4
- `network` 1.20.0
- `networkDynamic` 0.12.0
- `psych` 2.6.5
- `scales` 1.4.0
- `tensorTS` 1.0.2
- `tidygraph` 1.3.1
- `timetk` 2.9.1
- `tsibble` 1.1.6
- `zoo` 1.8-14




