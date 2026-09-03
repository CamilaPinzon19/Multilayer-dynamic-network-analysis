# Multilayer-dynamic-network-analysis

## Description

Analysis of the structure and temporal evolution of a multilayer dynamic network of international interactions using ICEWS data. The analysis considers four interaction layers: Material-, Material+, Verbal-, and Verbal+.

The study includes network centrality analysis, structural network features, Matrix Autoregressive (MAR) modeling, impulse response analysis, residual diagnostics, and simulation experiments.

## Project Structure

- `analysis.R`: Data analysis, MAR model estimation, diagnostics, impulse response analysis, and simulation experiments
- `visualizations.Rmd`: Static and dynamic network visualizations
- `icews_tensor_proyecto.RData`: Processed ICEWS data
- `README.md`: Project documentation

## Methodology

1. Selection of the 25 most central countries using PageRank centrality
2. Computation of structural network features across layers and time
3. Trend and stationarity analysis joint with standardization
4. Estimation of a Matrix Autoregressive (MAR) model
5. Impulse response and residual diagnostics
6. Monte Carlo simulations to evaluate model performance

## Data Source

Integrated Crisis Early Warning System (ICEWS) event data, used to construct multilayer dynamic networks of international interactions.
