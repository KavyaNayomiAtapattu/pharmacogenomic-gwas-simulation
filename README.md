# Pharmacogenomic GWAS Simulation

R simulation code for the master's thesis:
**"Addressing Statistical Complexities in Modeling Treatment Effects 
from Repeated Measurements: A Simulation Study in Pharmacogenomics"**

University of Helsinki, Faculty of Science — May 2026

## Requirements

R packages required (install before running):

```r
install.packages(c("dplyr", "tidyr", "ggplot2", 
                   "patchwork", "gridExtra", "grid"))

# SlopeHunter must be installed from GitHub:
install.packages("remotes")
remotes::install_github("Osmahmoud/SlopeHunter")
```

## How to run

Open `SimulationMay12v1.R` in RStudio and run the entire script.
All figures and tables will be written to `simulation_outputs_May12/`.
Runtime: approximately 2–4 hours for 1000 iterations on a standard laptop.

## Outputs

- Fig01–Fig22: all figures appearing in the thesis
- Table_Bias_M1_M2.csv: bias table for Research Objective 2
- Table_Bias_Obj3.csv: bias table for Research Objective 3
