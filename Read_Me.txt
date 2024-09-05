

Doubly-robust machine learning based estimation methods for instrumental variables with an application to surgical care for cholecystitis
Kenta Takatsu, Alexander W. Levis, Edward Kennedy, Rachel Kelz, Luke Keele
Aug 21, 2024

This document describes the replication files for the above paper. The archive contains all of the programs used in the analysis. 

Application

For this application data are publicly available, but require obtaining a data use agreement from the states of PA, FL, and NY. Data must also be purchased from each state. As such, we do not include the data, but include all the R scripts used in the analysis.

1. Functions.R -- R scripts used in various analyses.
2. Nuisance-Est --- Estimates the nuisance functions.
3. Balance-Table.R -- Covariate balance results contained in the appendix.
4. Bootstrap.R -- Does the bootstrap resampling.
5. Main-Analysis.R -- Produces main results in the paper.

Simulation

1. sims-late -- Simulations targeting the LATE parameter
	In this directory there are three scripts for each DGP scenario. The summary file generates the figures in the paper.
2. sims-clate -- Simulations targeting the Conditional LATE parameter
	In this directory there are three scripts for each DGP scenario. The summary file generates the figures in the paper.
	There is one additional script that varies the instrument strength.

 