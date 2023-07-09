# CueR-Switch

# Biophysical model of a cell-free copper-responsive biosensor

### Introduction
This repository contains code written in the [Julia](https://www.julialang.org) programming language for the cell-free copper responsive biosensor: 

Paper Title: "Quantitative model for Copper biosensor using CueR transcription factor in cell-free systems"


### Installation and Requirements

To get the model codes, you can download this model repository as a zip file, clone or pull it using the command:

	git pull https://github.com/AravindSuresh12/CueR-Switch
or

	git clone https://github.com/AravindSuresh12/CueR-Switch

The ``src`` directory contains the code for the model, the ``data`` directory contains the experimental data, and the ``simulated`` directory contains the simulated model files. plots- contains the plots generated. ``misc`` directory contains graphs and processed data for controls. Each folder has a README.md to better understand what files are there in it. 

### Scripts- Essential
Script | Description
---: | ---
RUN_ALL.jl | Runs the entire model simulation and gives outputs. 
``parameter_estimation_W_splined.jl`` | Solves the model equations for the ensemble of parameters sets for the test case, 10mM gluconate. Saves solutions in the ``simulated/POETS`` directory


