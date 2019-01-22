# water-resource-emulation
Supplementary materials for 
Owen N.E. &amp; Liuzzo L. "Impact of land use on water resources via a Gaussian process emulator with dimension reduction"

Order for running code:

1. Create experimental design in R using create_clhs_design.R OR use existing design saved in text file clhs_design.txt

2. Run Rushton Model in MATLAB using run_rushton_model.m code. This code requires the following files:
- experimental design in clhs_design.txt
- rushton model code in rushton_model.m
- basin data in basin_data.txt
- crop data in crop_data.txt
- parameters in parameters.txt
- weather data in weather_data.txt
OR use existing model runs saved in text file model_runs.txt

3. Build Gaussian process emulator and perform leave-one-out cross validation in R using build_emulator.R. 
This code requires the following files:
- model runs in model_runs.txt
- R code for leave-one-out cross validation for Gaussian process emulator in loo_km.R 
- R code for transformation and propagation of uncertainty of Gaussian process emulator in transform_emulator.R

Additional code:

1. Create plot of Rushton model output (Figure 3) in R using rushton_model_plot.R. This code requires data in plot_model_runs.txt.
(plot_model_runs.txt can be created by running MATLAB code run_rushton_model_for_plot.m using data plot_design.txt)

2. Record the execution time of the Rushton Model experiment in MATLAB using time_rushton_model.m.
