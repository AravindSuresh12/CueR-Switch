using Shell

# # DELETE copper DYNAMICS SIMULATION FILES 
# Shell.run("rm -f simulated/copper_dynamics/0uM/*");
# Shell.run("rm -f simulated/copper_dynamics/5uM/*");
# Shell.run("rm -f simulated/copper_dynamics/10uM/*");
# Shell.run("rm -f simulated/copper_dynamics/50uM/*");
# Shell.run("rm -f simulated/copper_dynamics/100uM/*");
# Shell.run("rm -f simulated/Sensitivity_matrix_copper.dat") #refers to the sensitivity array from Morris
# Shell.run("rm -f simulated/copper_dynamics/Sensitivity_matrix_copper.dat") #refers to the sensitivity array for dose response ensemble


# # DELETE POETS ENSEMBLE FILES
# Shell.run("rm -f simulated/POETS/*");

#script takes about 45 min to run appropos
# ESTIMATE PARAMETERS
include("Parameter_Estimation_W_splined.jl");
include("Updated_Driver.jl");
include("Parameter_set.jl");

# SENSITIVITY
include("sensitivity-copper-cuer-updated.jl");
include("visualize-sensitivity-array.jl");


# AFTER ESTIMATING PARAMETERS
include("copper_dynamics_automated.jl") #MAKE SURE to confirm the simulation file being used, eg. PC_T5.dat
include("dose_response_ensemble.jl") #MAKE SURE to confirm the simulation file being used, eg. PC_T5.dat

#dose response plot
include("plot_dr_venus.jl") 

# Species trajectories
include("plot_CueR_mrna.jl")
include("plot_venus_mrna.jl")
include("plot_prot_CueR.jl")
include("plot_prot_Venus.jl")

#parameter file stored in simulated

include("Ensemble_plot_W.jl");


#Parameter Ensemble distributions for N=50


