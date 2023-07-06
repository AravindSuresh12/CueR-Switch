using Shell

# DELETE GLUCONATE DYNAMICS SIMULATION FILES 
Shell.run("rm -f simulated/copper_dynamics/0uM/*");
Shell.run("rm -f simulated/copper_dynamics/5uM/*");
Shell.run("rm -f simulated/copper_dynamics/10uM/*");
Shell.run("rm -f simulated/copper_dynamics/50uM/*");
Shell.run("rm -f simulated/copper_dynamics/100uM/*");


# DELETE POETS ENSEMBLE FILES
Shell.run("rm -f simulated/POETS/*");

# ESTIMATE PARAMETERS
include("Parameter_Estimation_W_splined.jl");
include("Updated_Driver.jl");


# AFTER ESTIMATING PARAMETERS
include("copper_dynamics_automated.jl") #MAKE SURE to confirm the simulation file being used, eg. PC_T5.dat
include("dose_response_ensemble.jl") #MAKE SURE to confirm the simulation file being used, eg. PC_T10.dat
include("plot_dose_response_venus.jl")
include("plot_dose_response_bfp.jl")
include("plot_venus_protein.jl")
include("plot_venus_mrna.jl")
include("plot_bfp_protein.jl")
include("plot_bfp_mrna.jl")

# SENSITIVITY
include("sensitivity-gntr-gluconate-updated.jl")
include("plot_sensitivity_array.jl")
