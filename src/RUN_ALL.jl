using Shell

# # DELETE copper DYNAMICS SIMULATION FILES - Uncomment me for deleting previous files 
# Shell.run("rm -f simulated/copper_dynamics/0uM/*");
# Shell.run("rm -f simulated/copper_dynamics/5uM/*");
# Shell.run("rm -f simulated/copper_dynamics/10uM/*");
# Shell.run("rm -f simulated/copper_dynamics/50uM/*");
# Shell.run("rm -f simulated/copper_dynamics/100uM/*");
# Shell.run("rm -f simulated/Sensitivity_matrix_copper.dat") #refers to the sensitivity array from Morris
# Shell.run("rm -f simulated/copper_dynamics/Sensitivity_matrix_copper.dat") #refers to the sensitivity array for dose response ensemble


# # DELETE POETS ENSEMBLE FILES
# Shell.run("rm -f simulated/POETS/*");
# Shell.run("rm -f simulated/POETS_N_Vary/N=5/*");
# Shell.run("rm -f simulated/POETS_N_Vary/N=100/*");

#script takes about 45 min to run appropos
# ESTIMATE PARAMETERS
include("Parameter_Estimation_W_splined.jl");
#include("Parameter_Estimation_W_splined_N=100.jl"); UNCOMMENT ME FOR RUNNING THIS ENSEMBLE FOR MODEL FIT
#include("Parameter_Estimation_W_splined_N=5.jl"); UNCOMMENT ME FOR RUNNING THIS ENSEMBLE FOR MODEL FIT
include("Updated_Driver.jl"); #uncomment specific graphs in the file
include("Parameter_set.jl");

# SENSITIVITY
include("sensitivity-copper-cuer-updated.jl");
include("visualize-sensitivity-array.jl");


# AFTER ESTIMATING PARAMETERS
include("copper_dynamics_automated.jl") #MAKE SURE to confirm the simulation file being used, eg. PC_T5.dat
include("dose_response_ensemble.jl") #MAKE SURE to confirm the simulation file being used, eg. PC_T5.dat

#dose response plot-  #MAKE SURE to confirm the simulation file being used, eg. PC_T5.dat
include("plot_dr_venus.jl") 

# Species trajectories-  #MAKE SURE to confirm the simulation file being used, eg. PC_T5.dat
include("plot_CueR_mrna.jl")
include("plot_venus_mrna.jl")
include("plot_prot_CueR.jl")
include("plot_prot_Venus.jl")

#parameter distribution ensembles-  #MAKE SURE to confirm the simulation file being used, eg. PC_T5.dat

include("Ensemble_plot_W.jl")
include("Ensemble_plot_tc.jl")
include("Ensemble_plot_Kc.jl")
include("Ensemble_plot_Kc.jl")
include("Ensemble_plot_n.jl")
include("Ensemble_misc.jl")

#include("Model_pred_N_vary.jl) #for comparisons based on the N files

include("Model_stats.jl") #calculates the overall model's predictive power as stat terms



