include("Include.jl")

PC_Ensemble=readdlm("./simulated/POETS/PC_T3.dat") #best fit
#PC_Ensemble=readdlm("./poets_ensemble_T_100/PC_T3.dat")
#PC_Ensemble=readdlm("./poets_ensemble_T_1000/PC_T10.dat")



W_T7RNAP=PC_Ensemble[1,:]
W_apo=PC_Ensemble[2,:]
W_holo=PC_Ensemble[3,:]

P1= PyPlot.hist(W_T7RNAP, bins=50)
P1=PyPlot.xlabel("Pseudo Energy (kbT)", fontsize=20)
P1=PyPlot.ylabel("count", fontsize=20)


# P2=hist(W_apo, bins=50, xlabel="Pseudo energy distribution W_apo_CueR_DNA (kT)", ylabel="Count", grid=false)

# P3=hist(W_holo, bins=50, xlabel="Pseudo energy distribution W_holo_CueR  (kT)", ylabel="Count", grid=false)


# PLOT1=plot(P1,P2,P3)

Plots.savefig("./plots/Energy_distribution.pdf")
#Plots.savefig("./Images/Pareto_100/pVec_no_perturb/Energy_distribution-100_runs.pdf")
#Plots.savefig("./Images/Pareto_1000/pVec_no_perturb/Energy_distribution-1000_runs.pdf")