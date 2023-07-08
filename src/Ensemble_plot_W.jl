include("Include.jl")

PC_Ensemble=readdlm("./simulated/POETS/PC_T5.dat") #best fit
#PC_Ensemble=readdlm("./poets_ensemble_T_100/PC_T3.dat")
#PC_Ensemble=readdlm("./poets_ensemble_T_1000/PC_T10.dat")



W_T7RNAP=PC_Ensemble[1,:]
W_apo=PC_Ensemble[2,:]
W_holo=PC_Ensemble[3,:]


fig=figure() 

subplot(1,3,1)
PyPlot.plot()
PyPlot.hist(W_T7RNAP, bins=50)
PyPlot.xlabel("Pseudo Energy - dE_T7RNAP (kbT)", fontsize=5)
PyPlot.ylabel("count", fontsize=5)

subplot(1,3,2)
PyPlot.plot()
PyPlot.hist(W_apo, bins=50)
PyPlot.xlabel("Pseudo Energy dE_Apo-CueR (kbT)", fontsize=5)
PyPlot.ylabel("count", fontsize=5)

subplot(1,3,3)
PyPlot.plot()
PyPlot.hist(W_holo, bins=50)
PyPlot.xlabel("Pseudo Energy dE_holo-CueR (kbT)", fontsize=5)
PyPlot.ylabel("count", fontsize=5)


PyPlot.savefig("./plots/Energy_distribution.pdf")
#Plots.savefig("./Images/Pareto_100/pVec_no_perturb/Energy_distribution-100_runs.pdf")
#Plots.savefig("./Images/Pareto_1000/pVec_no_perturb/Energy_distribution-1000_runs.pdf")