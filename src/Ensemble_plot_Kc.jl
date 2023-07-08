include("Include.jl")

PC_Ensemble=readdlm("./simulated/POETS/PC_T6.dat") #best fit

K_CueR_OCueR_apo=PC_Ensemble[5,:]
K_CueR_OCueR_holo=PC_Ensemble[6,:]
half_max_induction_CuSO4=PC_Ensemble[18,:]
fig=figure() 


subplot(1,3,1)
PyPlot.hist(K_CueR_OCueR_apo, bins=50)
PyPlot.xlabel("K_CueR_OCueR_apo (μM)", fontsize=5)
PyPlot.ylabel("count", fontsize=5)

subplot(1,3,2)
PyPlot.hist(K_CueR_OCueR_holo, bins=50)
PyPlot.xlabel("K_CueR_OCueR_holo (μM)", fontsize=5)
PyPlot.ylabel("count", fontsize=5)

subplot(1,3,3)
PyPlot.hist(half_max_induction_CuSO4 , bins=50)
PyPlot.xlabel("half_max_induction_CuSO4 (μM)", fontsize=5)
PyPlot.ylabel("count", fontsize=5)

PyPlot.savefig("./plots/Dissociation_constants.pdf")



