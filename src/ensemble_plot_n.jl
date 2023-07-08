include("Include.jl")
PC_Ensemble=readdlm("./simulated/POETS/PC_T6.dat") #best fit

n_CueR_OCueR=PC_Ensemble[4,:]
n_half_max_induction_CuSO4=PC_Ensemble[19,:]

fig=figure() 

subplot(1,2,1)
PyPlot.hist(n_CueR_OCueR, bins=50)
PyPlot.xlabel("n_CueR_OCueR", fontsize=5)
PyPlot.ylabel("count", fontsize=5)

subplot(1,2,2)
PyPlot.hist(n_half_max_induction_CuSO4 , bins=50)
PyPlot.xlabel("half_max_induction_CuSO4", fontsize=5)
PyPlot.ylabel("count", fontsize=5)

PyPlot.savefig("./plots/n_values.pdf")
