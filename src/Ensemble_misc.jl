include("Include.jl")

PC_Ensemble=readdlm("./simulated/POETS/PC_T5.dat") #best fit for N=20

saturation_constant_TL=PC_Ensemble[16,:]
saturation_constant_TX=PC_Ensemble[17,:]
half_life_TL=PC_Ensemble[15,:]
fig=figure() 


subplot(1,3,1)
PyPlot.hist(saturation_constant_TL, bins=50)
PyPlot.xlabel("saturation_constant_TL(KL)(s)", fontsize=5)
PyPlot.ylabel("count", fontsize=5)


subplot(1,3,2)
PyPlot.hist(saturation_constant_TX, bins=50)
PyPlot.xlabel("saturation_constant_TX(KX)(s)", fontsize=5)
PyPlot.ylabel("count", fontsize=5)

subplot(1,3,3)
PyPlot.hist(half_life_TL, bins=50)
PyPlot.xlabel("half_life translation(h)", fontsize=5)
PyPlot.ylabel("count", fontsize=5)


PyPlot.savefig("./plots/Kx,l_and_half_life.pdf")



