using DelimitedFiles
using Plots

PC_Ensemble=readdlm("./poets_ensemble_W/PC_T6.dat") #best fit



n_CueR_OCueR=PC_Ensemble[4,:]
K_CueR_OCueR_apo=PC_Ensemble[5,:]
K_CueR_OCueR_holo=PC_Ensemble[6,:]
K_diss_CuSO4=PC_Ensemble[18,:]

P15= Plots.histogram(n_CueR_OCueR, bins=50, xlabel="n_CueR_OCueR", ylabel="Count", xguidefontsize=4, yguidefontsize=4,legend=:false, tickfontsize=4, grid=false)

P16= Plots.histogram(K_CueR_OCueR_apo, bins=50, xlabel="K_CueR_OCueR_apo (uM)", ylabel="Count", xguidefontsize=4, yguidefontsize=4,legend=:false, tickfontsize=4, grid=false)

P17= Plots.histogram(K_CueR_OCueR_holo, bins=50, xlabel="K_CueR_OCueR_holo (uM)", ylabel="Count", xguidefontsize=4, yguidefontsize=4,legend=:false, tickfontsize=4, grid=false)


P18= Plots.histogram(K_diss_CuSO4, bins=50, xlabel="Half max induction CuS04 (uM)", ylabel="Count", xguidefontsize=4, yguidefontsize=4,legend=:false, tickfontsize=4, grid=false)

plot(P15,P16,P17,P18, layout=(2,2))
Plots.savefig("K_values.pdf")
#Plots.savefig("./Images/Pareto_100/Kc.pdf")
#Plots.savefig("./Images/Pareto_1000/Kc.pdf")


