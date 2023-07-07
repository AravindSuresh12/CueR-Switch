using DelimitedFiles
using Plots

PC_Ensemble=readdlm("./poets_ensemble_T_20/PC_T3.dat") #best fit
#PC_Ensemble=readdlm("./poets_ensemble_T_100/PC_T3.dat")
#PC_Ensemble=readdlm("./poets_ensemble_T_1000/PC_T10.dat")



mRNA_tc_GntR=PC_Ensemble[11,:]
mRNA_tc_Venus=PC_Ensemble[12,:]
protein_tc_GntR=PC_Ensemble[13,:]
protein_tc_Venus=PC_Ensemble[14,:]


P6= Plots.histogram(mRNA_tc_GntR, bins=50, xlabel="GntR mRNA time constant (h)", ylabel="Count", xguidefontsize=4, yguidefontsize=4,legend=:false, tickfontsize=4, grid=false)

P7= Plots.histogram(mRNA_tc_Venus, bins=50, xlabel="Venus mRNA time constant (h)", ylabel="Count", xguidefontsize=4, yguidefontsize=4,legend=:false, tickfontsize=4, grid=false)

P8= Plots.histogram(protein_tc_GntR, bins=50, xlabel="GntR protein time constant (h)", ylabel="Count", xguidefontsize=4, yguidefontsize=4,legend=:false, tickfontsize=4, grid=false)

P9= Plots.histogram(protein_tc_Venus, bins=50, xlabel="Venus mRNA time constant (h)", ylabel="Count", xguidefontsize=4, yguidefontsize=4,legend=:false, tickfontsize=4, grid=false)

plot(P6,P7,P8,P9, layout=(2,2))


Plots.savefig("./Images/Pareto_20/pVec_no_perturb/tc-20_runs.pdf")
#Plots.savefig("./Images/Pareto_100/pVec_no_perturb/tc-100_runs.pdf")
#Plots.savefig("./Images/Pareto_1000/pVec_no_perturb/tc-1000_runs.pdf")