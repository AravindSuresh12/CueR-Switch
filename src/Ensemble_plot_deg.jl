using DelimitedFiles
using Plots

#The ensembles chosen are with the dual criterion of best fit + best prediction of dose response


#PC_Ensemble=readdlm("./poets_ensemble_T_20/PC_T3.dat") 
PC_Ensemble=readdlm("./poets_ensemble_W/PC_T6.dat")
#PC_Ensemble=readdlm("./poets_ensemble_T_1000/PC_T10.dat")


mRNA_deg_CueR=log(2)./PC_Ensemble[11,:] #calculating half life from deg constants
mRNA_deg_Venus=log(2)./PC_Ensemble[12,:]
protein_deg_CueR=log(2)./PC_Ensemble[13,:]
protein_deg_Venus=log(2)./PC_Ensemble[14,:]
#protein_deg_sigma_70=log(2)./PC_Ensemble[19,:] #extra, use if needed
 
P1= Plots.histogram(mRNA_deg_CueR, bins=50, xlabel="CueR mRNA half life (min)", ylabel="Count", xguidefontsize=4, yguidefontsize=4,legend=:false, tickfontsize=4, grid=false)

P2= Plots.histogram(mRNA_deg_Venus, bins=50, xlabel="Venus mRNA half life (min)", ylabel="Count", xguidefontsize=4, yguidefontsize=4,legend=:false, tickfontsize=4, grid=false)

P3= Plots.histogram(protein_deg_CueR, bins=50, xlabel="CueR protein half life (days)", ylabel="Count", xguidefontsize=4, yguidefontsize=4,legend=:false, tickfontsize=4, grid=false)

P4= Plots.histogram(protein_deg_Venus, bins=50, xlabel="Venus degradation constant (days)", ylabel="Count", xguidefontsize=4, yguidefontsize=4,legend=:false, tickfontsize=4, grid=false)


plot(P1,P2,P3,P4, layout=(2,2))


Plots.savefig("./deg_constants.pdf")
#Plots.savefig("./Images/Pareto_1000/Degradation_constants.pdf")