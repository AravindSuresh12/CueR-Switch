include("Include.jl")

venus_data = readdlm("./simulated/copper_dynamics/Sensitivity_matrix_copper.dat")

# first column contains the copper concentration
copper_concentration = venus_data[:,1]
# the remaining columns contain the data
venus_response = venus_data[:,2:end]


μ_venus = mean(venus_response,dims=2)
σ_venus = std(venus_response,dims=2)


# read dose response data
dose_response_exp = CSV.read("./data/dose_response.csv",DataFrame)
copper_concentration_exp = dose_response_exp[!,"Copper_conc(uM)"] # mM

venus_mean_exp = dose_response_exp[!,"Venus(uM)"] # μM
venus_stderr_exp = dose_response_exp[!,"STERR(uM)"] # μM

LB_venus = μ_venus .- (1.96/sqrt(1))*σ_venus
UB_venus = μ_venus .+ (1.96/sqrt(1))*σ_venus


PyPlot.clf()
PyPlot.figure(1)
fill_between(copper_concentration, vec(UB_venus), vec(LB_venus), color="#5C995B", alpha=0.45)

PyPlot.plot(copper_concentration,μ_venus,"-",color="black",lw=1.5)
yerr_array_venus = transpose([venus_stderr_exp venus_stderr_exp])


PyPlot.errorbar(copper_concentration_exp, venus_mean_exp,yerr=yerr_array_venus,fmt="o",mfc="#252525",mec="#252525",color="#252525", lw=1.5)




# labels - 
PyPlot.plot(figsize=(4,4))
PyPlot.xlabel("[Copper](μM)", fontsize=16)
PyPlot.ylabel("[Venus] (μM)", fontsize=16)
# PyPlot.axis([-0.5,16.5,-50,1600])
PyPlot.xticks(fontsize=14)
PyPlot.yticks(fontsize=14)
PyPlot.tight_layout()
PyPlot.savefig("./plots/Dose_response.pdf")