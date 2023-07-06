include("Include.jl")

venus_data = readdlm("./simulated/dose_response_simulations/results_matrix_venus.dat")
bfp_data = readdlm("./simulated/dose_response_simulations/results_matrix_bfp.dat")

# first column contains the gluconate concentration
gluconate_concentration = venus_data[:,1]
# the remaining columns contain the data
venus_response = venus_data[:,2:end]
bfp_response = bfp_data[:,2:end]

μ_venus = mean(venus_response,dims=2)
σ_venus = std(venus_response,dims=2)
μ_bfp = mean(bfp_response,dims=2)
σ_bfp = std(bfp_response,dims=2)


# LB = μ .- (1.96/sqrt(1))*σ
# UB = μ .+ (1.96/sqrt(1))*σ

# read dose response data
dose_response_exp = CSV.read("./data/FINAL_FULL_DOSE_RESPONSE_POOLED.csv",DataFrame)
gluconate_concentration_exp = dose_response_exp[!,"logGlucose_mM"] # mM
gluconate_concentration_exp = 10.0 .^ gluconate_concentration_exp # mM

venus_mean_exp = dose_response_exp[!,"Venus_µM"] # μM
venus_stderr_exp = dose_response_exp[!,"STDERR_V_µM"] # μM
bfp_mean_exp = dose_response_exp[!,"BFP_µM"] # μM
bfp_stderr_exp = dose_response_exp[!,"STDERR_BFP_µM"] # μM


# # Plot simulated protein
# Plots.plot()
# # p1 = Plots.plot!(log.(gluconate_concentration), μ, ribbon=(LB,UB),label = "Dose Response Simulated", legend = :topright,xlabel="Gluconate Concentration (mM)",ylabel = "Venus Concentration (μM)", lw=3,fillalpha=0.35,c=:orange)
# # p1 = Plots.plot!(log.(gluconate_concentration), μ, fillbetween=[LB UB],label = "Dose Response Simulated", legend = :topleft,xlabel="Gluconate Concentration (mM)",ylabel = "Venus Concentration (μM)", lw=3,fillalpha=0.35,c=:orange)
# p1 = Plots.plot!(log.(gluconate_concentration), μ, ribbon=2.576*(σ),label = "Dose Response Simulated", legend = :topright,xlabel="log[ Gluconate Concentration (mM) ]",ylabel = "Venus Concentration (μM)", lw=3,fillalpha=0.35,grid=false)


# # Plot experimental protein
# p2 = Plots.scatter!(log.(gluconate_concentration_exp),venus_mean_exp,label = "Dose Response Exp", markercolor = "black", legend = :topleft, yerror=venus_stderr_exp)
# Plots.plot(p2)
# # Plots.savefig("simulated/dose_response_simulations_5/Dose_response_ABHI_smooth_normalized.pdf")
# Plots.savefig("simulated/dose_response_simulations_5/Dose_response_ABHI_smooth.pdf")

LB_venus = μ_venus .- (1.96/sqrt(1))*σ_venus
UB_venus = μ_venus .+ (1.96/sqrt(1))*σ_venus
LB_bfp = μ_bfp .- (1.96/sqrt(1))*σ_bfp
UB_bfp = μ_bfp .+ (1.96/sqrt(1))*σ_bfp

PyPlot.clf()
PyPlot.figure(1)
fill_between(log10.(gluconate_concentration), vec(UB_venus), vec(LB_venus), color="#5C995B", alpha=0.45)
# PyPlot.fill_between(log10.(gluconate_concentration), vec(UB_bfp), vec(LB_bfp), color="#CDE1BF", alpha=0.60)

# for i in 1:length(venus_response[1,:])
    # PyPlot.plot(log10.(gluconate_concentration),venus_response[:,i],"-",color="#5C995B",lw=0.5,alpha=0.10)
# end

PyPlot.plot(log10.(gluconate_concentration),μ_venus,"-",color="black",lw=1.5)
# PyPlot.plot(log10.(gluconate_concentration),μ_bfp,"-",color="black",lw=1.5)
yerr_array_venus = transpose([venus_stderr_exp venus_stderr_exp])
# yerr_array_bfp = transpose([bfp_stderr_exp bfp_stderr_exp])

PyPlot.errorbar(log10.(gluconate_concentration_exp), venus_mean_exp,yerr=yerr_array_venus,fmt="o",mfc="#252525",mec="#252525",color="#252525", lw=1.5)
# PyPlot.errorbar(log10.(gluconate_concentration_exp), bfp_mean_exp,yerr=yerr_array_bfp,fmt="v",mfc="#CDE1BF",mec="#252525",color="#252525", lw=1.5)



# labels - 
PyPlot.plot(figsize=(4,4))
PyPlot.xlabel("log [Glucose](mM)", fontsize=16)
PyPlot.ylabel("[Venus] (μM)", fontsize=16)
# PyPlot.axis([-0.5,16.5,-50,1600])
PyPlot.xticks(fontsize=14)
PyPlot.yticks(fontsize=14)
PyPlot.tight_layout()
PyPlot.savefig("./plots/dose_response_plot_venus_2.pdf")