include("Include.jl")

function main(path_to_simulation_dir::String, path_to_plot_file::String, concentration)
    # what index is prot venus?
    state_index = 6
    # how many files?
    searchdir(path,key) = filter(x->contains(x,key),readdir(path)) #function
    file_name_array = searchdir(path_to_simulation_dir, ".dat")
    number_of_trials = length(file_name_array)

	dt = 0.1 # hr
	tEND = convert(Int64,16/dt)
	t_intervals = collect(0:dt:tEND*dt)
	data_array = zeros(length(t_intervals),number_of_trials)

    # read the simulation dir -
    for (file_index,file_name) in enumerate(file_name_array)
        # load -
        file_path = "$(path_to_simulation_dir)/$(file_name)"
        sim_data_array = readdlm(file_path)

        x = sim_data_array[:,state_index+1]
        t = sim_data_array[:,1]

        spline_obj = DataInterpolations.CubicSpline(x,t)
		prot_values = spline_obj.(t_intervals)
		data_array[:, file_index] = prot_values

    
    end
    # plot -
	μ = mean(data_array,dims=2)
    σ = std(data_array,dims=2)
    LB = μ .- (1.96/sqrt(1))*σ
    UB = μ .+ (1.96/sqrt(1))*σ

    fill_between(t_intervals, vec(UB), vec(LB), color="#5c995b", alpha=0.40)
    # Plot mean -
    PyPlot.plot(t_intervals,μ,"-",color="black",lw=2)

    prot_data = CSV.read("./data/protein_data.csv",DataFrame)

	# plot the experimemtal data -
	TEXP = prot_data[!,"time(h)"]   
	DATA = prot_data[!,"Mean_0uM(uM)"] #change this for the column name for plotting specific data
    STD = prot_data[!,"Sterr_0uM(uM)"] 



	yerr_array = transpose([STD STD])
	PyPlot.errorbar(TEXP, DATA,yerr=yerr_array,fmt="o",mfc="#252525",mec="#252525",color="#252525", lw=1.5,ms=8)

    # labels -
	PyPlot.plot(figsize=(5,4))
	PyPlot.xlabel("Time (hr)", fontsize=20)
    PyPlot.ylabel("[Venus] (μM)", fontsize=20)
    PyPlot.yticks(fontsize=22)
    PyPlot.tight_layout()
    PyPlot.savefig("$(path_to_plot_file)/prot-Venus_$(concentration)uM.pdf")
end


for conc_ in Any[5]
concentration = conc_
local path_to_simulation_dir = "$(pwd())/simulated/copper_dynamics/$(conc_)uM"
local path_to_plot_file = "$(pwd())/plots"
clf()
main(path_to_simulation_dir, path_to_plot_file,concentration)
end
    
      
  