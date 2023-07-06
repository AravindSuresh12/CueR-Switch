include("Include.jl")

function main(path_to_simulation_dir::String, path_to_plot_file::String, concentration)
    # what index is mRNA venus?
    state_index = 4 #should be one less
    # how many files?
    searchdir(path,key) = filter(x->contains(x,key),readdir(path))
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

        x = sim_data_array[:,state_index+1]*1000 #bc the first array is always the time
        t = sim_data_array[:,1]

        spline_obj = DataInterpolations.CubicSpline(x,t)
		mRNA_values = spline_obj.(t_intervals)
		data_array[:, file_index] = mRNA_values

        #PyPlot.plot(t,x,color="dimgrey",alpha=0.80,lw=0.5)
    end
    # plot -
	μ = mean(data_array,dims=2)
    σ = std(data_array,dims=2)
    LB = μ .- (1.96/sqrt(1))*σ
    UB = μ .+ (1.96/sqrt(1))*σ

 
    fill_between(t_intervals, vec(UB), vec(LB), color="#5c995b", alpha=0.4)
    # Plot mean -
    PyPlot.plot(t_intervals,μ,"-",color="black",lw=2)

	# load the experimemtal data -
	# experimental_data_dictionary = load_experimental_data_dictionary(pwd())
	# plot experimental data
	mRNA_data = CSV.read("./data/mRNA_data.csv",DataFrame)

	# plot the experimemtal data -
	TEXP = mRNA_data[!,"Average_time(h)"]
	DATA = mRNA_data[!,"<Venus+CueR+Cu>(nM)"]
    STD = mRNA_data[!,"SE3(nM)"] #ABHI changed to stderr

    yerr_array = transpose([STD STD])
	PyPlot.errorbar(TEXP, DATA,yerr=yerr_array,fmt="o",mfc="#252525",mec="#252525",color="#252525", lw=1.5,ms=8)

    # labels -
	PyPlot.plot(figsize=(5,4))
	PyPlot.xlabel("Time (hr)", fontsize=20)
    PyPlot.ylabel("[Venus] (nM)", fontsize=20)
    # PyPlot.axis([-0.5,16.5,-50,1600])
    # PyPlot.xticks([0,1,2,4,8,16], fontsize=12)
    PyPlot.yticks(fontsize=22)
    PyPlot.tight_layout()
    PyPlot.savefig("$(path_to_plot_file)/mrna-venus_$(concentration)uM.pdf")
end


for conc_ in Any[100]

concentration = conc_
local path_to_simulation_dir = "$(pwd())/simulated/copper_dynamics/$(conc_)uM"
local path_to_plot_file = "$(pwd())/plots"
clf()
main(path_to_simulation_dir, path_to_plot_file,concentration)
end

  