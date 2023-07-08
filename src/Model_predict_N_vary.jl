#model prediction capability for N=5, 20 and 100 - without ensemble fits
include("Include.jl")

# Script to solve the balance equations -
time_start = 0.0
time_stop = 16.0
time_step_size = 0.01

# what is the host_type?
host_type = :cell_free

# path to parameters -
path_to_biophysical_constants_file = "./CellFree.json"

# Load the data dictionary (uses the default biophysical_constants file)
data_dictionary = build_data_dictionary((time_start,time_stop,time_step_size), path_to_biophysical_constants_file, host_type)
#print(data_dictionary)

R = data_dictionary["R"]
T_K = data_dictionary["T_K"]
# compute W -

PC_N5=readdlm("./simulated/POETS_N_Vary/N=5/PC_T4.dat") 
PC_N20=readdlm("./simulated/POETS/PC_T2.dat") 
PC_N100=readdlm("./simulated/POETS_N_Vary/N=100/PC_T10.dat") 

average_params_N5 = mean(PC_N5,dims=2)
poets_params_N5=average_params_N5

average_params_N20 = mean(PC_N20,dims=2)
poets_params_N20=average_params_N20

average_params_N100 = mean(PC_N100,dims=2)
poets_params_N100=average_params_N100

new_array=[poets_params_N5,poets_params_N20,poets_params_N100]


storage_array=[] #used for storing the dose response values for range

for element in 1:length(new_array)
    
# Update the data dictionary
local control_parameter_dictionary = data_dictionary["control_parameter_dictionary"]

control_parameter_dictionary["W_T7RNAP"] =exp(-1*(new_array[element][1]))
control_parameter_dictionary["W_apo"]=exp(-1*(new_array[element][2]))
control_parameter_dictionary["W_holo"]=exp(-1*(new_array[element][3]))

local binding_parameter_dictionary = data_dictionary["binding_parameter_dictionary"]
binding_parameter_dictionary["n_CueR_OCueR"] =new_array[element][4]
binding_parameter_dictionary["K_CueR_OCueR_apo"]=new_array[element][5]
binding_parameter_dictionary["K_CueR_OCueR_holo"]=new_array[element][6]
data_dictionary["binding_parameter_dictionary"] = binding_parameter_dictionary

time_constant_modifier_array = [
    0.0							;# 1	CueR
    0.0							;# 2	Venus
    new_array[element][7]	        ;	# 4	mRNA_CueR
    new_array[element][8]       ;	    # 5	mRNA_Venus
    new_array[element][9]	        ;	# 7	protein_CueR
    new_array[element][10]	        ;# 8 protein_Venus
]

data_dictionary["time_constant_modifier_array"] = time_constant_modifier_array

local degradation_modifier_array = [
    0.0	;	# 1	GntR
    0.0	;	# 2	Venus
    new_array[element][11]	;	# 4	mRNA_GntR
    new_array[element][12]	;	# 5	mRNA_Venus
    new_array[element][13]	;	# 7	protein_CueR
    new_array[element][14]	;	# 8	protein_Venus
]

local data_dictionary["degradation_modifier_array"] = degradation_modifier_array

# update the translation time -
local data_dictionary["half_life_translation_capacity"] = new_array[element][15]

local biophysical_constants_dictionary = data_dictionary["biophysical_constants_dictionary"]
biophysical_constants_dictionary["translation_saturation_constant"] = new_array[element][16]
biophysical_constants_dictionary["transcription_saturation_constant"] = new_array[element][17]
data_dictionary["biophysical_constants_dictionary"] = biophysical_constants_dictionary

# Copper CueR binding parameters
local cuprous_parameter_dictionary = data_dictionary["cuprous_parameter_dictionary"]
cuprous_parameter_dictionary["K_diss_CuSO4"] = new_array[element][18]
cuprous_parameter_dictionary["n_copper_cuer"] = new_array[element][19]


local data_dictionary["cuprous_parameter_dictionary"] = cuprous_parameter_dictionary

local species_symbol_type_array = data_dictionary["species_symbol_type_array"]
local protein_coding_length_array = data_dictionary["protein_coding_length_array"]
local gene_coding_length_array = data_dictionary["gene_coding_length_array"]
local time_constant_modifier_array = data_dictionary["time_constant_modifier_array"]
local initial_condition_array = data_dictionary["initial_condition_array"]

# # get gene IC -
local idx_gene = findall(x->x==:gene,species_symbol_type_array)
local gene_abundance_array = initial_condition_array[idx_gene]

# Precompute the translation parameters -
local translation_parameter_array = precompute_translation_parameter_array(biophysical_constants_dictionary, protein_coding_length_array, time_constant_modifier_array,host_type)
data_dictionary["translation_parameter_array"] = translation_parameter_array

# Precompute the kinetic limit of transcription -
local transcription_kinetic_limit_array = precompute_transcription_kinetic_limit_array(biophysical_constants_dictionary, gene_coding_length_array, gene_abundance_array, time_constant_modifier_array, host_type)
data_dictionary["transcription_kinetic_limit_array"] = transcription_kinetic_limit_array

# Dilution degrdation matrix -
local dilution_degradation_matrix = build_dilution_degradation_matrix(biophysical_constants_dictionary, species_symbol_type_array, degradation_modifier_array)
data_dictionary["dilution_degradation_matrix"] = dilution_degradation_matrix

local copper_sim=[0,5,10,50,100]

for copper_conc in copper_sim
    local cuprous_parameter_dictionary["copper_salt_conc"]=copper_conc
    local (T,X) = SolveBalances(time_start,time_stop,time_step_size,data_dictionary)
    local end_value=last(X[:,6])
    push!(storage_array,end_value)
end

end

N_array_5=storage_array[1:5]
N_array_20=storage_array[6:10]
N_array_100=storage_array[11:15]

dose_data = CSV.read("./data/dose_response.csv",DataFrame)
exp_val= dose_data[!,"Venus(uM)"]
exp_std= dose_data[!,"STERR(uM)"]
yerr_array = transpose([exp_std exp_std])

fig=figure()

subplot(1,3,1)
PyPlot.scatter(copper_sim, exp_val)
PyPlot.plot(copper_sim, N_array_5)
PyPlot.xlabel("Copper_concentration(μM)", fontsize=5)
PyPlot.ylabel("Venus_output(μM)", fontsize=5)
PyPlot.yticks(fontsize=5)


subplot(1,3,2)
PyPlot.scatter(copper_sim, exp_val)
PyPlot.plot(copper_sim, N_array_20)
PyPlot.xlabel("Copper_concentration(μM)", fontsize=5)
PyPlot.ylabel("Venus_output(μM)", fontsize=5)
PyPlot.yticks(fontsize=5)


subplot(1,3,3)
PyPlot.scatter(copper_sim, exp_val)
PyPlot.plot(copper_sim, N_array_20)
PyPlot.xlabel("Copper_concentration(μM)", fontsize=5)
PyPlot.ylabel("Venus_output(μM)", fontsize=5)
PyPlot.yticks(fontsize=5)

PyPlot.savefig("./plots/Model_pred.pdf")
