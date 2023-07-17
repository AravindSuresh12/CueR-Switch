include("Include.jl")

#AUC AND MODEL PERFORMANCE CHARECTARISTICS


function predict(path_to_ensemble_file,copper_conc)

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


R = data_dictionary["R"]
T_K = data_dictionary["T_K"]
# compute W -

poets= readdlm(path_to_ensemble_file)
parameter_guess_array = mean(poets,dims=2)
poets_params= parameter_guess_array


control_parameter_dictionary = data_dictionary["control_parameter_dictionary"]

control_parameter_dictionary["W_T7RNAP"] =exp(-1*(poets_params[1]))
control_parameter_dictionary["W_apo"]=exp(-1*(poets_params[2]))
control_parameter_dictionary["W_holo"]=exp(-1*(poets_params[3]))

binding_parameter_dictionary = data_dictionary["binding_parameter_dictionary"]
binding_parameter_dictionary["n_CueR_OCueR"] =poets_params[4]
binding_parameter_dictionary["K_CueR_OCueR_apo"]=poets_params[5]
binding_parameter_dictionary["K_CueR_OCueR_holo"]=poets_params[6]
data_dictionary["binding_parameter_dictionary"] = binding_parameter_dictionary

time_constant_modifier_array = [
    0.0							;# 1	CueR
    0.0							;# 2	Venus
    poets_params[7]	        ;	# 4	mRNA_CueR
    poets_params[8]       ;	    # 5	mRNA_Venus
    poets_params[9]	        ;	# 7	protein_CueR
    poets_params[10]	        ;# 8 protein_Venus- should be 0.5
]

data_dictionary["time_constant_modifier_array"] = time_constant_modifier_array

degradation_modifier_array = [
    0.0	;	# 1	GntR
    0.0	;	# 2	Venus
    poets_params[11]	;	# 4	mRNA_GntR
    poets_params[12]	;	# 5	mRNA_Venus
    poets_params[13]	;	# 7	protein_CueR
    poets_params[14]	;	# 8	protein_Venus
]

data_dictionary["degradation_modifier_array"] = degradation_modifier_array

# update the translation time -
data_dictionary["half_life_translation_capacity"] = poets_params[15]

biophysical_constants_dictionary = data_dictionary["biophysical_constants_dictionary"]
biophysical_constants_dictionary["translation_saturation_constant"] = poets_params[16]
biophysical_constants_dictionary["transcription_saturation_constant"] = poets_params[17]
data_dictionary["biophysical_constants_dictionary"] = biophysical_constants_dictionary

# Copper CueR binding parameters
cuprous_parameter_dictionary = data_dictionary["cuprous_parameter_dictionary"]
cuprous_parameter_dictionary["copper_salt_conc"]=copper_conc
cuprous_parameter_dictionary["K_diss_CuSO4"] = poets_params[18]
cuprous_parameter_dictionary["n_copper_cuer"] = poets_params[19]


data_dictionary["cuprous_parameter_dictionary"] = cuprous_parameter_dictionary

species_symbol_type_array = data_dictionary["species_symbol_type_array"]
protein_coding_length_array = data_dictionary["protein_coding_length_array"]
gene_coding_length_array = data_dictionary["gene_coding_length_array"]
time_constant_modifier_array = data_dictionary["time_constant_modifier_array"]
initial_condition_array = data_dictionary["initial_condition_array"]

# # get gene IC -
idx_gene = findall(x->x==:gene,species_symbol_type_array)
gene_abundance_array = initial_condition_array[idx_gene]

# Precompute the translation parameters -
translation_parameter_array = precompute_translation_parameter_array(biophysical_constants_dictionary, protein_coding_length_array, time_constant_modifier_array,host_type)
data_dictionary["translation_parameter_array"] = translation_parameter_array

# Precompute the kinetic limit of transcription -
transcription_kinetic_limit_array = precompute_transcription_kinetic_limit_array(biophysical_constants_dictionary, gene_coding_length_array, gene_abundance_array, time_constant_modifier_array, host_type)
data_dictionary["transcription_kinetic_limit_array"] = transcription_kinetic_limit_array

# Dilution degrdation matrix -
dilution_degradation_matrix = build_dilution_degradation_matrix(biophysical_constants_dictionary, species_symbol_type_array, degradation_modifier_array)
data_dictionary["dilution_degradation_matrix"] = dilution_degradation_matrix

# Solve the model equations -
(T,X) = SolveBalances(time_start,time_stop,time_step_size,data_dictionary)

return (last(X[:,6]))

end

path_to_ensemble_file = "./simulated/POETS/PC_T5.dat"


