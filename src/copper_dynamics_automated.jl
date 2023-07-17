# ----------------------------------------------------------------------------------- #
# Copyright (c) 2019 Varnerlab
# Robert Frederick School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ----------------------------------------------------------------------------------- #

# include -
include("Include.jl")


time_start = 0.0
time_stop = 16.0
time_step_size = 0.01

function update_model_dictionary(parameter_array,default_model_dictionary)

    # what is the host_type?
	host_type = :cell_free

	# path to parameters -
	path_to_biophysical_constants_file = "./CellFree.json"

	# Load the data dictionary (uses the default biophysical_constants file)
	data_dictionary = build_data_dictionary((time_start,time_stop,time_step_size), path_to_biophysical_constants_file, host_type)
	#print(data_dictionary)

	R = data_dictionary["R"]
	T_K = data_dictionary["T_K"]

	control_parameter_dictionary = data_dictionary["control_parameter_dictionary"]

	control_parameter_dictionary["W_T7RNAP"] = exp(-1*(parameter_array[1]))
	control_parameter_dictionary["W_apo"]=  exp(-1*(parameter_array[2]))
	control_parameter_dictionary["W_holo"]= exp(-1*(parameter_array[3]))
	
	binding_parameter_dictionary = data_dictionary["binding_parameter_dictionary"]
	binding_parameter_dictionary["n_CueR_OCueR"] =parameter_array[4]
	binding_parameter_dictionary["K_CueR_OCueR_apo"] =parameter_array[5]
	binding_parameter_dictionary["K_CueR_OCueR_holo"] =parameter_array[6]
	data_dictionary["binding_parameter_dictionary"] = binding_parameter_dictionary
	
	time_constant_modifier_array = [
		0.0							;	# 1	CueR
		0.0							;	# 2	Venus
		parameter_array[7]	        ;	# 4	mRNA_CueR
		parameter_array[8]       ;	    # 5	mRNA_Venus
		parameter_array[9]	        ;	# 7	protein_CueR
		parameter_array[10]	        ;	# 8	protein_Venus
	]
	
	data_dictionary["time_constant_modifier_array"] = time_constant_modifier_array
	
	degradation_modifier_array = [
		0.0	;	# 1	GntR
		0.0	;	# 2	Venus
		parameter_array[11]	;	# 4	mRNA_GntR
		parameter_array[12]	;	# 5	mRNA_Venus
		parameter_array[13]	;	# 7	protein_CueR
		parameter_array[14]	;	# 8	protein_Venus
	]
	
	data_dictionary["degradation_modifier_array"] = degradation_modifier_array
	
	# update the translation time -
	data_dictionary["half_life_translation_capacity"] = parameter_array[15]
	biophysical_constants_dictionary = data_dictionary["biophysical_constants_dictionary"]
	
    biophysical_constants_dictionary["translation_saturation_constant"] = parameter_array[16]
    biophysical_constants_dictionary["transcription_saturation_constant"]=parameter_array[17]
 

	# gluconate GntR binding parameters
	cuprous_parameter_dictionary = data_dictionary["cuprous_parameter_dictionary"]
	cuprous_parameter_dictionary["K_diss_CuSO4"] = parameter_array[18]
    cuprous_parameter_dictionary["n_copper_cuer"] = parameter_array[19]

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

    # return -
    return data_dictionary
end

copper_array=[0,5,10,50,100]

l=length(copper_array)
for k = 1:l # for now just do it for 5 


function main(path_to_ensemble_file::String, path_to_sim_dir::String; sample_fraction::Float64=1.0)

    # load the default data_dictionary -
    time_start = 0.0
    time_stop = 16.0
    time_step_size = 0.01

    # what is the host_type?
    host_type = :cell_free

    # path to parameters -
    path_to_biophysical_constants_file = "./CellFree.json"

    # Load the data dictionary (uses the default biophysical_constants file)
    default_data_dictionary = build_data_dictionary((time_start,time_stop,time_step_size), path_to_biophysical_constants_file, host_type)

    # customize -
    customized_data_dictionary = deepcopy(default_data_dictionary)

    # load the ensemble file -
    parameter_ensemble_array = readdlm(path_to_ensemble_file)

    # sort the array -
    parameter_ensemble_array = sort_parameter_ensemble_array(parameter_ensemble_array)

    # what is the size of the ensemble -
    (number_of_parameters, number_of_samples) = size(parameter_ensemble_array)

    # take the top sample fraction -
    number_of_samples = round(Int,sample_fraction*number_of_samples)
    for sample_index = 1:number_of_samples

        # grab the parameter array -
        local_parameter_array = parameter_ensemble_array[:,sample_index]
		# print(local_parameter_array)
        # update the default data_dictionary to reflect the new parameters -
        model_data_dictionary = update_model_dictionary(local_parameter_array, customized_data_dictionary)

		cuprous_parameter_dictionary = model_data_dictionary["cuprous_parameter_dictionary"]
		cuprous_parameter_dictionary["copper_salt_conc"]=copper_array[k]
		model_data_dictionary["cuprous_parameter_dictionary"] = cuprous_parameter_dictionary

		# solve the model equations -
        (T,X) = SolveBalances(time_start,time_stop,time_step_size,model_data_dictionary) #the dictionary needs to update here- not updating

        # dump -
        data_array = [T X]
        filename = "$(path_to_sim_dir)/simulation-P$(sample_index).dat"
        writedlm(filename, data_array)

        # give user some notification -
        @show sample_index
    end
end

# setup paths -
local path_to_sim_file = "$(pwd())/simulated/copper_dynamics/$(copper_array[k])uM"
local path_to_ensemble_file = "$(pwd())/simulated/POETS/PC_T5.dat"




# call -
main(path_to_ensemble_file, path_to_sim_file; sample_fraction = 1.0)

end

