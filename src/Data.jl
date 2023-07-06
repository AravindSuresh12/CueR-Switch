# ----------------------------------------------------------------------------------- #
# Copyright (c) 2021 Varnerlab
# Robert Frederick Smith School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850
#
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
#
# ----------------------------------------------------------------------------------- #
# Function: DataDictionary
# Description: Holds simulation and model parameters as key => value pairs in a Julia Dict()
# Generated on: 2021-04-29T20:40:53.886
#
# Input arguments:
# time_start::Float64 => Simulation start time value (scalar)
# time_stop::Float64 => Simulation stop time value (scalar)
# time_step::Float64 => Simulation time step (scalar)
#
# Output arguments:
# data_dictionary::Dict{String,Any} => Dictionary holding model and simulation parameters as key => value pairs
# ----------------------------------------------------------------------------------- #
function build_data_dictionary(time_span::Tuple{Float64,Float64,Float64}, path_to_biophysical_constants_file::String = "./CellFree.json", host_type::Symbol = :cell_free)::Dict{String,Any}

	# load the biophysical_constants dictionary
	biophysical_constants_dictionary = build_biophysical_dictionary(path_to_biophysical_constants_file, host_type)

	# stoichiometric_matrix and dilution_matrix -
	stoichiometric_matrix = readdlm("./Network.dat")

	# number of states, and rates -
	(number_of_states,number_of_rates) = size(stoichiometric_matrix)

	# array of species types -
	species_symbol_type_array = [
		:gene	;	# 1	CueR
		:gene	;	# 2	Venus
		:mrna	;	# 3	mRNA_CueR
		:mrna	;	# 4	mRNA_venus
		:protein;	# 5	protein_CueR
		:protein;	# 6	protein_venus
	]

	# we need to store the species symbol array for later -
	biophysical_constants_dictionary["species_symbol_type_array"] = species_symbol_type_array

	# array of gene lengths -
	gene_coding_length_array = [
		413.0	;	# 1	CueR- will be more with our modified T7  
		730.0	;	# 2	Venus
	]

	# array of mRNA coding lengths -
	mRNA_coding_length_array = [
		gene_coding_length_array[1]	;	# 4	1	mRNA_CueR
		gene_coding_length_array[2]	;	# 5	2	mRNA_venus
	]

	# array of mRNA coding lengths -
	protein_coding_length_array = [
		round((0.33)*mRNA_coding_length_array[1])	;	# 7	1	protein_GntR
		round((0.33)*mRNA_coding_length_array[2])	;	# 8	2	protein_Venus
	
	]

	# array of intracellular gene copy numbers -
	gene_abundance_array = [
		0.007	;	# (number/cell) 1	CueR
		0.005	;	# (number/cell) 2	venus
	]

	# initial condition array -
	initial_condition_array = [
		gene_abundance_array[1]	;	# 1	CueR
		gene_abundance_array[2]	;	# 2	Venus
		0.0	;	# 3	mRNA_CueR
		0.0	;	# 4	mRNA_Venus
		0.0	;	# 5	protein_CueR
		0.0	;	# 6	protein_Venus
		# translation capacity -
		100.0	;	#7 translation capacity 
	]

	binding_parameter_dictionary = Dict{String,Float64}()
	binding_parameter_dictionary["n_CueR_OCueR"] = 1 #
	binding_parameter_dictionary["K_CueR_OCueR_apo"] = 0.106  #uM - Martella- open complex- 0.106
	binding_parameter_dictionary["K_CueR_OCueR_holo"] = 0.04  #UM  -Martella closed complex 0.4

	
	
	# Alias the control function parameters -
	control_parameter_dictionary = Dict{String,Float64}()
	control_parameter_dictionary["W_T7RNAP"] =1e4 # baseline reference of T7RNAP binding to T7promoter
	control_parameter_dictionary["W_apo"] =  1e7 # refers to apo-CueR protein binding on OCueR downstream of T7promoter
	control_parameter_dictionary["W_holo"] = 3e7 #w refers to the holo-CueR protein binding on to OCueR downstream of T7 promoter

	# Copper parameter values
	cuprous_parameter_dictionary = Dict{String,Float64}()
	cuprous_parameter_dictionary["K_diss_CuSO4"] = 18 #uM half maximal induction of the system to CuS04 in S12 extract. 
	cuprous_parameter_dictionary["copper_salt_conc"]=100 #testing of model- copper conc- uM
	cuprous_parameter_dictionary["n_copper_cuer"]=1.5 #copper binding to cuer protein

	# time constant modifiers - 
	time_constant_modifier_array = [
		0.0	;	# 1	CueR
		0.0	;	# 2	Venus
		1.0;	# 3	mRNA_CueR-
		1.0	;	# 4	mRNA_Venus
		1.0	;	# 5	protein_CueR #
		1.0;	# 6	protein_Venus- 
	]

	# degradation modifiers - #
	degradation_modifier_array = [
		0.0	;	# 1	CueR
		0.0	;	# 2	sfgfp
		1.0	;	# 3	mRNA_CueR 
		1.0	;# 4 mRNA_Venus
		1.0	;	# 5	protein_CueR 
		1.0	;	# 6	protein_Venus
	]



	# Dilution degrdation matrix -
	dilution_degradation_matrix = build_dilution_degradation_matrix(biophysical_constants_dictionary,species_symbol_type_array,degradation_modifier_array)

	# Precompute the translation parameters -
	translation_parameter_array = precompute_translation_parameter_array(biophysical_constants_dictionary, protein_coding_length_array, time_constant_modifier_array, host_type)

	# Precompute the kinetic limit of transcription -
	transcription_kinetic_limit_array = precompute_transcription_kinetic_limit_array(biophysical_constants_dictionary, gene_coding_length_array, gene_abundance_array, time_constant_modifier_array, host_type)

	# Parameter name index array -
	parameter_name_mapping_array = [
		"W_T7RNAP";  #1
		"W_apo";  #2 
		"W_holo"  #3
		"n_CueR_OCueR"; #4
		"K_CueR_OCueR_apo"; #5
		"K_CueR_OCueR_holo"; #6
		"K_diss_CuSO4" #7
		"n_copper_cuer"
		"RNAP concentration" #8
		"ribosome_concentration"	;	# 9
		"degradation_constant_mRNA"	;	# 10
		"degradation_constant_protein"	;	# 11
		"kcat_transcription"	;	# 12
		"kcat_translation"	;	# 13
		"saturation_constant_transcription"	;	# 14
		"saturation_constant_translation"; #15
		"tuning_factor" #16
	]

	# =============================== DO NOT EDIT BELOW THIS LINE ============================== #
	data_dictionary = Dict{String,Any}()
	data_dictionary["number_of_states"] = number_of_states
	data_dictionary["species_symbol_type_array"] = species_symbol_type_array
	data_dictionary["initial_condition_array"] = initial_condition_array
	data_dictionary["gene_coding_length_array"] = gene_coding_length_array
	data_dictionary["mRNA_coding_length_array"] = mRNA_coding_length_array
	data_dictionary["protein_coding_length_array"] = protein_coding_length_array
	data_dictionary["stoichiometric_matrix"] = stoichiometric_matrix
	data_dictionary["dilution_degradation_matrix"] = dilution_degradation_matrix
	data_dictionary["binding_parameter_dictionary"] = binding_parameter_dictionary
	data_dictionary["control_parameter_dictionary"] = control_parameter_dictionary
	data_dictionary["parameter_name_mapping_array"] = parameter_name_mapping_array
	data_dictionary["transcription_kinetic_limit_array"] = transcription_kinetic_limit_array
	data_dictionary["translation_parameter_array"] = translation_parameter_array
	data_dictionary["degradation_modifier_array"] = degradation_modifier_array
	data_dictionary["time_constant_modifier_array"] = time_constant_modifier_array
	data_dictionary["biophysical_constants_dictionary"] = biophysical_constants_dictionary
	data_dictionary["cuprous_parameter_dictionary"] = cuprous_parameter_dictionary

	# =============================== DO NOT EDIT ABOVE THIS LINE ============================== #
	data_dictionary["R"] = 8.314 			# J mol^-1 K^-1
	data_dictionary["T_K"] = 273.15 + 30.0 	# K
	data_dictionary["half_life_translation_capacity"] = 8 # hr

	return data_dictionary

end
