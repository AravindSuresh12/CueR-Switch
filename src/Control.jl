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
# Function: calculate_transcription_control_array
# Description: Calculate the transcriptional control array at time t
# Generated on: 2021-04-29T20:40:54.294
#
# Input arguments:
# t::Float64 => Current time value (scalar)
# x::Array{Float64,1} => State array (number_of_species x 1)
# data_dictionary::Dict{String,Any} => Dictionary holding model parameters
#
# Output arguments:
# control_array::Array{Float64,1} => Transcriptional control array (number_of_genes x 1) at time t
# ----------------------------------------------------------------------------------- #
function calculate_transcription_control_array(t::Float64,x::Array{Float64,1},data_dictionary::Dict{String,Any})

	# initialize the control -
	control_array = zeros(2)

	# Alias the species -
	CueR = x[1]
	Venus = x[2]
	mRNA_CueR = x[3]
	mRNA_Venus = x[4]
	protein_CueR = x[5]
	protein_Venus = x[6]


	# Alias the binding parameters - 
	binding_parameter_dictionary = data_dictionary["binding_parameter_dictionary"]
	n_CueR_OCueR=binding_parameter_dictionary["n_CueR_OCueR"] 
	K_CueR_OCueR_apo=binding_parameter_dictionary["K_CueR_OCueR_apo"] 
	K_CueR_OCueR_holo=binding_parameter_dictionary["K_CueR_OCueR_holo"] 




	# Alias the control function parameters -
	control_parameter_dictionary = data_dictionary["control_parameter_dictionary"]

	W_T7RNAP= control_parameter_dictionary["W_T7RNAP"] 
	W_apo=control_parameter_dictionary["W_apo"]
	W_holo=control_parameter_dictionary["W_holo"]


	cuprous_parameter_dictionary=data_dictionary["cuprous_parameter_dictionary"]

	Cu_initial=cuprous_parameter_dictionary["copper_salt_conc"]
	K_diss_CuSO4=cuprous_parameter_dictionary["K_diss_CuSO4"] 
	n_copper_cuer=cuprous_parameter_dictionary["n_copper_cuer"]


	f_CueR= ((Cu_initial)^n_copper_cuer)/((Cu_initial)^n_copper_cuer + (K_diss_CuSO4^n_copper_cuer))


	num1= (W_T7RNAP)
	den1= 1 +num1 


	control_array[1]= num1/den1

	actor_set= [
		protein_CueR
	]

	actor=prod(actor_set)

	holo_CueR= f_CueR*actor
	apo_CueR= (1-f_CueR)*actor



	g_holo_CueR= abs(holo_CueR)^(n_CueR_OCueR)/ (abs(holo_CueR)^(n_CueR_OCueR) + (K_CueR_OCueR_holo)^(n_CueR_OCueR))

	g_apo_CueR= abs(apo_CueR)^(n_CueR_OCueR)/(abs(apo_CueR)^(n_CueR_OCueR) + abs(K_CueR_OCueR_apo)^(n_CueR_OCueR)) #THIS IS A REPRESSOR
	

	num= W_T7RNAP + W_holo*g_holo_CueR

	den= num + W_apo*g_apo_CueR


	control_array[2]=num/den

	correction_term = (x[7]/1000.0)
    control_array = control_array*correction_term

end

#
# ----------------------------------------------------------------------------------- #
# Function: calculate_translation_control_array
# Description: Calculate the translation control array at time t
# Generated on: 2021-04-29T20:40:54.447
#
# Input arguments:
# t::Float64 => Current time value (scalar)
# x::Array{Float64,1} => State array (number_of_species x 1)
# data_dictionary::Dict{String,Any} => Dictionary holding model parameters
#
# Output arguments:
# control_array::Array{Float64,1} => Translation control array (number_of_genes x 1) at time t
# ----------------------------------------------------------------------------------- #
function calculate_translation_control_array(t::Float64,x::Array{Float64,1},data_dictionary::Dict{String,Any})

	# initialize the control -
	control_array = ones(2)
	

	correction_term = (x[7]/1000.0)

    control_array = control_array*correction_term
	# # # return -
	return control_array

end
