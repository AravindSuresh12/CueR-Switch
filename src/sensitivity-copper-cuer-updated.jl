include("Include.jl")
# mean center -
function mean_center_array(results_array::Array{Float64,2})::Array{Float64,2}

    # get the size -
    (NR,NC) = size(results_array)
    scaled_array = zeros(NR,NC)

    for col_index = 1:NC

        data_col = results_array[:,col_index]
        mu_value = mean(data_col)
        std_value = std(data_col)

        for row_index = 1:NR
            scaled_array[row_index,col_index] = (data_col[row_index] - mu_value)/(std_value)
        end
    end

    return scaled_array
end

# computes the model performance -
function model_performance(parameter_guess_array,index)

    # what is the host_type?
    host_type = :cell_free

    # load the default data_dictionary -
    time_start = 0.0
    time_stop = 16.0
    time_step_size = 0.01

    # what is the host_type?
    host_type = :cell_free

    path_to_biophysical_constants_file = "./CellFree.json"

    # Load the data dictionary (uses the default biophysical_constants file)
    default_data_dictionary = build_data_dictionary((time_start,time_stop,time_step_size), path_to_biophysical_constants_file, host_type)

    # Set some values for the weights, gene lengths etc -
    model_data_dictionary = deepcopy(default_data_dictionary)

     # compute W -
     tmp_W_array = Float64[]
     for index = 1:3
         parameter_guess = parameter_guess_array[index]
         value = exp(-1*(parameter_guess))
         push!(tmp_W_array,value)
     end
 
     # update the control W's -
     control_parameter_dictionary = model_data_dictionary["control_parameter_dictionary"]
     control_parameter_dictionary["W_T7RNAP"] = tmp_W_array[1]
     control_parameter_dictionary["W_apo"]=tmp_W_array[2]
     control_parameter_dictionary["W_holo"]=tmp_W_array[3]
     model_data_dictionary["control_parameter_dictionary"] = control_parameter_dictionary
     #print(control_parameter_dictionary)
 
     binding_parameter_dictionary = model_data_dictionary["binding_parameter_dictionary"]
 
     binding_parameter_dictionary["n_CueR_OCueR"] =parameter_guess_array[4]
     binding_parameter_dictionary["K_CueR_OCueR_apo"] =parameter_guess_array[5]
     binding_parameter_dictionary["K_CueR_OCueR_holo"] =parameter_guess_array[6]
 
     model_data_dictionary["binding_parameter_dictionary"] = binding_parameter_dictionary
     #print(binding_parameter_dictionary)
     # time constant modifier -
 
     time_constant_modifier_array = [
         0.0							;	# 1	CueR
         0.0							;	# 2	Venus
         parameter_guess_array[7]	;	# 4	mRNA_CueR
         parameter_guess_array[8]	;	# 5	mRNA_Venus
         parameter_guess_array[9]	;	# 7	protein_CueR
         parameter_guess_array[10]	;	# 8	protein_Venus
     ]
 
     model_data_dictionary["time_constant_modifier_array"] = time_constant_modifier_array
 
     # setup degradation_modifier_array -
 
     degradation_modifier_array = [
         0.0	;	# 1	GntR
         0.0	;	# 2	Venus
         parameter_guess_array[11]	;	# 4	mRNA_CueR
         parameter_guess_array[12]	;	# 5	mRNA_Venus
         parameter_guess_array[13]	;	# 7	protein_CueR
         parameter_guess_array[14]	;	# 8	protein_Venus
     ]
 
     model_data_dictionary["degradation_modifier_array"] = degradation_modifier_array
 
     # update the translation time -
     model_data_dictionary["half_life_translation_capacity"] = parameter_guess_array[15]
 
     # lastly, update KL -
     biophysical_constants_dictionary = model_data_dictionary["biophysical_constants_dictionary"]
     biophysical_constants_dictionary["translation_saturation_constant"] = parameter_guess_array[16]
     biophysical_constants_dictionary["transcription_saturation_constant"] = parameter_guess_array[17]
     model_data_dictionary["biophysical_constants_dictionary"] = biophysical_constants_dictionary
 
     # gluconate GntR binding parameters
     cuprous_parameter_dictionary = model_data_dictionary["cuprous_parameter_dictionary"]
     cuprous_parameter_dictionary["K_diss_CuSO4"] = parameter_guess_array[18]
     cuprous_parameter_dictionary["n_copper_cuer"] = parameter_guess_array[19]
 
     # grab defaults -
     species_symbol_type_array = model_data_dictionary["species_symbol_type_array"]
     protein_coding_length_array = model_data_dictionary["protein_coding_length_array"]
     gene_coding_length_array = model_data_dictionary["gene_coding_length_array"]
     time_constant_modifier_array = model_data_dictionary["time_constant_modifier_array"]
     initial_condition_array = model_data_dictionary["initial_condition_array"]
 
     # # get gene IC -
     idx_gene = findall(x->x==:gene,species_symbol_type_array)
     gene_abundance_array = initial_condition_array[idx_gene]
 
     # Precompute the translation parameters -
     translation_parameter_array = precompute_translation_parameter_array(biophysical_constants_dictionary, protein_coding_length_array, time_constant_modifier_array,host_type)
     model_data_dictionary["translation_parameter_array"] = translation_parameter_array
 
     # Precompute the kinetic limit of transcription -
     transcription_kinetic_limit_array = precompute_transcription_kinetic_limit_array(biophysical_constants_dictionary, gene_coding_length_array, gene_abundance_array, time_constant_modifier_array, host_type)
     model_data_dictionary["transcription_kinetic_limit_array"] = transcription_kinetic_limit_array
 
     # Dilution degrdation matrix -
     dilution_degradation_matrix = build_dilution_degradation_matrix(biophysical_constants_dictionary, species_symbol_type_array, degradation_modifier_array)
     model_data_dictionary["dilution_degradation_matrix"] = dilution_degradation_matrix
   # ===================================================================================================== #
   #print(model_data_dictionary)
   # Phase 2:  solve model equations ===================================================================== #
   # solve the balance equations -
   (TSIM,XSIM) = SolveBalances(time_start,time_stop,time_step_size,model_data_dictionary)
    # ===================================================================================================== #

    # Phase 3: compute the model performance metrics ====================================================== #
    p_Venus_AUC = integrate(TSIM,XSIM[:,index])
    # ===================================================================================================== #

    # return the performance_array -
    return p_Venus_AUC
end

function main(path_to_ensemble_file::String,index)

    # setup the sensitivity function -
    SF(P) = model_performance(P,index)

    # setup ranges -
    sample_bounds_array = Array{Tuple,1}()
    ensemble_array = readdlm(path_to_ensemble_file)
    (number_of_parameters,number_of_trials) = size(ensemble_array)
    for parameter_index = 1:(number_of_parameters)

        # get row of parameters -
        parameter_row = ensemble_array[parameter_index,:]
        min_value = minimum(parameter_row)
        max_value = maximum(parameter_row)

        # create the tuple -
        tmp_tuple = (min_value,max_value)

        # cache -
        push!(sample_bounds_array,tmp_tuple)
    end

    @show index

    # do the global sensitivity analysis -
    sensitivity_results = GlobalSensitivity.gsa(SF,Morris(total_num_trajectory=10000,num_trajectory=1000),sample_bounds_array)

    # return -
    return sensitivity_results
end

# setup paths -
path_to_ensemble_file = "./simulated/POETS/PC_T5.dat"

# compute a sensitivity array for the AUC of each species -
species_index_array = [3 5 4 6]  #CueR mRNA, CueR Protein, Venus mRNA, Venus protein
number_of_species = length(species_index_array)
number_of_parameters = 19
results_array = zeros(number_of_parameters,1)
for species_index in species_index_array

    global results_array

    # conduct senstivity analysis -
    sensitivity_results = main(path_to_ensemble_file,species_index)

    # get the μ and σ^2
    mu = sensitivity_results.means_star
    
    var = sensitivity_results.variances

    @show mu, var

    results_array = [results_array transpose(mu) transpose(var)]

end

results_array = results_array[:,2:end]

fname = "./simulated/Sensitivity_matrix_copper.dat"
writedlm(fname,results_array)
