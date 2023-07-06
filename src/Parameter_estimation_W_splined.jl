##

include("Include.jl")
##
function objective_function(parameter_guess_array,time_start,time_step_size,time_stop,model_data_dictionary,exp_data_dictionary)

    # what is the host_type?
    host_type = :cell_free

    # Phase 1: parameter update =========================================================================== #
    # update the paramaters in the model data dictionary -
    # for now - lets only search over dG's -
    R = model_data_dictionary["R"]
    T_K = model_data_dictionary["T_K"]

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

    # lastly, update KL - keeping kx constant
    biophysical_constants_dictionary = model_data_dictionary["biophysical_constants_dictionary"]
    biophysical_constants_dictionary["translation_saturation_constant"] = parameter_guess_array[16]
    biophysical_constants_dictionary["transcription_saturation_constant"]=parameter_guess_array[17]
    model_data_dictionary["biophysical_constants_dictionary"] = biophysical_constants_dictionary

	# gluconate GntR binding parameters
	cuprous_parameter_dictionary = model_data_dictionary["cuprous_parameter_dictionary"]
	cuprous_parameter_dictionary["K_diss_CuSO4"] = parameter_guess_array[18]
    cuprous_parameter_dictionary["n_copper_cuer"] = parameter_guess_array[19]
	model_data_dictionary["cuprous_parameter_dictionary"] = cuprous_parameter_dictionary

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
	#print(TSIM, XSIM)
    # Phase 3:  compute simulation error ================================================================== #
    # compute the error - we need to do a bunch of interpolation -


    error_term_array = zeros(3,1) 

	tsim_exp_protein = exp_data_dictionary["prot_data_array"][:,1]

    # mRNA Venus -


	tsim_exp_mRNA = exp_data_dictionary["mRNA_data_array"][:,1]
    mRNA_Venus_exp= exp_data_dictionary["mRNA_data_array"][:,2]
    itp_Venus_mRNA =  Interpolations.LinearInterpolation(TSIM, 1000*XSIM[:,4]);
    A = DataInterpolations.CubicSpline(mRNA_Venus_exp,tsim_exp_mRNA)
	mRNA_Venus_exp = A.(tsim_exp_protein)
    mRNA_Venus_sim = itp_Venus_mRNA[tsim_exp_protein]

    error_vector_1 = mRNA_Venus_exp .- mRNA_Venus_sim
    error_term_array[1] = (transpose(error_vector_1)*error_vector_1)


   #mRNA CueR


   mRNA_CueR_exp= exp_data_dictionary["mRNA_data_array"][:,4]
   itp_CueR_mRNA =  Interpolations.LinearInterpolation(TSIM, 1000*XSIM[:,3]);
   A = DataInterpolations.CubicSpline(mRNA_CueR_exp,tsim_exp_mRNA)
   mRNA_CueR_exp = A.(tsim_exp_protein)
   mRNA_CueR_sim = itp_CueR_mRNA[tsim_exp_protein]

   error_vector_2 = mRNA_CueR_exp .- mRNA_CueR_sim
   error_term_array[2] = (transpose(error_vector_2)*error_vector_2)


    #create fake CueR protein data = venus protein data
    # Venus protein -



    # Venus protein -
    itp_Venus_protein =  Interpolations.LinearInterpolation(TSIM, XSIM[:,6]);
    protein_Venus_sim = itp_Venus_protein[tsim_exp_protein]
    protein_Venus_exp=  exp_data_dictionary["prot_data_array"][:,2]
    itp_CueR_protein =  Interpolations.LinearInterpolation(TSIM, XSIM[:,5]);
    protein_CueR_sim = itp_CueR_protein[tsim_exp_protein]



    error_vector_3a= 1*protein_Venus_exp .- protein_Venus_sim
    error_vector_3b= 1*protein_Venus_exp .- protein_CueR_sim
   
    error_term_array[3] = transpose(error_vector_3a)*error_vector_3a + transpose(error_vector_3b)*error_vector_3b
 
    # return -
    return error_term_array
end

# Evaluates the objective function values -
function local_refienment_step(path_to_data_dir, parameter_array; sigma=0.05, iteration_max=100)

    # inner functions -
    function _compute_error_total(objective_array,W)
        value = transpose(objective_array)*W*objective_array
        return value[1]
    end

    # initialize -
    number_of_parameters = length(parameter_array)
    BIG = 1e10

    # load the default data_dictionary -
    time_start = 0.0
    time_stop = 16.0
    time_step_size = 0.01

    # what is the host_type?
    host_type = :cell_free

    # path to parameters -
    path_to_biophysical_constants_file = "./CellFree.json"

    # wght array -
    W = diagm(ones(3)) 


    # Load the data dictionary (uses the default biophysical_constants file)
    default_data_dictionary = build_data_dictionary((time_start,time_stop,time_step_size), path_to_biophysical_constants_file, host_type)

    # Set some values for the weights, gene lengths etc -
    # model_data_dictionary = customize_data_dictionary(default_data_dictionary,host_type)
	model_data_dictionary = default_data_dictionary

    # load the experimental data -
    exp_data_dictionary = load_experimental_data_dictionary(path_to_data_dir)

    # setup the functions -
    OF(P) = objective_function(P,time_start,time_step_size,time_stop,model_data_dictionary,exp_data_dictionary)

    # calculate the starting error -
    parameter_array_best = parameter_array
    error_array = BIG*ones(4)
    error_array[1] = _compute_error_total(OF(parameter_array_best), W)

    # main refinement loop -
    iteration_counter = 1
    while (iteration_counter<iteration_max)

        # take a step up -
        parameter_up = parameter_array_best.*(1 .+ sigma*rand(number_of_parameters))
        parameter_up = check_parameter_bounds(parameter_up)

        # take a step down -
        parameter_down = parameter_array_best.*(1 .- sigma*rand(number_of_parameters))
        parameter_down = check_parameter_bounds(parameter_down)

        # Evaluate the obj function -
        error_array[2] = _compute_error_total(OF(parameter_up),W)
        error_array[3] = _compute_error_total(OF(parameter_down),W)

        # Calculate a correction factor -
        a = error_array[2]+error_array[3] - 2.0*error_array[1]
        parameter_corrected = parameter_array_best
        if (a>0.0)
            amda = -0.5*(error_array[3] - error_array[2])/a
            parameter_corrected = parameter_array_best .+ amda*rand(number_of_parameters)
            parameter_corrected = check_parameter_bounds(parameter_corrected)
            error_array[4] = _compute_error_total(OF(parameter_corrected), W)
        end

        # Which step has the min error?
        min_index = argmin(error_array)
        if (min_index == 1)
            parameter_array_best = parameter_array_best
        elseif (min_index == 2)
            parameter_array_best = parameter_up
        elseif (min_index == 3)
            parameter_array_best = parameter_down
        elseif (min_index == 4)
            parameter_array_best = parameter_corrected
        end

        # Update the local error
        error_array[1] = error_array[min_index]

        @show iteration_counter,error_array[min_index]

        # update local counter -
        iteration_counter = iteration_counter + 1
    end

    return parameter_array_best
end

function check_parameter_bounds(parameter_array)

    pvec_bounds= [

     
        # dG's -
        -10  -8    	;   # 2     W_t7rnap 8.5 to 8
        -17  -16    	;   # 2     W_apo 16.5 to 16
        -18  -17;           #w_holo 

        # binding parameters -
	    2.0  2.5            ;# 4     n_CueR_OCueR 2 to 2.3
        0.1  0.25        ;   # 5     K_CueR_OCueR_apo -Um 0.1 to 0.13
        0.04 0.08        ;   # 6     K_CueR_OCueR_holo um 0.035 to 0.045

        # time constants - units h
		0.001 100.0         ;	# 7	 mRNA_CueR 0.045 to 0.06
	    0.001 100.0         ;	# 8	 mRNA_Venus 18 to 22
	    0.001 100.0         ;	# 9	 protein_CueR 55 to 70
		0.001 100.0         ;	# 10 protein_Venus- CHANGED this to 1.8 to 2.2

        # degradation mods - unts h
		0.001  100.0 ;	# 11	        mRNA_CueR 1.8 to 2 
	    0.001  100.0	    ;	# 12	        mRNA_Venus 3.8 to 4.5
        0.001  100.0;	# 13	        protein_CueR 2.8 to 3
        0.001  100.0    ;	# 14	        protein_Venus 1.4 to 1.6

         # w -
        4.0  10.0           ;   # 15  units h  translation capacity half-life 5 to 7.5 

        # KL value -
        1 200.0  #change it to 18 to 25 uM

        #KX value

        0.05 0.10      ; #17 KX  in mUm  0.05 to 0.06

		# K_copper_CueR-
		10 50				; # 

        1.5 2.0  #N CUER 18 2.6 to 3.2
    ];
	
    pvec_initial = parameter_array

    # check bounds -
    number_of_parameters = length(pvec_initial)
    for parameter_index = 1:number_of_parameters

        # what is the parameter value?
        p_i = pvec_initial[parameter_index]

        # is p_i outside of the bounds?
        lb_value = pvec_bounds[parameter_index,1]
        ub_value = pvec_bounds[parameter_index,2]

        if (p_i<lb_value)
            pvec_initial[parameter_index,1] = lb_value
        end

        if (p_i>ub_value)
            pvec_initial[parameter_index,1] = ub_value
        end
    end

    # return -
    return pvec_initial
end

function neighbor_function(parameter_array; sigma=0.05)

    # setup -
    number_of_parameters = length(parameter_array)

    # calculate new parameter array -
    new_parameter_array = parameter_array.*(1 .+ sigma*randn(number_of_parameters))

    # check the bounds and return -
    return check_parameter_bounds(new_parameter_array)
end

function cooling_function(temperature)

  # define my new temperature -
  alpha = 0.9
  return alpha*temperature
end

function acceptance_probability_function(rank_array,temperature)
    return (exp(-rank_array[end]/temperature))
end

function main(path_to_data_dir::String, initial_parameter_array::Array{Float64,1}; rank_cutoff::Int64=4, maximum_number_of_iterations::Int64=100)

    # load the default data_dictionary -
    time_start = 0.0
    time_stop = 16.0
    time_step_size = 0.01

    # what is the host_type?
    host_type = :cell_free
    number_of_parameters = length(initial_parameter_array)

    # path to parameters -
    path_to_biophysical_constants_file = "./CellFree.json"

    # Load the data dictionary (uses the default biophysical_constants file)
    default_data_dictionary = build_data_dictionary((time_start,time_stop,time_step_size), path_to_biophysical_constants_file, host_type)

    # Set some values for the weights, gene lengths etc -
    # model_data_dictionary = customize_data_dictionary(default_data_dictionary,host_type)
	model_data_dictionary = deepcopy(default_data_dictionary)
    # load the experimental data -
    exp_data_dictionary = load_experimental_data_dictionary(path_to_data_dir)

    # setup the functions -
    OF(P) = objective_function(P,time_start,time_step_size,time_stop,model_data_dictionary,exp_data_dictionary)
    NF(P) = neighbor_function(P;sigma=0.05)

    # make call to POETs -
    (EC,PC,RA) = estimate_ensemble(OF,NF,acceptance_probability_function,cooling_function,initial_parameter_array;rank_cutoff=rank_cutoff,maximum_number_of_iterations=maximum_number_of_iterations)

    # return -
    return (EC,PC,RA)
end

##

 # setup initial condition vector -
pvec_initial = [

	# dG's -
	-8.0    	;   # 1   W_T7RNAP -8 
	-16.0   ;   # 2        W_apo -16 
	-17.0  	;   # 3       W_holo -17
	

	# binding parameters -
	2           ;   # 4     n_CueR_OCueR 2
	106E-3         ;   # 5    K_CueR_OCueR_apo uM 106e-3
    40E-3         ;   # 6    K_CueR_OCueR_holo uM 40e-3


	# time constants -
	1         ;	# 7	    mRNA_CueR 0.05
	1         ;	# 8	    mRNA_Venus 20
	1         ;	# 9	    protein_CueR 60
	1        ;	# 10	    protein_Venus 2

	# degradation mods -
	1	        ;	# 11	    mRNA_CueR 2
	1 	    ;	# 12	    mRNA_Venus 4
	1 	    ;	# 13	    protein_CueR 3
	1 	    ;	# 14	    protein_Venus 1.5

	 # w -
	6.0           ;   # 15   translation capacity half-life 6

	# KL value -
	22.8         ;   # 16   KL in muM 22.8

    #Kx value

    0.05  # 17


	# K_copper_CueR-
	25				; # 18 uM

    1.5 #N COPPER CUER 19


    # # transcription capacity terms-
    # 10; # decay (hours)
    # 0.5; # slope

    # # translation capacity terms-
    # 3.0; # decay (hours)
    # 0.5; # slope
];

# setup -
path_to_data_dir = "$(pwd())"
pV = neighbor_function(pvec_initial; sigma=0.25)
EC = 0
PC = 0
RA = 0

##

# execute -
number_of_trials = 10
for trial_index = 1:number_of_trials

    global pV
    global EC
    global PC
    global RA


   # do a local step -
    if (mod(trial_index,2) == 0)

        # find the lowest score pV -
        sum_error_array = sum(EC,dims=1)
        best_p_index = argmin(vec(sum_error_array))
        pV_best = PC[:,best_p_index]

        # local refine -
        pV = local_refienment_step(path_to_data_dir, pV_best; iteration_max=50)
    end

    # main -
    (EC,PC,RA) = main(path_to_data_dir, vec(pV); rank_cutoff=4,maximum_number_of_iterations=50)

    # dump results to disk -
    local  fname = "./simulated/POETS/RA_T$(trial_index).dat"
    writedlm(fname,RA)
    fname = "./simulated/POETS/EC_T$(trial_index).dat"
    writedlm(fname,EC)
    fname = "./simulated/POETS/PC_T$(trial_index).dat"
    writedlm(fname,PC)

    @show trial_index
end