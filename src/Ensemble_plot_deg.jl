include("Include.jl")


PA=readdlm("./simulated/POETS/PC_T5.dat") #best fit
#PA=readdlm("./poets_ensemble_T_100/PC_T3.dat")
#PA=readdlm("./poets_ensemble_T_1000/PC_T10.dat")

time_start = 0.0
time_stop = 16.0
time_step_size = 0.01

# what is the host_type?
host_type = :cell_free

# path to parameters -
path_to_biophysical_constants_file = "./CellFree.json"

default_data_dictionary = build_data_dictionary((time_start,time_stop,time_step_size), path_to_biophysical_constants_file, host_type)

# Set some values for the weights, gene lengths etc -
custom_data_dictionary = deepcopy(default_data_dictionary)

# get elongation rates et al -
RNAP_elongation_rate = default_data_dictionary["biophysical_constants_dictionary"]["transcription_elongation_rate"]
Ribsome_elongation_rate = default_data_dictionary["biophysical_constants_dictionary"]["translation_elongation_rate"]
characteristic_initiation_time_transcription = default_data_dictionary["biophysical_constants_dictionary"]["characteristic_initiation_time_transcription"]
characteristic_initiation_time_translation = default_data_dictionary["biophysical_constants_dictionary"]["characteristic_initiation_time_translation"]
default_mRNA_half_life_in_hr = default_data_dictionary["biophysical_constants_dictionary"]["mRNA_half_life_in_hr"]
default_prot_half_life_in_hr = default_data_dictionary["biophysical_constants_dictionary"]["protein_half_life_in_hr"]

 # compute half lifes for mRNA_cuer -
 deg_mod_mRNA = PA[11,:]
 number_of_mRNA = 1
 (number_of_samples) = length(deg_mod_mRNA)
 mRNA_half_life_array_cuer = zeros(number_of_samples)
 kD = log(2)/(default_mRNA_half_life_in_hr)
 for sample_index = 1:number_of_samples
     kA = kD*deg_mod_mRNA[sample_index]
     half_life_value = (log(2)/kA)*60    # convert to min
     mRNA_half_life_array_cuer[sample_index] = half_life_value
 end

 # compute half lifes for mRNA_venus -
 deg_mod_mRNA = PA[12,:]
 number_of_mRNA = 1
 (number_of_samples) = length(deg_mod_mRNA)
 mRNA_half_life_array_venus = zeros(number_of_samples)
 kD = log(2)/(default_mRNA_half_life_in_hr)
 for sample_index = 1:number_of_samples
     kA = kD*deg_mod_mRNA[sample_index]
     half_life_value = (log(2)/kA)*60    # convert to min
     mRNA_half_life_array_venus[sample_index] = half_life_value
 end


 # compute half lifes for protein - cuer
 deg_mod_prot = PA[13:14,:]
 (number_of_prot,number_of_samples) = size(deg_mod_prot)
 prot_half_life_array = zeros(number_of_prot, number_of_samples)
 for prot_index = 1:number_of_prot

     kD = log(2)/(default_prot_half_life_in_hr)

     for sample_index = 1:number_of_samples
         kA = kD*deg_mod_prot[prot_index, sample_index]
         half_life_value = (log(2)/kA)*(1/24)    # convert to days
         prot_half_life_array[prot_index, sample_index] = half_life_value
     end
 end

#  fig=figure() 

subplot(2,2,1)
PyPlot.hist( (mRNA_half_life_array_cuer)', bins=50)
PyPlot.xlabel("half life mRNA CueR (min)", fontsize=5)
PyPlot.ylabel("count", fontsize=5)

subplot(2,2,2)
PyPlot.hist( (mRNA_half_life_array_venus)', bins=50)
PyPlot.xlabel("tau_factor mRNA Venus (min)", fontsize=5)
PyPlot.ylabel("count", fontsize=5)

subplot(2,2,3)
PyPlot.hist((prot_half_life_array[1,:])', bins=50)
PyPlot.xlabel("half life protein CueR (days)", fontsize=5)
PyPlot.ylabel("count", fontsize=5)

subplot(2,2,4)
PyPlot.hist((prot_half_life_array[2,:])', bins=50)
PyPlot.xlabel("half life protein Venus (days)", fontsize=5)
PyPlot.ylabel("count", fontsize=5)

PyPlot.savefig("./plots/Ensemble_plot_deg.pdf")