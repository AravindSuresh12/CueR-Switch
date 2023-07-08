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


tc_mods_tx_ensemble = PA[7,:]
number_of_genes = 1
number_of_samples = length(tc_mods_tx_ensemble)

# compute the time constants for TX -
tau_tx_ensemble_cuer = zeros(number_of_genes, number_of_samples)
gene_coding_length_array = custom_data_dictionary["gene_coding_length_array"]

# what is the length -
gene_length = gene_coding_length_array[1]

# what is kE?
kE = (1/gene_length)*RNAP_elongation_rate
kI = (1/characteristic_initiation_time_transcription)

for sample_index = 1:number_of_samples
    tau_factor = (kE/kI)*tc_mods_tx_ensemble[sample_index]
    tau_tx_ensemble_cuer[sample_index] = tau_factor
end

# get TX_venus time constant mods -
tc_mods_tx_ensemble = PA[8,:]
number_of_genes = 1
number_of_samples = length(tc_mods_tx_ensemble)

# compute the time constants for TX -
tau_tx_ensemble_venus = zeros(number_of_genes, number_of_samples)
gene_coding_length_array = custom_data_dictionary["gene_coding_length_array"]

# what is the length -
gene_length = gene_coding_length_array[2]

# what is kE?
kE = (1/gene_length)*RNAP_elongation_rate
kI = (1/characteristic_initiation_time_transcription)

for sample_index = 1:number_of_samples
    tau_factor = (kE/kI)*tc_mods_tx_ensemble[sample_index]
    tau_tx_ensemble_venus[sample_index] = tau_factor
end

# compute the time constants for TL cuer -
tc_mods_tl_ensemble = PA[9,:]
number_of_prots = 1
number_of_samples = length(tc_mods_tl_ensemble)
tau_tl_ensemble_cuer = zeros(number_of_prots, number_of_samples)
prot_coding_length_array = custom_data_dictionary["protein_coding_length_array"]

# what is the length -
prot_length = prot_coding_length_array[1]

# what is kE?
kE = (1/prot_length)*Ribsome_elongation_rate
kI = (1/characteristic_initiation_time_translation)

for sample_index = 1:number_of_samples
    tau_factor = (kE/kI)*tc_mods_tl_ensemble[sample_index]
    tau_tl_ensemble_cuer[sample_index] = tau_factor
end

# compute the time constants for TL venus -
tc_mods_tl_ensemble = PA[10,:]
number_of_prots = 1
number_of_samples = length(tc_mods_tl_ensemble)
tau_tl_ensemble_venus = zeros(number_of_prots, number_of_samples)
prot_coding_length_array = custom_data_dictionary["protein_coding_length_array"]

# what is the length -
prot_length = prot_coding_length_array[2]

# what is kE?
kE = (1/prot_length)*Ribsome_elongation_rate
kI = (1/characteristic_initiation_time_translation)

for sample_index = 1:number_of_samples
    tau_factor = (kE/kI)*tc_mods_tl_ensemble[sample_index]
    tau_tl_ensemble_venus[sample_index] = tau_factor
end

fig=figure() 

subplot(2,2,1)
PyPlot.hist( (tau_tx_ensemble_cuer)', bins=50)
PyPlot.xlabel("tau_factor mRNA CueR (s)", fontsize=5)
PyPlot.ylabel("count", fontsize=5)

subplot(2,2,2)
PyPlot.hist( (tau_tx_ensemble_venus)', bins=50)
PyPlot.xlabel("tau_factor mRNA Venus (s)", fontsize=5)
PyPlot.ylabel("count", fontsize=5)

subplot(2,2,3)
PyPlot.hist((tau_tl_ensemble_cuer)', bins=50)
PyPlot.xlabel("tau_factor protein CueR (s)", fontsize=5)
PyPlot.ylabel("count", fontsize=5)

subplot(2,2,4)
PyPlot.hist( (tau_tl_ensemble_venus)', bins=50)
PyPlot.xlabel("tau_factor protein Venus (s)", fontsize=5)
PyPlot.ylabel("count", fontsize=5)


PyPlot.savefig("./plots/Ensemble_plot_tc.pdf")
