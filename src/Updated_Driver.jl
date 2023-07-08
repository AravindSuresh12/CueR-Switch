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

PC = readdlm("./simulated/POETS/PC_T5.dat") #N=20
#PC= readdlm("./simulated/POETS_N_Vary/N=5/PC_T5.dat") 
#PC = readdlm("./simulated/POETS_N_Vary/N=100/PC_T5.dat") 
average_params = mean(PC,dims=2)
poets_params=average_params

# Update the data dictionary
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

## Step1- Plot Protein simulated data- Venus and also plot experimental data

P1 = Plots.plot(T,X[:,6], label = "Venus_Model_Prot",xlabel="Time (hr)",ylabel = "uM",linewidth=3)

# # plot experimental data- protein Venus
prot_data = CSV.read("./data/protein_data.csv",DataFrame)
Tx= prot_data[!,"time(h)"]
mean_Venus = prot_data[!,"Mean_100uM(uM)"]
stdev_Venus = prot_data[!,"Sterr_100uM(uM)"]

P1= Plots.scatter!(Tx, mean_Venus, label="Venus_Expt_Protein",xlabel="Time (hr)",ylabel = "uM",linewidth=3, legendfontsize=4)

#Step 2- Plot mRNA data of CueR and Venus- experimental and simulated and compare

mRNA_data = CSV.read("./data/mRNA_data.csv",DataFrame)
T1 = mRNA_data[!,"Average_time(h)"]
mRNA_Venus = mRNA_data[!,"<Venus+CueR+Cu>(nM)"]
stdev_Venus = mRNA_data[!,"SE3(nM)"]
mRNA_CueR = mRNA_data[!,"<CueR>(nM)"]
stdev_CueR = mRNA_data[!,"SE_1(nM)"]

P2= Plots.scatter(T1, mRNA_CueR, label="CueR_mRNA",xlabel="Time (hr)",ylabel = "nM",linewidth=3, yerror=stdev_CueR)
P2= Plots.plot!(T, 1000*(X[:,3]), label="CueR_Simulated_mRNA",xlabel="Time (hr)",ylabel = "nM",linewidth=3, legend=:bottomright, legendfontsize=4)


P3= Plots.scatter(T1, mRNA_Venus, label="Venus_mRNA",xlabel="Time (hr)",ylabel = "nM",linewidth=3, yerror=stdev_Venus)
P3= Plots.plot!(T, 1000*(X[:,4]), label="Venus_Simulated_mRNA",xlabel="Time (hr)",ylabel = "nM",linewidth=3, legend=:bottomright, legendfontsize=4)

#step3- plot dose response experimental versus model simulated

dose_data = CSV.read("./data/dose_response.csv",DataFrame)
exp_val= dose_data[!,"Venus(uM)"]
exp_std= dose_data[!,"STERR(uM)"]

empty_test=[]

copper_sim=[0,5,10,50,100]

for element in copper_sim

    local cuprous_parameter_dictionary["copper_salt_conc"]=element
    local (T,X) = SolveBalances(time_start,time_stop,time_step_size,data_dictionary)
    local end_value=last(X[:,6])
    push!(empty_test,end_value)
end


P4= Plots.plot(copper_sim,empty_test, xlabel="copper salt concentration (μM)", ylabel="Venus (μM)", label="simulated dose response")

P4= Plots.scatter!(copper_sim,exp_val, xlabel="copper salt concentration (μM)", ylabel="Venus (μM)",label="experimental dose response", legend=:bottomright, legendfontsize=4)

Plots.plot(P1,P2,P3,P4, layout=(2,2))

#Plots.savefig("./plots/Ensemble.pdf")
#Plots.savefig("./plots/Model_fit_comp_N=5.pdf")
Plots.savefig("./plots/Model_fit_comp_N=100.pdf")