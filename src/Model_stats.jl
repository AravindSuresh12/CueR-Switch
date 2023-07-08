#Model stats

include("AUC_model.jl")


#path_to_ensemble_file = "./simulated/POETS/PC_T5.dat" 
#path_to_ensemble_file = "./simulated/POETS_N_Vary/N=5/PC_T5.dat" 
path_to_ensemble_file = "./simulated/POETS_N_Vary/N=100/PC_T5.dat"
#Uncomment us for calculating model performance for increasing N


poets= readdlm(path_to_ensemble_file);
average_params = mean(PC,dims=2);
parameter_guess_array=average_params;

dose_data = CSV.read("./data/dose_response.csv",DataFrame)
exp_val= dose_data[!,"Venus(uM)"]
exp_std= dose_data[!,"STERR(uM)"]

empty_test=[]

copper_sim=[0,5,10,50,100]

for element in copper_sim
    push!(empty_test,predict(path_to_ensemble_file, element))
end


##MAE Parameter 

mae = mean(abs.(empty_test .- exp_val))
target_range = maximum(exp_val) - minimum(exp_val)
mae_percentage = (mae / target_range) * 100.0


##RMSE Parameter
squared_diff = (empty_test .- exp_val).^2
mean_squared_diff = mean(squared_diff)
rmse = sqrt(mean_squared_diff)

print("The Mean Absolute error percentage is $(mae_percentage) and The RMSE value is $(rmse)")

