using DelimitedFiles
using Plots
using PyPlot
using CSV
using DataFrames

#S12

Test = CSV.read("Dose-Plasmid-Test.csv",DataFrame)


Time=Test[:,1] 
Operator_case1=Test[:,2] #OCueR
Operator_SD_case1=Test[:,3] 
Case2= Test[:,4] #Venus
Case2_SD=Test[:,5]
Case3 =Test[:,6] #ocuer + cuer
Case3_SD=Test[:,7]


PyPlot.plot(Time,Operator_case1, label="OCueR 2.5nM")
PyPlot.errorbar(Time, Operator_case1,yerr=Operator_SD_case1,fmt="o",
mfc="#252525",mec="#252525",color="#252525", lw=1.5,ms=2, alpha=0.4)


PyPlot.plot(Time,Case2, label="OCueR 5nM")
PyPlot.errorbar(Time, Case2,yerr=Case2_SD,fmt="o",mfc="#252525",mec="#254525",color="#253325", lw=1.5,ms=2, alpha=0.4)

PyPlot.plot(Time,Case3, label="OCueR 10nM")
PyPlot.errorbar(Time, Case3, yerr= Case3, fmt="o",mfc="#252525",mec="#252525",color="#252325", lw=1.5,ms=2, alpha=0.4)


PyPlot.plot(figsize=(5,4))
PyPlot.xlabel("Time (hr)", fontsize=10)
PyPlot.ylabel("[Concentration] (μM)", fontsize=10)
PyPlot.tight_layout()
PyPlot.legend(fontsize=10, loc="upper left")

PyPlot.savefig("dose_plasmid_var_conc.pdf")