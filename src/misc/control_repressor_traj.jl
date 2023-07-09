using DelimitedFiles
using Plots
using PyPlot
using CSV
using DataFrames

#S12

Test = CSV.read("Control-S12-Graph1-2.csv",DataFrame)

Time=Test[:,1] 
Operator_case1=Test[:,8] #Venus + CueR
Operator_SD_case1=Test[:,9]
Venus= Test[:,12] # Venus + GntR
Venus_SD=Test[:,13]
Operator_case2=Test[:,4] #Venus
Operator_SD_case2=Test[:,5]



PyPlot.plot(Time,Operator_case1, label="+CueR -Cu Venus")
PyPlot.errorbar(Time, Operator_case1,yerr=Operator_SD_case1,fmt="o",
mfc="#252525",mec="#252525",color="#252525", lw=1.5,ms=2, alpha=0.4)


PyPlot.plot(Time,Operator_case2, label="-CueR -Cu Venus")
PyPlot.errorbar(Time, Operator_case2,yerr=Operator_SD_case2,fmt="o",mfc="#252525",mec="#254525",color="#253325", lw=1.5,ms=2, alpha=0.4)

PyPlot.plot(Time,Venus, label="+CueR +Cu Venus")
PyPlot.errorbar(Time, Venus, yerr= Venus_SD, fmt="o",mfc="#252525",mec="#252525",color="#252325", lw=1.5,ms=2, alpha=0.4)


PyPlot.plot(figsize=(5,4))
PyPlot.xlabel("Time (hr)", fontsize=10)
PyPlot.ylabel("[Concentration] (μM)", fontsize=10)
PyPlot.tight_layout()
PyPlot.legend(fontsize=10, loc="upper left")

PyPlot.savefig("Venus_control_traj_repressor.pdf")

