using DelimitedFiles
using Plots
using PyPlot
using CSV
using DataFrames

#S12

Test = CSV.read("./data/mRNA_data.csv",DataFrame)

Time=Test[:,1] 
mRNA_base=Test[:,2]
STDERR_base=Test[:,5]
mRNA_case1=Test[:,3]
STDERR_case1=Test[:,6]
mRNA_case3=Test[:,4]
STDERR_case3=Test[:,7]



PyPlot.plot(Time,mRNA_base, label="-CueR -Cu")
PyPlot.errorbar(Time, mRNA_base , yerr=STDERR_base,fmt="o",
mfc="#252525",mec="#252525",color="#252525", lw=1.5,ms=2, alpha=0.4)

PyPlot.plot(Time,mRNA_case1, label="+CueR -Cu")
PyPlot.errorbar(Time, mRNA_case1 , yerr=STDERR_case1,fmt="o",
mfc="#252525",mec="#252525",color="#252525", lw=1.5,ms=2, alpha=0.4)

PyPlot.plot(Time,mRNA_case3, label="+CueR +Cu")
PyPlot.errorbar(Time, mRNA_case3 , yerr=STDERR_case3,fmt="o",
mfc="#252525",mec="#252525",color="#252525", lw=1.5,ms=2, alpha=0.4)

PyPlot.plot(figsize=(5,4))
PyPlot.xticks([0,2,4,8,16], fontsize=10)
PyPlot.xlabel("Time (hr)", fontsize=10)
PyPlot.ylabel("[Concentration] (nM)", fontsize=10)
PyPlot.tight_layout()
PyPlot.legend(fontsize=10, loc="upper right")
PyPlot.savefig("./plots/mRNA-traj.pdf")

