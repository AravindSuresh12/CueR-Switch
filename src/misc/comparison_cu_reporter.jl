include("plot_master.jl")

#Graph for controls- with and without copper for reporter and T7 control

Test = CSV.read("Control-S12-Graph1-2.csv",DataFrame)

Time=Test[:,1] 
Operator_case1=Test[:,16]
Operator_SD_case1=Test[:,17]
Venus= Test[:,18]
Venus_SD=Test[:,19]
Operator_case2=Test[:,2]
Operator_SD_case2=Test[:,3]
Venus2= Test[:,4]
Venus_SD2=Test[:,5]


PyPlot.plot(Time,Operator_case1, label="5nM T7-Venus")
PyPlot.errorbar(Time, Operator_case1,yerr=Operator_SD_case1,fmt="o",
mfc="#252525",mec="#252525",color="#252525", lw=1.5,ms=2, alpha=0.4)

PyPlot.plot(Time,Venus, label="5nM T7-Venus + 100uM Cu2+")
PyPlot.errorbar(Time, Venus, yerr= Venus_SD, fmt="o",mfc="#252525",mec="#252525",color="#252525", lw=1.5,ms=2, alpha=0.4)

PyPlot.plot(Time,Operator_case2, label="5nM T7-OCueR-Venus")
PyPlot.errorbar(Time, Operator_case2,yerr=Operator_SD_case2,fmt="o",mfc="#252525",mec="#254525",color="#253325", lw=1.5,ms=2, alpha=0.4)

PyPlot.plot(Time,Venus2, label="5nM T7-Venus")
PyPlot.errorbar(Time, Venus2, yerr= Venus_SD2, fmt="o",mfc="#252525",mec="#252525",color="#252325", lw=1.5,ms=2, alpha=0.4)


PyPlot.plot(figsize=(5,4))
PyPlot.xlabel("Time (hr)", fontsize=10)
PyPlot.ylabel("[Concentration] (μM)", fontsize=10)
PyPlot.tight_layout()
PyPlot.legend()

PyPlot.savefig("comparison_cu_reporter.pdf")

