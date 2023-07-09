include("plot_master.jl")
#S12

Test = CSV.read("Control-S12-Graph1-2.csv",DataFrame)

Time=Test[:,1] 
Operator_case1=Test[:,2] #OCueR
Operator_SD_case1=Test[:,3] 
Venus= Test[:,4] #Venus
Venus_SD=Test[:,5]
Case2 =Test[:,6] #ocuer + cuer
Case2_SD=Test[:,7]
Venus_case2= Test[:,8] #Venus +cuer
Venus_case2_SD=Test[:,9]

PyPlot.plot(Time,Operator_case1, label="OCueR")
PyPlot.errorbar(Time, Operator_case1,yerr=Operator_SD_case1,fmt="o",mfc="#252525",mec="#252525",color="#252525", lw=1.5,ms=2, alpha=0.4)
PyPlot.plot(Time,Venus, label="Venus")
PyPlot.errorbar(Time, Venus, yerr= Venus_SD, fmt="o",mfc="#252525",mec="#252525",color="#252525", lw=1.5,ms=2, alpha=0.4)


PyPlot.plot(Time,Case2, label="OCueR + CueR")
PyPlot.errorbar(Time, Case2,yerr=Case2,fmt="o",mfc="#252525",mec="#252525",color="#252525", lw=1.5,ms=2, alpha=0.4)

PyPlot.plot(Time,Venus_case2, label="Venus+CueR")
PyPlot.errorbar(Time, Venus_case2, yerr= Venus_case2_SD, fmt="o",mfc="#252525",mec="#252525",color="#252525", lw=1.5,ms=2, alpha=0.4)


PyPlot.plot(figsize=(5,4))


PyPlot.xlabel("Time (hr)", fontsize=10)
PyPlot.ylabel("[Concentration] (μM)", fontsize=10)
PyPlot.tight_layout()
PyPlot.legend()

PyPlot.savefig("Reporter_OCueR_comp_with_CueR.pdf")

