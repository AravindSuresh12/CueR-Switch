include("plot_master.jl")

#S12

Test = CSV.read("NOR-CTRL1.csv",DataFrame)

Time=Test[:,1] 
Operator_case1=Test[:,2]
Operator_SD_case1=Test[:,3]
Operator_case2=Test[:,4]
Operator_SD_case2=Test[:,5]
Operator_case3=Test[:,6]
Operator_SD_case3=Test[:,7]
Operator_case4=Test[:,8]
Operator_SD_case4=Test[:,9]



PyPlot.plot(Time,Operator_case1, label="OCueR-EcRNAP")
PyPlot.errorbar(Time, Operator_case1,yerr=Operator_SD_case1,fmt="o",
mfc="#252525",mec="#252525",color="#252525", lw=1.5,ms=2, alpha=0.4)


PyPlot.plot(Time,Operator_case2, label="OCueR+EcRNAP")
PyPlot.errorbar(Time, Operator_case2,yerr=Operator_SD_case2,fmt="o",mfc="#252525",mec="#254525",color="#253325", lw=1.5,ms=2, alpha=0.4)

PyPlot.plot(Time,Operator_case3, label="Venus-EcRNAP")
PyPlot.errorbar(Time, Operator_case3,yerr=Operator_SD_case3,fmt="o",mfc="#252525",mec="#254525",color="#253325", lw=1.5,ms=2, alpha=0.4)


PyPlot.plot(Time,Operator_case4, label="Venus+EcRNAP")
PyPlot.errorbar(Time, Operator_case4,yerr=Operator_SD_case3,fmt="o",mfc="#252525",mec="#254525",color="#253325", lw=1.5,ms=2, alpha=0.4)



PyPlot.plot(figsize=(5,4))
PyPlot.xlabel("Time (hr)", fontsize=10)
PyPlot.ylabel("Normalized RFU units wrt T7-OCueR-Venus", fontsize=10)
PyPlot.tight_layout()
PyPlot.legend(fontsize=10, loc="upper left")

PyPlot.savefig("Control_PURE_Venus_comparison.pdf")

