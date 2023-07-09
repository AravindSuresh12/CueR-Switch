
include("plot_master.jl")

#graph for appendix images- 

Test = CSV.read("Venus-ECRNAP.csv",DataFrame)

Time=Test[:,1] 
Operator_case_1=Test[:,2]
Operator_case1_SD=Test[:,3]
Operator_case_2=Test[:,4]
Operator_case2_SD=Test[:,5]
Operator_case_3=Test[:,6]
Operator_case3_SD=Test[:,7]
Operator_case_4=Test[:,8]
Operator_case4_SD=Test[:,9]
Operator_case_5=Test[:,10]
Operator_case5_SD=Test[:,11]
Operator_case_6=Test[:,12]
Operator_case6_SD=Test[:,13]




PyPlot.plot(Time,Operator_case_1, label="-CueR -Cu -EcRNAP")
PyPlot.errorbar(Time, Operator_case_1,yerr=Operator_case1_SD,fmt="o",
mfc="#252525",mec="#252525",color="#252525", lw=1.5,ms=2, alpha=0.4)


PyPlot.plot(Time,Operator_case_2, label="+CueR -Cu -EcRNAP")
PyPlot.errorbar(Time, Operator_case_2,yerr=Operator_case2_SD,fmt="o",mfc="#252525",mec="#254525",color="#253325", lw=1.5,ms=2, alpha=0.4)

PyPlot.plot(Time,Operator_case_3, label="+CueR +Cu -EcRNAP")
PyPlot.errorbar(Time, Operator_case_3,yerr=Operator_case3_SD,fmt="o",mfc="#252525",mec="#254525",color="#253325", lw=1.5,ms=2, alpha=0.4)

PyPlot.plot(Time,Operator_case_4, label="-CueR -Cu +EcRNAP")
PyPlot.errorbar(Time, Operator_case_4,yerr=Operator_case4_SD,fmt="o",
mfc="#252525",mec="#252525",color="#252525", lw=1.5,ms=2, alpha=0.4)

PyPlot.plot(Time,Operator_case_5, label="+CueR -Cu -EcRNAP")
PyPlot.errorbar(Time, Operator_case_5,yerr=Operator_case5_SD,fmt="o",mfc="#252525",mec="#254525",color="#253325", lw=1.5,ms=2, alpha=0.4)


PyPlot.plot(Time,Operator_case_6, label="+CueR +Cu +EcRNAP")
PyPlot.errorbar(Time, Operator_case_6,yerr=Operator_case6_SD,fmt="o",mfc="#252525",mec="#254525",color="#253325", lw=1.5,ms=2, alpha=0.4)

PyPlot.plot(figsize=(5,4))
PyPlot.xlabel("Time (hr)", fontsize=10)
PyPlot.ylabel("Normalized RFU wrt T7-OCueR-Venus -EcRNAP", fontsize=10)
PyPlot.tight_layout()
PyPlot.legend(fontsize=6, loc="upper left")

PyPlot.savefig("EcRNAP_comp_full_switch_PURE.pdf")
