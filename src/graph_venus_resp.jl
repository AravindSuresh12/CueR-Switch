include("Include.jl")

prot_array= CSV.read("./data/protein_data.csv",DataFrame)

Time=prot_array[!,"time(h)"]
Venus_only=prot_array[!,"Mean_no_CueR(uM)"]
Sterr_Venus_only=prot_array[!,"Sterr_no_CueR(uM)"]
Venus_CueR=prot_array[!,"Mean_0uM(uM)"]
Sterr_Venus_CueR=prot_array[!,"Sterr_0uM(uM)"]
Venus_CueR_Cu=prot_array[!,"Mean_100uM(uM)"]
Sterr_Venus_CueR_Cu=prot_array[!,"Sterr_100uM(uM)"]



PyPlot.plot(Time,Venus_only, label="OCueR-CueR-Cu")
PyPlot.errorbar(Time, Venus_only,yerr=Sterr_Venus_only,fmt="o",
mfc="#252525",mec="#252525",color="#252525", lw=1.5,ms=2, alpha=0.4)

PyPlot.plot(Time,Venus_CueR, label="OCueR+CueR-Cu")
PyPlot.errorbar(Time, Venus_CueR, yerr= Sterr_Venus_CueR, fmt="o",mfc="#252525",mec="#252525",color="#252525", lw=1.5,ms=2, alpha=0.4)

PyPlot.plot(Time,Venus_CueR_Cu, label="OCueR+CueR+Cu")
PyPlot.errorbar(Time, Venus_CueR_Cu,yerr=Sterr_Venus_CueR_Cu,fmt="o",mfc="#252525",mec="#254525",color="#253325", lw=1.5,ms=2, alpha=0.4)


PyPlot.plot(figsize=(5,4))


PyPlot.xlabel("Time (hr)", fontsize=10)
PyPlot.ylabel("[Concentration] (μM)", fontsize=10)
PyPlot.tight_layout()
PyPlot.legend()
PyPlot.savefig("./plots/visualize_prot_OCueR_traj.pdf")
