include("interpolation.jl")


using Plots
ENV["GKSwstype"] = "100"
using DelimitedFiles


xk_DA = readdlm("xk_dFBA_iJ.csv")
ts_DA = readdlm("tsn_dFBA_iJ.csv")


xk_lab = readdlm("dFBAlab_iJ.csv",',')
xk_lab=xk_lab'
ts_lab = readdlm("dFBAlab_iJ_ts.csv")


p14 = plot(tsn_i,xk_i[2,:],linewidth=2,color=:green,xaxis="time [h]",yaxis="Ac [mM]",label="NLP-PF")

p14 = plot!(ts,zeros(size(ts)),seriestype = :scatter,linewidth=2,color=:black,xaxis="time [h]",yaxis="Ac [mM]",label="NLP-PF")

p14 = plot!(ts_lab,xk_lab[2,:],linewidth=2,linestyle = :dot,color=:green,xaxis="time [h]",yaxis="Ac [mM]",label="NLP-PF",legend=:false)
p14 = plot!(ts_DA,xk_DA[2,:],linewidth=2,color=:green,linestyle = :dash,xaxis="time [h]",yaxis="Acetate [mM]",label="DA",legend=:false)


# #C_D
p14 = plot!(tsn_i,xk_i[1,:],linewidth=2,xaxis="time [h]",color=:red,yaxis="Glucose [mM]",label="NLP-PF",legend=:false)


p14 = plot!(ts_lab,xk_lab[1,:],linewidth=2,linestyle = :dot,color=:red,xaxis="time [h]",yaxis="Concentration [mmol/L]",label="DA",legend=:false)
p14 = plot!(ts_DA,xk_DA[1,:],linewidth=2,linestyle = :dash,color=:red,xaxis="time [h]",yaxis="Concentration [mmol/L]",label="DA",legend=:false,grid=:off,ylims=(0,10.5))


p14 = plot!(twinx(),tsn_i,xk_i[3,:],linewidth=2,color=:blue,legend=:false,yticks = 0:0.2:0.9,ylims=(0,0.95))
p14 = plot!(twinx(),ts_lab,xk_lab[3,:],linewidth=2,color=:blue,linestyle = :dot,legend=:false,yticks = 0:0.2:0.9,ylims=(0,0.95))
p14 = plot!(twinx(),ts_DA,xk_DA[3,:],linewidth=2,color=:blue,linestyle = :dash,xaxis="time [h]",yaxis="Biomass [g/L]",label="DA",legend=:false,yticks = 0:0.2:0.95,ylims=(0,0.95))


l = @layout [
    a    
    b  c 
]

lj = @layout [
    a{0.80w} _ b{0.01w}
    c{0.80w} _ d{0.01w}
]

lk = @layout [
    a{0.80w} _ b{0.01w}
]

p19=plot()

g2 = plot(p14,p19,layout=lk,reuse = true)


savefig(g2,"test2.pdf")

