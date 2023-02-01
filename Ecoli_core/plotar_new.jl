include("interpolation.jl")


using Plots
ENV["GKSwstype"] = "100"
using DelimitedFiles


xk_DA = readdlm("xk_DA.csv")
xk_DA=xk_DA'
xk_DA_=xk_DA[:,1:69]
ts_DA = readdlm("ts_DA.csv")

ts_DA_ = ts_DA[1:69]

xk_lab = readdlm("dFBAlab_core.csv",',')
xk_lab=xk_lab'
xk_lab_=xk_lab[:,1:69]
ts_lab = readdlm("dFBAlab_core_ts.csv")
ts_lab_ = ts_lab[1:69]


p14 = plot(tsn_i,xk_i[2,:],linewidth=2,color=:green,xaxis="time [h]",yaxis="Ac [mM]",label="NLP-PF")


p14 = plot!(ts_lab,xk_lab[2,:],linewidth=2,linestyle = :dot,color=:green,xaxis="time [h]",yaxis="Ac [mM]",label="NLP-PF",legend=:false)
p14 = plot!(ts_DA,xk_DA[2,:],linewidth=2,color=:green,linestyle = :dash,xaxis="time [h]",yaxis="Acetate [mM]",label="DA",legend=:false)


# #C_D
p14 = plot!(tsn_i,xk_i[1,:],linewidth=2,xaxis="time [h]",color=:red,yaxis="Glucose [mM]",label="NLP-PF",legend=:false)


p14 = plot!(ts_lab,xk_lab[1,:],linewidth=2,linestyle = :dot,color=:red,xaxis="time [h]",yaxis="Concentration [mmol/L]",label="DA",legend=:false)
p14 = plot!(ts_DA,xk_DA[1,:],linewidth=2,linestyle = :dash,color=:red,xaxis="time [h]",yaxis="Conc. [mmol/L]",label="DA",legend=:false,grid=:off,ylims=(0,10.5))


p14 = plot!(twinx(),tsn_i,xk_i[3,:],linewidth=2,color=:blue,legend=:false,yticks = 0:0.2:0.9,ylims=(0,0.9))
p14 = plot!(twinx(),ts_lab,xk_lab[3,:],linewidth=2,color=:blue,linestyle = :dot,legend=:false,yticks = 0:0.2:0.9,ylims=(0,0.9))
p14 = plot!(twinx(),ts_DA,xk_DA[3,:],linewidth=2,color=:blue,linestyle = :dash,xaxis="time [h]",yaxis="Biomass [g/L]",label="DA",legend=:false,yticks = 0:0.2:0.9,ylims=(0,0.9))



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


savefig(g2,"test4.pdf")

