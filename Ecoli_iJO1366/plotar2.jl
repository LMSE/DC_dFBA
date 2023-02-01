#Ra32353235
#  daviddeo@sobol.biozone.utoronto.ca:9833
#  scp daviddeo@sobol.biozone.utoronto.ca:/nfs/homes/daviddeo/ts.csv ./
#  cd("Documents/UT/dFBA/code/ecoli/")

using Plots
ENV["GKSwstype"] = "100"
using DelimitedFiles

#include("interpolation.jl")
#xk = readdlm("xk.csv")
#tsn = readdlm("tsn.csv")

#xk_DA = readdlm("xk_DA.csv")
#xk_DA=xk_DA'
#xk_DA = readdlm("dFBAlab_core.csv",',')
#xk_DA=xk_DA'
#xk_DA=xk_DA[:,1:69]
#ts_DA = readdlm("ts_DA.csv")
#ts_DA = readdlm("dFBAlab_core_ts.csv")
#ts_DA = ts_DA[1:69]
#vk = readdlm("vk.csv")
#lStar = readdlm("lStar.csv")
#auStar = readdlm("auStar.csv")
#alStar = readdlm("alStar.csv")
#agStar = readdlm("agStar.csv")



p12 = plot(tsn_i,xk_i[1,:],linewidth=2,xaxis="time [h]",color=:red,yaxis="Glucose [mM]",legend=:false)
p12 = plot!(tsn_i_red,xk_i_red[1,:],linewidth=2,linestyle = :dot,xaxis="time [h]",color=:red,yaxis="Glucose [mM]",legend=:false)

p12 = plot!(tsn_i,xk_i[4,:],linewidth=2,color=:green,grid=:off,xaxis="time [h]",yaxis="Concentration [mmol/L]",legend=:false)
p12 = plot!(tsn_i_red,xk_i_red[4,:],linewidth=2,linestyle = :dot,color=:green,grid=:off,xaxis="time [h]",yaxis="Concentration [mmol/L]",legend=:false)

p12 = plot!(twinx(),tsn_i,xk_i[3,:],linewidth=2,color=:blue,grid=:off,xaxis="time [h]",yaxis="Biomass [g/L]",legend=:false,ylims=(0,26.0))
p12 = plot!(twinx(),tsn_i_red,xk_i_red[3,:],linewidth=2,linestyle = :dot,color=:blue,grid=:off,xaxis="time [h]",yaxis="Biomass [g/L]",legend=:false,ylims=(0,26.0))

####

lk = @layout [
    a{0.80w} _ b{0.01w}
]

p19=plot()

g2 = plot(p12,p19,layout=lk,reuse = true)


savefig(g2,"test3.pdf")

