using Plots
ENV["GKSwstype"] = "100"
using DelimitedFiles

# include("interpolation.jl")
#xk = readdlm("xk.csv")
#tsn = readdlm("tsn.csv")

p12 = plot(tsn_d,xk_d[2,:],linewidth=2,color=:green,xaxis="time [h]",yaxis="Acetate [mM]",label="on-off",legend=:topleft)
p12 = plot!(tsn_pta,xk_pta[2,:],linewidth=2,color=:green,linestyle = :dot,xaxis="time [h]",yaxis="Acetate [mM]",label="pta",legend=:topleft)
p12 = plot!(tsn_wild,xk_wild[2,:],linewidth=2,color=:green,linestyle = :dash,xaxis="time [h]",yaxis="Acetate [mM]",label="wild",legend=:false)
p12 = plot!(ts,zeros(size(ts)),seriestype = :scatter,linewidth=2,color=:black)

p12 = plot!(tsn_d,xk_d[1,:],linewidth=2,xaxis="time [h]",color=:red,yaxis="Glucose [mM]",label="on-off",legend=:bottomleft)
p12 = plot!(tsn_pta,xk_pta[1,:],linewidth=2,color=:red,linestyle = :dot,xaxis="time [h]",yaxis="Glucose [mM]",label="pta",legend=:bottomleft)
p12 = plot!(tsn_wild,xk_wild[1,:],linewidth=2,color=:red,linestyle = :dash,xaxis="time [h]",yaxis="Glucose [mM]",label="wild",legend=:false)


p12 = plot!(tsn_d,xk_d[4,:],linewidth=2,color=:black,xaxis="time [h]",yaxis="Ethanol [mM]",label="on-off",legend=:topleft)
p12 = plot!(tsn_pta,xk_pta[4,:],linewidth=2,color=:black,linestyle = :dot,xaxis="time [h]",yaxis="Ethanol [mM]",label="pta",legend=:topleft)
p12 = plot!(tsn_wild,xk_wild[4,:],linewidth=2,color=:black,linestyle = :dash,xaxis="time [h]",yaxis="Concentration [mM]",label="wild",legend=:false)


p12 = plot!(twinx(),tsn_d,xk_d[3,:],linewidth=2,color=:blue,xaxis="time [h]",yaxis="Biomass [g]",label="on-off",legend=:false,yticks = 0:0.1:0.9,ylims=(0,0.3))

p12 = plot!(twinx(),tsn_pta,xk_pta[3,:],linewidth=2,color=:blue,linestyle = :dot,xaxis="time [h]",yaxis="Biomass [g]",label="pta",legend=:false,yticks = 0:0.1:0.9,ylims=(0,0.3))
p12 = plot!(twinx(),tsn_wild,xk_wild[3,:],linewidth=2,color=:blue,linestyle = :dash,xaxis="time [h]",yaxis="Biomass [g]",label="wild",legend=:false,yticks = 0:0.1:0.9,ylims=(0,0.3))


lk = @layout [
    a{0.80w} _ b{0.01w}
]

p19=plot()


g2 = plot(p12,p19,layout=lk,reuse = true)



savefig(g2,"test3.pdf")

