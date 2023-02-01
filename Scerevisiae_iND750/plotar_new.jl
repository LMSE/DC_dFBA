include("interpolation.jl")


using Plots
ENV["GKSwstype"] = "100"
using DelimitedFiles


xk_DA = readdlm("xk_DA.csv")
xk_DA=xk_DA'

vg_DA= Vector{Float64}(undef,size(xk_DA,2))
for i in 1:size(xk_DA,2)
vg_DA[i] = -vg_max*(  ((xk_DA[2,i])/(xk_DA[4,i]))/(Kg+((xk_DA[2,i])/(xk_DA[4,i]))+
            ((((xk_DA[2,i])/(xk_DA[4,i]))^2)/Kig)) )*(1/(1+(((xk_DA[3,i])/(xk_DA[4,i]))/Kie)))
end   

xk_DA[2,:]=xk_DA[2,:].*conv_glu
xk_DA[3,:]=xk_DA[3,:].*conv_eth

ts_DA = readdlm("ts_DA.csv")


xk_lab = readdlm("xk_DFBAlab.csv",',')
xk_lab=xk_lab'
xk_lab[3,:]=xk_lab[3,:].*conv_glu
xk_lab[4,:]=xk_lab[4,:].*conv_eth

ts_lab = readdlm("tsn_DFBAlab.csv")


p14 = plot(tsn_i,xk_i[3,:]./xk_i[4,:],linewidth=2,color=:green,xaxis="time [h]",yaxis="Ethanol [mM]",label="NLP-PF")


p14 = plot!(ts_lab,xk_lab[4,:]./xk_lab[1,:],linewidth=2,linestyle = :dot,color=:green,xaxis="time [h]",yaxis="Ac [mM]",label="NLP-PF",legend=:false)
p14 = plot!(ts_DA,xk_DA[3,:],linewidth=2,color=:green,linestyle = :dash,xaxis="time [h]",yaxis="Acetate [mM]",label="DA",legend=:false)


# #C_D
p14 = plot!(tsn_i,xk_i[2,:]./xk_i[4,:],linewidth=2,xaxis="time [h]",color=:red,yaxis="Conc. [g/L]",label="NLP-PF",legend=:false)

p14 = plot!(ts_lab,xk_lab[3,:]./xk_lab[1,:],linewidth=2,linestyle = :dot,color=:red,xaxis="time [h]",yaxis="Concentration [mmol/L]",label="DA",legend=:false)
p14 = plot!(ts_DA,xk_DA[2,:],linewidth=2,linestyle = :dash,color=:red,xaxis="time [h]",yaxis="Conc. [g/L]",label="DA",legend=:false,grid=:off,ylims=(0,10.5),xlims=(0.0,tsn_i[end]))

p14 = plot!(twinx(),tsn_i,xk_i[1,:]./xk_i[4,:],linewidth=2,color=:blue,legend=:false,yaxis="Biomass [g/L]",ylims=(0.2,0.8),xlims=(0.0,tsn_i[end]))
p14 = plot!(twinx(),ts_lab,xk_lab[2,:]./xk_lab[1,:],linewidth=2,color=:blue,linestyle = :dot,legend=:false,ylims=(0.2,0.8),xlims=(0.0,tsn_i[end]))
p14 = plot!(twinx(),ts_DA,xk_DA[1,:],linewidth=2,color=:blue,linestyle = :dash,xaxis="time [h]",yaxis="Biomass [g/L]",label="DA",legend=:false,ylims=(0.2,0.8),xlims=(0.0,tsn_i[end]))


p15 = plot(tsn,vk[eth,:],linewidth=2,xaxis="time [h]",yaxis="fluxes",label="Ethanol",legend=:bottomleft)
p16 = plot(tsn,vk[obj,:],linewidth=2,xaxis="time [h]",yaxis="fluxes",label="mi",legend=:bottomleft)
p17 = plot(tsn,vk[glu,:],linewidth=2,xaxis="time [h]",yaxis="fluxes",label="glu",legend=:topleft)
p17 = plot!(ts_DA,vg_DA[:],linewidth=2,xaxis="time [h]",yaxis="fluxes",label="DA",legend=:topleft)
p18 = plot(tsn,vk[o2,:],linewidth=2,xaxis="time [h]",yaxis="fluxes",label="o2",legend=:bottomleft)


p20 = plot(tsn_i,xk_i[4,:],linewidth=2,grid=:off,color=:black,xaxis="time [h]",yaxis="Volume [L]",legend=:false)
p21 = plot(tsn_i,uk_i[2,:],linewidth=2,grid=:off,color=:black,xaxis="time [h]",yaxis="F [L/min]",legend=:false)
p22 = plot(tsn_i,uk_i[1,:]*100/Osat,linewidth=2,grid=:off,color=:black,xaxis="time [h]",yaxis="DO [%]",legend=:false)

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

ljk = @layout [
    a{0.40w}  b{0.4w} _ e{0.01w}
    c{0.4w}  d{0.4w} _ f{0.01w}
]

p19=plot()

g2 = plot(p20,p14,p19,p21,p22,p19,layout=ljk,reuse = true)
g3 = plot(p15,p16,p17,p18,layout=(2,2),reuse = true)
#display(g2)

savefig(g2,"test2.pdf")

savefig(g3,"test3.pdf")