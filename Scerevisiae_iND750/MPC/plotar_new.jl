
using Plots
ENV["GKSwstype"] = "100"
using DelimitedFiles

xk_lab = readdlm("xk_dFBAlab.csv",',')
xk_lab=xk_lab'
xk_lab[3,:]=xk_lab[3,:].*conv_glu
xk_lab[4,:]=xk_lab[4,:].*conv_eth

uk_lab = readdlm("uk_dFBAlab.csv",',')

ts_lab = readdlm("tsn_dFBAlab.csv")



include("interpolation.jl")



p14 = plot(tsn_i,xk_i[3,:]./xk_i[4,:],linewidth=2,color=:green,xaxis="time [h]",yaxis="Ethanol [mM]",label="NLP-PF")


p14 = plot!(ts_lab,xk_lab[4,:]./xk_lab[1,:],linewidth=2,linestyle = :dot,color=:green,xaxis="time [h]",yaxis="Ac [mM]",label="NLP-PF",legend=:false)



# #C_D
p14 = plot!(tsn_i,xk_i[2,:]./xk_i[4,:],linewidth=2,xaxis="time [h]",color=:red,yaxis="Conc. [g/L]",label="NLP-PF",legend=:false)

p14 = plot!(ts_lab,xk_lab[3,:]./xk_lab[1,:],linewidth=2,linestyle = :dot,color=:red,xaxis="time [h]",yaxis="Conc. [g/L]",label="DA",legend=:false)

p14 = plot!(twinx(),tsn_i,xk_i[1,:]./xk_i[4,:],linewidth=2,grid=:off,color=:blue,legend=:false,yaxis="Biomass [g/L]",ylims=(0.2,42.0))
p14 = plot!(twinx(),ts_lab,xk_lab[2,:]./xk_lab[1,:],linewidth=2,color=:blue,linestyle = :dot,legend=:false,ylims=(0.2,42.0))


p20 = plot(tsn_i,xk_i[4,:],linewidth=2,grid=:off,color=:black,xaxis="time [h]",yaxis="Volume [L]",legend=:false)
Vmax = [1.2] .*ones(1,nfe+1)
p20 = plot!(ts_lab,xk_lab[1,:],linewidth=2,linestyle = :dot,grid=:off,color=:black,xaxis="time [h]",yaxis="Volume [L]",legend=:false)
p20 = plot!(ts,Vmax[1,:],linewidth=5,linestyle = :dot,linecolor = :red)

p21 = plot(tsn_i,uk_i[2,:],linewidth=2,grid=:off,color=:black,xaxis="time [h]",yaxis="F [L/min]",legend=:false)
p21 = plot!(ts_lab_ii,uk_lab_i[1,:],linewidth=2,linestyle = :dot,grid=:off,color=:black,xaxis="time [h]",yaxis="F [L/min]",legend=:false)


p22 = plot(tsn_i,uk_i[1,:]*100/Osat,linewidth=2,grid=:off,color=:black,xaxis="time [h]",yaxis="DO [%]",legend=:false)
p22 = plot!(ts_lab_ii,uk_lab_i[2,:]*100,linewidth=2,linestyle = :dot,grid=:off,color=:black,xaxis="time [h]",yaxis="DO [%]",legend=:false)

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
g3 = plot(p14,p19,layout=lk,reuse = true)


savefig(g2,"test4.pdf")

savefig(g3,"test3.pdf")