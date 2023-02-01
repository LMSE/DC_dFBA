#=
      simulation of dFBA by orthogonal collocation
=#

#−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−#
#Formulating the discrete dynamics
using TickTock
using JuMP
using Ipopt
using LinearAlgebra
#reading GSM model
using DelimitedFiles

include("pFBA_KKT_flux.jl")

#stoichiometric matrix
S = readdlm("Si.csv", ',');

#E.coli iJR904
obj=269  #objective function index
glu=398  #glucose uptake
o2=500  #oxygen uptake
ac=289
ATP=212

nc = 3
nv=size(S,2)
nm=size(S,1)
#fluxes constraints
vlb2 = readdlm("lbi.csv", ',');
vlb=vlb2[:,1]
vub2 = readdlm("ubi.csv", ',');
vub=vub2[:,1]

#...Model fixed Parameters
Kg = 0.015         #
Ko=0.01
vg_max= 10.5
vlb[o2]= -19.0
vlb[ac] = -2.5

vlb[glu] = -vg_max
d=0.0*Vector{Float64}(undef,nv) 
d[obj]=-1
up=0.0*Vector{Float64}(undef,nv) 
up[glu] = 1

#initial conditions
b0=0.003        #                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
g0=10.4        #
ac0=0.3         #
eth0=0.0
c0=[g0,ac0,b0]
v0=zeros(nv)


#scaling
cs=Vector{Float64}(undef,nc) 
cs[1]=1.0
cs[2]=1.0
cs[3]=1.0

vs=Vector{Float64}(undef,nv) 
for i in 1:nv

    if vlb[i] == 0.0 && vub[i] == 0.0

         vs[i] = 1.0
    else

        vs[i] = 1.0
      
    end

end



#integration parameters
nfe    = 8   # number of control intervals
ncp    = 3           # number of collocation points
th     = 7.1        # time horizon
h = th/nfe           # length of one finite element on the time horizon
hm=h*ones(nfe)'
var_h=1.0

w=1e-20        #weigth for flux sum minimization on the OF
phi1=1e0
phi2=1e0
phi3=1e0

ts     = Vector{Float64}(undef,nfe) # time series for plotting
tsn     = Vector{Float64}(undef,nfe*ncp+1)
xk = Matrix{Float64}(undef,nc,nfe*ncp+1)
vk = Matrix{Float64}(undef,nv,nfe*ncp+1)
MPEC = Matrix{Float64}(undef,nv,nfe*ncp+1)
MPEC_vg = Vector{Float64}(undef,nfe*ncp+1)

ts = pushfirst!(ts,0)


    

    #tick()
       
        sol, solv, solcd, soll, solal, solag, solh = pFBA_KKT_flux(c0,v0)

        for i in 2:nfe+1
            
            ts[i] = ts[i-1] + solh[i-1]
          
        end
        
        #tock()
   
#building vectors
xk[:,1]=c0
vk[:,1]=v0
MPEC[:,1]=v0
tsn[1] = 0.0
global kk = 2

radau  = [0.15505 0.64495 1.00000]

for i in 1:nfe


   for j in 1:ncp
    xk[:,kk]=sol[:,i,j]
     tsn[kk] = ts[i] + radau[j]*solh[i]

    global kk=kk+1
   end


end

writedlm("xk.csv",xk)
writedlm("tsn.csv",tsn)

