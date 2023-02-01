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
S = readdlm("S.csv", ',');

#iND750
eth=420 #ethanol reaction
obj=1266 #objective function index
glu=428  #glucose uptake
o2=458   #oxygen uptake

zymst=500
ergst=419
hdcea=439
ocdca=459
ocdcya=461
ocdcea=460

nc = 4
nv=size(S,2)
nm=size(S,1)
#fluxes constraints
vlb2 = readdlm("lb.csv", ',');
vlb=vlb2[:,1]
vub2 = readdlm("ub.csv", ',');
vub=vub2[:,1]

#...Model fixed Parameters
conv_eth=0.04607 #mmol/g
conv_glu=0.180156  #mmol/g

Gf = 50/conv_glu           #mmol/L
Kie= 10/conv_eth        #...mmol/L  ok
Kg = 0.5/conv_glu           #...,mmol/L  ok
Ko = 3.00e-6       #...mol/L   ok
Kig = 10/conv_glu      # mmol/L     ok
vg_max= 20.0 #mmol/gdw h    ok
vo_max= 8.0  #mmol/gdw h   ok
Osat=3.0e-4       #mol/L  ok


vlb[zymst] = -1000.0
vlb[ergst] = -1000.0
vlb[hdcea] = -1000.0
vlb[ocdca] = -1000.0
vlb[ocdcya] = -1000.0
vlb[ocdcea] = -1000.0

vlb[glu] = -vg_max
vlb[o2]= -vo_max
d=0.0*Vector{Float64}(undef,nv) 
d[obj]=-1
up=0.0*Vector{Float64}(undef,nv) 
up[glu] = 1

#initial conditions
x0=0.2         #biomas initial condition g
g0=2.0/conv_glu       #glucose initial condition mmol    
e0=0.0/conv_eth        #ethanol initial condition mmol
v0=1.0         #L

c0 = [x0,g0,e0,v0]

v0=zeros(nv)

#scaling
cs=Vector{Float64}(undef,nc) 

for i in 1:nc
cs[i]=1.0
end

vs=Vector{Float64}(undef,nv) 
for i in 1:nv

    if vlb[i] == 0.0 && vub[i] == 0.0

         vs[i] = 1.0
    else

        vs[i] = 1.0
      
    end

end



#integration parameters
nfe    = 4   # number of control intervals
ncp    = 3           # number of collocation points
th     = 3.0        # time horizon
h = th/nfe           # length of one finite element on the time horizon
hm=h*ones(nfe)'
var_h=1.0
v_dot = 0.1

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

#control variables
ncv=2
uk = Matrix{Float64}(undef,ncv,nfe)
uk[1,:].=Osat/2
uk[2,:].= 0.0

uk[1,3:4].=0.0
uk[2,3:4].= 0.01

    #tick()
       
        sol, solv, solcd, soll, solal, solag, solh = pFBA_KKT_flux(c0,v0)

        for i in 2:nfe+1
          
            ts[i] = ts[i-1] + solh[i-1]
          
        end
        
        #tock()
   
#building vectors
xk[:,1]=c0
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

