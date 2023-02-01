#=
DC_dFBA
=#

#−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−#
#Formulating the discrete dynamics
using JuMP
using Ipopt
using LinearAlgebra
using TickTock
using Plots
#reading GSM model
using FileIO, JLD2
data = FileIO.load("data.jld2","data");
using DelimitedFiles

#stoichiometric matrix
S = readdlm("S.csv", ',');

#Yeast 8.3
nm=2666   #metabolites number
nr=3928   #reactions number
eth=2630  #ethanol reaction
obj=3414  #objective function index
glu=2588  #glucose uptake
o2=2816   #oxygen uptake
ATP = 3415
xyl=2592;  #xylose uptake

nc = 4
nv=size(S,2)
nm=size(S,1)
#fluxes constraints
vlb2 = readdlm("lb.csv", ',');
vlb=vlb2[:,1]
vub2 = readdlm("ub.csv", ',');
vub=vub2[:,1]

vlb[o2]= 0.0
vub[o2]= 0.0

vlb[ATP]= 0.0

#estimation
teta_opt=[7.3,1.03,32.0,14.85,0.5]
#teta0=teta_opt
teta0=[7.5,1.0,35.0,15.00,0.5]
teta0=log.(teta0)
UB=[8.0,1.13,33.0,15.85,0.6]
#UB=teta_opt
#UB=[7.4,2.0,40.0,20.0,1.0]
UB=log.(UB)

LB=[7.0,0.9,31.0,13.85,0.4]
#LB=teta_opt
#LB=[7.3,0.5,30.0,10.0,0.1]
LB=log.(LB)

#...Model fixed Parameters
#vg_max=7.3;
#Kg=1.026;
#vz_max=32.0;
#Kz=14.85;
#Kig=0.5;

#initial conditions
x0=0.2         #biomas initxkial condition g
g0=4.0        #glucose initial condition g
z0=2.0         #xylose initial condition g
e0=0.0         #ethanol initial condition g


#integration parameters
nc = 4
nfe    = 12        # number of control intervals
ncp    = 3           # number of collocation points
th     = 22          # time horizon
h = th/nfe           # length of one finite element on the time horizon
ph = nfe          #prediction horizon
np=5                #number of parameter to be estimate

hm=h*ones(nfe)'
var_h=1.0

w=1e-20        #weigth for flux sum minimization on the OF
omega = 1e2
phi1=1e0
phi2=1e0
phi3=1e0

d=0.0*Vector{Float64}(undef,nv) 
d[obj]=-1
up=0.0*Vector{Float64}(undef,nv) 
up[glu] = 1
up2=0.0*Vector{Float64}(undef,nv) 
up2[xyl] = 1
n_up=2

ts     = Vector{Float64}(undef,nfe) # time series for plotting
tsn     = Vector{Float64}(undef,nfe*ncp+1)
xk = Matrix{Float64}(undef,nc,nfe*ncp+1)
vk = Matrix{Float64}(undef,nv,nfe*ncp+1)

v0=zeros(nv)

#scaling
cs=Vector{Float64}(undef,nc) 

cs[1]=1.0
cs[2]=1.0
cs[3]=1.0
cs[4]=1.0

vs=Vector{Float64}(undef,nv) 
for i in 1:nv

    if vlb[i] == 0.0 && vub[i] == 0.0

         vs[i] = 1.0
    else

        vs[i] = 1.0
      
    end

end

for i in 1:nfe
    ts[i] = h*i
end
ts = pushfirst!(ts,0)


    c0 = [x0,g0,z0,e0]

        colmat = [0.19681547722366  -0.06553542585020 0.02377097434822;
                  0.39442431473909  0.29207341166523 -0.04154875212600;
                  0.37640306270047  0.51248582618842 0.11111111111111]
        radau  = [0.15505 0.64495 1.00000]

        #−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−#
        # JuMP model
        
       m = Model(with_optimizer(Ipopt.Optimizer, warm_start_init_point = "yes", print_level = 5, linear_solver = "ma27",tol=1e-4,acceptable_iter=5,acceptable_tol=1e-2) )

  

        # Set up variables for JuMP
        @variables(m, begin
            c[1:nc, 1:ph, 1:ncp]     #States: [mol/L] concentration (states)
            cdot[1:nc, 1:ph, 1:ncp]  #Differential equations: [mol/L/min]
            FO                       #Objective Function
            teta[1:np]                #parameters
            
            v[1:nv, 1:nfe] #fluxes                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
            lambda[1:nm, 1:nfe]  #multipliers
            alpha_U[1:nv, 1:nfe]
            alpha_L[1:nv, 1:nfe]
            alpha_upt[1:n_up,1:nfe]

            FO_U[1:nv, 1:nfe]
            FO_L[1:nv, 1:nfe]
            FO_upt[1:n_up,1:nfe]

            hv[1:nfe]
        end)

        # Set up initial guesses for solver
        for i in 1:ph
            for j in 1:ncp
                set_start_value(c[1,i,j], c0[1])      #v0
                set_start_value(c[2,i,j], c0[2])      #x0
                set_start_value(c[3,i,j], c0[3])      #g0
                set_start_value(c[4,i,j], c0[4])      #e0
               # set_start_value(FO, 1e4)
              #  set_start_value(teta[1], teta0[1])
              #  set_start_value(teta[2], teta0[2])
              #  set_start_value(teta[3], teta0[3])
              #  set_start_value(teta[4], teta0[4])
              #  set_start_value(teta[5], teta0[5])
            end

            #for mc in 1:nu
            #    set_start_value(u[mc,i], u0[mc])
                #set_start_value(du[mc,i], 0)
            #end
        end
        
        for i in 1:nfe
            set_start_value(hv[i], hm[i])
        
          end
        
           #scaling
         for i in 1:nc
          c0[i]=c0[i]/cs[i]
        
        end

        # Set up objective function
        @NLobjective(m, Min, omega*FO + sum( sum(-phi1*FO_L[mc,i] + -phi3*FO_U[mc,i] for mc in 1:nv) + phi2*FO_upt[1,i] + phi2*FO_upt[2,i] for i in 1:nfe))

            

        @NLexpressions(m, begin
        vg[i=1:ph], exp(teta[1])*(  (c[2,i,3])/(exp(teta[2])+(c[2,i,3])
        ) )
        vz[i=1:ph], exp(teta[3])*(  (c[3,i,3])/(exp(teta[4])+(c[3,i,3])
        ))*(1/ (1 + (c[2,i,3]/exp(teta[5]))  )  )
    
        end)

        #Set up the constraints
        @constraints(m, begin
            
            
            # set up collocation equations - 2nd-to-nth point
            coll_c_n[l=1:nc, i=2:ph, j=1:ncp], c[l,i,j] == c[l,i-1,ncp]+h*sum(colmat[j,k]*cdot[l,i,k] for k in 1:ncp)
            # set up collocation equations - 1st point
            coll_c_0[l=1:nc, i=1, j=1:ncp], c[l,i,j] == c0[l] + h*sum(colmat[j,k]*cdot[l,i,k] for k in 1:ncp)
            #parameter constraints
            teta_LB[p=1:np], teta[p] >= LB[p]
            teta_UB[p=1:np], teta[p] <= UB[p]

            Sc[mc=1:nm,i=1:nfe],  sum(S[mc,k]*v[k,i]*vs[k] for k in 1:nv) == 0
            v_UB[mc=1:nv,i=1:nfe], v[mc,i]*vs[mc] - vub[mc]  <= 0
            v_LB[mc=1:nv,i=1:nfe], -v[mc,i]*vs[mc]  + vlb[mc] <= 0

            c_LB[mc=1:nc,i=1:nfe, j=1:ncp], -c[mc,i,j]  <= 0

            MFE1, sum(hv[i] for i in 1:nfe) == th
            MFE3[i=1:nfe], hv[i]  >= 0.0
            MFE4[i=1:nfe], hv[i]  >= (1-var_h)*hm[1]
            MFE5[i=1:nfe], hv[i]  <= (1+var_h)*hm[1]

            #pFBA
            Lagr[mc=1:nv,i=1:nfe], +d[mc] + w*v[mc,i]*vs[mc]  + alpha_L[mc,i] + alpha_U[mc,i] + up[mc]*alpha_upt[1,i] + up2[mc]*alpha_upt[2,i]  + sum(S[k,mc]*lambda[k,i] for k in 1:nm) == 0
            
            alpha1_LB[mc=1:nv,i=1:nfe], alpha_L[mc,i] <= 0
            alpha4_LB[mc=1:n_up,i=1:nfe], alpha_upt[mc,i] <= 0
            alpha1_UB[mc=1:nv,i=1:nfe], alpha_U[mc,i] >= 0   


            m8, FO == sum( sum( sum( (data[i,j,mc]-c[i,j,mc])^2 for i in 1:nc)   for j in 1:ph)   for mc in 1:ncp)
        
        end)
        
        @NLconstraints(m, begin
          # set up differential equations
          m1[i=1:ph, j=1:ncp], cdot[1,i,j] == v[obj,i]*c[1,i,j]
          m2[i=1:ph, j=1:ncp], cdot[2,i,j] == - 0.180156*vg[i]*c[1,i,j]
          m3[i=1:ph, j=1:ncp], cdot[3,i,j] == - 0.15013*vz[i]*c[1,i,j]
          m4[i=1:ph, j=1:ncp], cdot[4,i,j] ==  0.04607*v[eth,i]*c[1,i,j] 

         v_LB_g[i=1:nfe], -v[glu,i]*vs[glu]  - vg[i] <= 0  
         v_LB_z[i=1:nfe], -v[xyl,i]*vs[xyl]  - vz[i] <= 0 
    
         FO1[mc=1:nv,i=1:nfe], FO_L[mc,i] == (v[mc,i]*vs[mc] -vlb[mc])*alpha_L[mc,i]
         FO2[mc=1:nv,i=1:nfe], FO_U[mc,i] == (v[mc,i]*vs[mc] -vub[mc])*alpha_U[mc,i]
         FO3_upt[i=1:nfe], FO_upt[1,i] == (-v[glu,i]*vs[glu] -vg[i])*alpha_upt[1,i]
         FO4_upt[i=1:nfe], FO_upt[2,i] == (-v[xyl,i]*vs[xyl] -vz[i])*alpha_upt[2,i]

        end)

        #−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−#
        # Solve the model
        tick()
        solveNLP = JuMP.optimize!       
        status = solveNLP(m)
        tock()

        # Get values for plotting
        xk = JuMP.value.(c[:,:,:]) # time series for plotting
        Residual=JuMP.value.(FO)
        teta_h=exp.(JuMP.value.(teta))
        

    
