function pFBA_KKT_flux(c0,v0)

    #−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−#
    # Collocation parameters and radau time series
    #N = M^-1 =
    colmat = [0.19681547722366  -0.06553542585020 0.02377097434822;
              0.39442431473909  0.29207341166523 -0.04154875212600;
              0.37640306270047  0.51248582618842 0.11111111111111]
    radau  = [0.15505 0.64495 1.00000]

    #−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−#
    # JuMP model
    m = Model(with_optimizer(Ipopt.Optimizer, warm_start_init_point = "yes", print_level = 5, linear_solver = "ma27",tol=1e-4,acceptable_iter=5,acceptable_tol=1e-2) )

    #m = Model(with_optimizer(Ipopt.Optimizer, warm_start_init_point = "yes", print_level = 5, linear_solver = "ma27") )
    
  

    # Set up variables for JuMP
    @variables(m, begin
        c[1:nc, 1:nfe, 1:ncp]     #States: 
        cdot[1:nc, 1:nfe, 1:ncp]  #Differential DifferentialEquations
        v[1:nv, 1:nfe] #fluxes                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
        lambda[1:nm, 1:nfe]  #multipliers
        alpha_U[1:nv, 1:nfe]
        alpha_L[1:nv, 1:nfe]
        alpha_upt[1:nfe]

        FO_U[1:nv, 1:nfe]
        FO_L[1:nv, 1:nfe]
        FO_upt[1:nfe]

        hv[1:nfe]
        reg[1]
   
    end)

    # Set up initial guesses for solver
    for i in 1:nfe
        for j in 1:ncp
            
            for k in 1:nc
            #set_start_value(c[k,i,j], c0[k])      #C_A
            end
             
            for k in 1:nv
            #set_start_value(v[k,i,j], v0[k])      #v
            end
            
        end

    end


    for i in 1:nfe
      set_start_value(hv[i], hm[i])

    end

    # Set up objective function

    @NLexpressions(m, begin
    
    #vg[i=1:nfe], vg_max*(c[1,i,3]/(Kg+(c[1,i,3])) )

    vg[i=1:nfe], qmax*(c[1,i,3]/(Kg+(c[1,i,3])) )
    
    #vg[i=1:nfe], qmin + (qmax - qmin)*(-1 + 2/(1 + exp(-Kup*v[obj,i])) )

    end)
  
  #@NLobjective(m, Min,  sum( sum(-phi1*FO_L[mc,i] + -phi3*FO_U[mc,i] for mc in 1:nv) + phi2*FO_upt[i] for i in 1:nfe))

  @NLobjective(m, Min, omega*-c[4,end,end] + sum( sum(-phi1*FO_L[mc,i] + -phi3*FO_U[mc,i] for mc in 1:nv) + phi2*FO_upt[i] for i in 1:nfe))

  #@NLobjective(m, Min, sum(sum(omega*-c[4,i,j] for i in 1:nfe) for j in 1:ncp) + sum( sum(-phi1*FO_L[mc,i] + -phi3*FO_U[mc,i] for mc in 1:nv) + phi2*FO_upt[i] for i in 1:nfe))

  #@NLobjective(m, Max, c[4,end,end] )
  
  #@objective(m, Min, 1)
  
  
    #Set up the constraints
    @constraints(m, begin
        # set up differential equations
        m1[i=1:nfe, j=1:ncp], cdot[1,i,j] == v[glu,i]*c[3,i,j]
        m2[i=1:nfe, j=1:ncp], cdot[2,i,j] == v[ac,i]*c[3,i,j]
        m3[i=1:nfe, j=1:ncp], cdot[3,i,j] == v[obj,i]*c[3,i,j]
        m4[i=1:nfe, j=1:ncp], cdot[4,i,j] == v[lac,i]*c[3,i,j]

        # set up collocation equations - 2nd-to-nth point
        coll_c_n[l=1:nc, i=2:nfe, j=1:ncp], c[l,i,j] == c[l,i-1,ncp]+hv[i]*sum(colmat[j,k]*cdot[l,i,k] for k in 1:ncp)
        
        # set up collocation equations - 1st point
        coll_c_0[l=1:nc, i=1, j=1:ncp], c[l,i,j] == c0[l] + hv[i]*sum(colmat[j,k]*cdot[l,i,k] for k in 1:ncp)


        #------------- system constraints --------------#
        Sc[mc=1:nm,i=1:nfe],  sum(S[mc,k]*v[k,i] for k in 1:nv) == 0
        v_UB[mc=1:nv,i=1:nfe], v[mc,i] - vub[mc] <= 0
        v_LB[mc=1:nv,i=1:nfe], -v[mc,i] + vlb[mc] <= 0

        c_LB[mc=1:nc,i=1:nfe, j=1:ncp], -c[mc,i,j]  <= 0

        MFE1, sum(hv[i] for i in 1:nfe) == th
        MFE3[i=1:nfe], hv[i]  >= 0.0
        MFE4[i=1:nfe], hv[i]  >= (1-var_h)*hm[1]
        MFE5[i=1:nfe], hv[i]  <= (1+var_h)*hm[1]

        Lagr[mc=1:nv,i=1:nfe], +d[mc] + w*v[mc,i] + alpha_L[mc,i] + alpha_U[mc,i] + up[mc]*alpha_upt[i] + sum(S[k,mc]*lambda[k,i] for k in 1:nm) == 0

        alpha1_LB[mc=1:nv,i=1:nfe], alpha_L[mc,i] <= 0
        alpha4_LB[i=1:nfe], alpha_upt[i] <= 0
        alpha1_UB[mc=1:nv,i=1:nfe], alpha_U[mc,i] >= 0    

        #regu[i=treg:nfe], v[obj,i]  <= reg[1]*0.229

        #regu[i=1:nfe], v[pta,i] - reg[i]*vub[pta] <= 0
        #regu2[1],   reg[1] <= 1
        #regu3[1],   reg[1] >= 0
        #regu4[i=2:nfe],   reg[i] - reg[i-1] <= 0

        regu[i=2:Int(nfe/2)],   v[obj,i] - v[obj,i-1] == 0
        regu2[i=Int((nfe/2)+2):nfe],   v[obj,i] - v[obj,i-1] == 0
       # v_reg[mc=obj,i=1:nfe], -v[mc,i] + 0.2 <= 0
    
        
    end)

    @NLconstraints(m, begin
    
      v_LB_g[i=1:nfe], -v[glu,i] - vg[i] <= 0  
      
      FO1[mc=1:nv,i=1:nfe], FO_L[mc,i] == (v[mc,i]-vlb[mc])*alpha_L[mc,i]
      FO2[mc=1:nv,i=1:nfe], FO_U[mc,i] == (v[mc,i]-vub[mc])*alpha_U[mc,i]
      FO3_upt[i=1:nfe], FO_upt[i] == (-v[glu,i]-vg[i])*alpha_upt[i]

    end)

    #−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−#
    # Solve the model
    #tick()
    solveNLP = JuMP.optimize!
    status = solveNLP(m)
    #tock()

    # Get values for plotting
    cStar = JuMP.value.(c[:,:,:]) # time series for plotting
    vStar = JuMP.value.(v[:,:])
    cdStar = JuMP.value.(cdot[:,:,:])
    lStar = JuMP.value.(lambda[:,:])
    alStar = JuMP.value.(alpha_L[:,:])
    agStar = JuMP.value.(alpha_upt[:])
    hStar = JuMP.value.(hv[:])
    rStar = JuMP.value.(reg[:])

return  cStar, vStar, cdStar, lStar, alStar, agStar, hStar, rStar

end