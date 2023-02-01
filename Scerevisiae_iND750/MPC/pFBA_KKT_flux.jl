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

  # Set up variables for JuMP
  @variables(m, begin
      c[1:nc, 1:nfe, 1:ncp]     #States: 
      cdot[1:nc, 1:nfe, 1:ncp]  #Differential DifferentialEquations
      v[1:nv, 1:nfe] #fluxes                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
      lambda[1:nm, 1:nfe]  #multipliers
      alpha_U[1:nv, 1:nfe]
      alpha_L[1:nv, 1:nfe]
      alpha_upt[1:nfe,1:nupt]

      FO_U[1:nv, 1:nfe]
      FO_L[1:nv, 1:nfe]
      FO_upt[1:nfe,1:nupt]

      hv[1:nfe]
      uk[1:ncv,1:nfe]
     # du[1:ncv,1:nfe]
 
  end)

  # Set up initial guesses for solver
  for i in 1:nfe
      for j in 1:ncp
          
          for k in 1:nc
          set_start_value(c[k,i,j], c0[k])      #C_A
          end
           
          for k in 1:nv
          #set_start_value(v[k,i,j], v0[k])      #v
          end
          
      end

  end


  for i in 1:nfe
    set_start_value(hv[i], hm[i])

  end

   #scaling
 for i in 1:nc
  c0[i]=c0[i]/cs[i]

end

  # Set up objective function

  @NLexpressions(m, begin


  vg[i=1:nfe], vg_max*(  ((c[2,i,3]*cs[2])/(c[4,i,3]*cs[4]))/(Kg+((c[2,i,3]*cs[2])/(c[4,i,3]*cs[4]))+
  ((((c[2,i,3]*cs[2])/(c[4,i,3]*cs[4]))^2)/Kig)) )*(1/(1+(((c[3,i,3]*cs[3])/(c[4,i,3]*cs[4]))/Kie)))


  vo[i=1:nfe], vo_max*((uk[1,i]*us[1])/(Ko+(uk[1,i]*us[1])) )

  
  end)


@NLobjective(m, Max, c[3,end,end] - sum( sum(-phi1*FO_L[mc,i] + -phi3*FO_U[mc,i] for mc in 1:nv) + sum(phi2*FO_upt[i,upt] for upt in 1:nupt) for i in 1:nfe))



  #Set up the constraints
  @constraints(m, begin
      # set up differential equations
      m1[i=1:nfe, j=1:ncp], cdot[1,i,j] == vs[obj]*v[obj,i]*c[1,i,j]
      m2[i=1:nfe, j=1:ncp], cdot[2,i,j] == (us[2]*uk[2,i]*Gf + vs[glu]*v[glu,i]*c[1,i,j]*cs[1])/cs[2]
      m3[i=1:nfe, j=1:ncp], cdot[3,i,j] == (vs[eth]*v[eth,i]*c[1,i,j]*cs[1])/cs[3]
      m4[i=1:nfe, j=1:ncp], cdot[4,i,j] == us[2]*uk[2,i]/cs[4]

      # set up collocation equations - 2nd-to-nth point
      coll_c_n[l=1:nc, i=2:nfe, j=1:ncp], c[l,i,j] == c[l,i-1,ncp]+hv[i]*sum(colmat[j,k]*cdot[l,i,k] for k in 1:ncp)
      
      # set up collocation equations - 1st point
      coll_c_0[l=1:nc, i=1, j=1:ncp], c[l,i,j] == c0[l] + hv[i]*sum(colmat[j,k]*cdot[l,i,k] for k in 1:ncp)


      #------------- system constraints --------------#
      Sc[mc=1:nm,i=1:nfe],  sum(S[mc,k]*v[k,i]*vs[k] for k in 1:nv) == 0
      v_UB[mc=1:nv,i=1:nfe], v[mc,i]*vs[mc] - vub[mc]  <= 0
      v_LB[mc=1:nv,i=1:nfe], -v[mc,i]*vs[mc]  + vlb[mc] <= 0

      c_LB[mc=1:nc,i=1:nfe, j=1:ncp], -c[mc,i,j]  <= 0.0
      c_UB_V[i=1:nfe, j=1:ncp], c[4,i,j]  <= 1.2

      MFE1, sum(hv[i] for i in 1:nfe) == th
      MFE3[i=1:nfe], hv[i]  >= 0.0
      MFE4[i=1:nfe], hv[i]  >= (1-var_h)*hm[1]
      MFE5[i=1:nfe], hv[i]  <= (1+var_h)*hm[1]

      #pFBA
      Lagr[mc=1:nv,i=1:nfe], +d[mc] + w*v[mc,i]*vs[mc]  + alpha_L[mc,i] + alpha_U[mc,i] + up[mc]*alpha_upt[i,1] + up2[mc]*alpha_upt[i,2] + sum(S[k,mc]*lambda[k,i] for k in 1:nm) == 0
      #FBA
     # Lagr[mc=1:nv,i=1:nfe], +d[mc] + alpha_L[mc,i] + alpha_U[mc,i] + up[mc]*alpha_upt[i] + sum(S[k,mc]*lambda[k,i] for k in 1:nm) == 0

      alpha1_LB[mc=1:nv,i=1:nfe], alpha_L[mc,i] <= 0
      alpha4_LB[i=1:nfe,upt=1:nupt], alpha_upt[i,upt] <= 0
      alpha1_UB[mc=1:nv,i=1:nfe], alpha_U[mc,i] >= 0    

      cons_uk[mc=1:ncv,i=1:nfe], uk[mc,i] <= uk_max[mc]
      cons_uk2[mc=1:ncv,i=1:nfe], uk[mc,i] >= uk_min[mc]
      #cons_uk3[mc=1:ncv,i=1], uk[mc,i] == u0[mc]

       # set up du - manipulated variable change
      # deltaU_n[mc=1:ncv, i=2:nfe], du[mc,i] == uk[mc,i] - uk[mc,i - 1]
       #deltaU_0[mc=1:ncv, i=1], du[mc,i] == uk[mc,i] - u0[mc]
       # set up du bounds
       #du_UB[mc=1:ncv, i=1:nfe], du[mc,i]/DuMax[mc] - 1 <= 0
       #du_LB[mc=1:ncv, i=1:nfe], -1 - du[mc,i]/DuMax[mc] <= 0
      
  end)

  @NLconstraints(m, begin
  
    v_LB_g[i=1:nfe], -v[glu,i]*vs[glu]  - vg[i] <= 0  
    v_LB_o[i=1:nfe], -v[o2,i]*vs[o2]  - vo[i] <= 0  
    
    FO1[mc=1:nv,i=1:nfe], FO_L[mc,i] == (v[mc,i]*vs[mc] -vlb[mc])*alpha_L[mc,i]
    FO2[mc=1:nv,i=1:nfe], FO_U[mc,i] == (v[mc,i]*vs[mc] -vub[mc])*alpha_U[mc,i]
    FO3_upt[i=1:nfe], FO_upt[i,1] == (-v[glu,i]*vs[glu] -vg[i])*alpha_upt[i,1]
    FO4_upt[i=1:nfe], FO_upt[i,2] == (-v[o2,i]*vs[o2] -vo[i])*alpha_upt[i,2]

  end)

  #−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−#
  # Solve the model
  #tick()
  solveNLP = JuMP.optimize!
  status = solveNLP(m)
  #tock()

  # Get values for plotting
  cStar = JuMP.value.(c[:,:,:]).*cs # time series for plotting
  vStar = JuMP.value.(v[:,:]).*vs
  cdStar = JuMP.value.(cdot[:,:,:])
  lStar = JuMP.value.(lambda[:,:])
  alStar = JuMP.value.(alpha_L[:,:])
  agStar = JuMP.value.(alpha_upt[:,:])
  hStar = JuMP.value.(hv[:])
  uStar = JuMP.value.(uk[:,:]).*us[:]

return  cStar, vStar, cdStar, lStar, alStar, agStar, hStar, uStar

end