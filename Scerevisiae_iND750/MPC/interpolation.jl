radau  = [0.15505 0.64495 1.00000]
#xnc = [0.1 0.15505 0.2 0.3 0.4 0.5 0.6 0.64495 0.7 0.8 0.9 1.0]
xnc = [0.1 0.15 0.15505 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.64495 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0]
xk_i = Matrix{Float64}(undef,nc,nfe*size(xnc,2)+1)
uk_i = Matrix{Float64}(undef,ncv,nfe*size(xnc,2)+1)
tsn_i = Vector{Float64}(undef,nfe*size(xnc,2)+1)
theta = Matrix{Float64}(undef,ncp,size(xnc,2))


tsn_i[1] = 0.0
xk_i[:,1]=c0.*cs
uk_i[:,1]=uk[:,1]

global k=2

for j in 1:size(xnc,2)
theta[1,j] = ((xnc[j]-radau[2])*(xnc[j]-radau[3])/((radau[1]-radau[2])*(radau[1]-radau[3])))
theta[2,j] = ((xnc[j]-radau[1])*(xnc[j]-radau[3])/((radau[2]-radau[1])*(radau[2]-radau[3])))
theta[3,j] = ((xnc[j]-radau[1])*(xnc[j]-radau[2])/((radau[3]-radau[1])*(radau[3]-radau[2])))
end

for i in 1:nfe
for j in 1:size(xnc,2)

for l in 1:nc
xk_i[l,k] = theta[1,j]*sol[l,i,1] + 
theta[2,j]*sol[l,i,2] +
theta[3,j]*sol[l,i,3] 

end

for ki in 1:ncv
  uk_i[ki,k] = uk[ki,i]
end

tsn_i[k] = ts[i] + xnc[j]*solh[i]

global k=k+1
end

end

xk_i[2,:] = xk_i[2,:].*conv_glu
xk_i[3,:] = xk_i[3,:].*conv_eth

uk_lab_i = Matrix{Float64}(undef,ncv,10*nfe)
ts_lab_i = Vector{Float64}(undef,nfe)
ts_lab_ii = Vector{Float64}(undef,10*nfe)

global kkl=1
for i in 1:nfe
  ts_lab_i[i] = i*th/nfe
  for j in 1:10

       uk_lab_i[:,kkl] = uk_lab[:,i]
       ts_lab_ii[kkl] = ts_lab_i[i] -h +  j*h/10
       global kkl=kkl+1
  end
 
  
end