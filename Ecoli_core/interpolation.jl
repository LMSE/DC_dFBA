radau  = [0.15505 0.64495 1.00000]
xnc = [0.1 0.15505 0.2 0.3 0.4 0.5 0.6 0.64495 0.7 0.8 0.9 1.0]
xk_i = Matrix{Float64}(undef,nc,nfe*size(xnc,2)+1)
tsn_i = Vector{Float64}(undef,nfe*size(xnc,2)+1)
theta = Matrix{Float64}(undef,ncp,size(xnc,2))

tsn_i[1] = 0.0
xk_i[:,1]=c0.*cs

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

tsn_i[k] = ts[i] + xnc[j]*solh[i]

global k=k+1
end

end
