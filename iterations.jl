include("right_hand_fun.jl")
using LinearAlgebra

function iterations(max_iter, k02, nu, n, E, h, reltol, epsil,al)

  iterations_count = 0
  sg = 6.0*(1 - nu)*k02/(h^2)
  (fi, ki) = right_hand_fun(k02, E, h, sg, nu, n, reltol, epsil, al)
  nn = collect(1:1:n)
  wn_m = zeros(n)  
  bt0 = 1.0 ./(1.0/k02 .+ 6.0*(1.0 - nu)./((h*pi*nn).^2))
  s = (4*(1 - nu)/pi)*((ki./nn).*bt0) - 2*(1 - nu^2)/(E*h*pi)*(fi./nn)  
  dt = zeros(n)
  al = -6*ki./(pi*nn)
  for iter = 1:max_iter
    wn_m0 = copy(wn_m)	
    for id = 1:n
		  nn1 = deleteat!(copy(nn), id)
		  wn_m1 = deleteat!(copy(wn_m0), id)
		  ki1 = deleteat!(copy(ki), id)
		  dt[id] = pi^2*dot(nn1.^2,wn_m1.^2) - 4.0*dot(ki1,wn_m1)
	  end    
    r = dt + (2*(1 - nu))*bt0 - 4*(ki./(pi*nn)).^2
    aa = sqrt.((s/2).^2 + (r/3).^3)
    sigma1 = cbrt.(aa - s/2)
    sigma2 = cbrt.(aa + s/2)
    wn_m = (-al/3 + sigma1 - sigma2)./(pi*nn)       
    if maximum(abs.(wn_m - wn_m0)) < epsil
      iterations_count = iter
      break
    end
    iterations_count = iter
  end
  return wn_m, iterations_count

end