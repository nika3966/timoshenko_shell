# Copyright (C) 2024 nikoloz katchakhidze
# Created: 2024-03-21
using QuadGK
using LinearAlgebra

function pexact_sol(xx, n, reltol, epsil, sg, al)

  w(x) = al*exp(pi*x)*sin(pi*x)
  k(x) = exp(-pi*x)  
  lx = length(xx) 
  pnww = zeros(n)
  for id = 1:n
    pw1(x) = 2*w(x)*sin(id*pi*x)	  
    pnww[id], _ = quadgk(pw1, 0, 1, rtol=reltol, atol=epsil)         
  end
  nn = collect(1:1:n)
  pnw(x) = dot(pnww, sin.(pi*nn*x))
  a = pi*nn
  b = pnww.*a
  dpnw(x) = dot(b, cos.(pi*nn*x))
  pnu = zeros(lx)
  pnpsi = zeros(lx)
  sgr = sqrt(sg)
  fu(x) = -2*k(x)*pnw(x) + (dpnw(x))^2
  aa1 = -sgr/sinh(sgr)
  fpsi1(x) = aa1*cosh(sgr*(x - 1))
  fpsi2(x) = aa1*cosh(sgr*x)
  fpsi3(x) = cosh(sgr*x)*dpnw(x)
  fpsi4(x) = cosh(sgr*(x - 1))*dpnw(x)  
  for id = 1:lx
    temp1, _ = quadgk(fu, 0, xx[id], rtol=reltol, atol=epsil)
    temp2, _ = quadgk(fu, xx[id], 1, rtol=reltol, atol=epsil)      
    pnu[id] = 0.5*(xx[id] - 1)*temp1 + 0.5*xx[id]*temp2
    temp3, _ = quadgk(fpsi3, 0, xx[id], rtol=reltol, atol=epsil)
    temp4, _ = quadgk(fpsi4, xx[id], 1, rtol=reltol, atol=epsil)      
    pnpsi[id] = fpsi1(xx[id])*temp3 + fpsi2(xx[id])*temp4	
  end  
  pnww = pnw.(xx)  
  return pnww, pnu, pnpsi

end
