using QuadGK
using LinearAlgebra

function app_sol_calc(xx, wn_m, n, reltol, epsil, sg)

  sgr = sqrt(sg)
  k(x) = exp(-pi*x)
  nn = collect(1:1:n)
  w_app(x) = dot(wn_m, sin.(pi*nn*x))
  dw_app(x) = dot(pi*(wn_m.*nn), cos.(pi*nn*x))  
  lx = length(xx)  
  u_app = zeros(lx)
  uu(x) = -2*k(x)*w_app(x) + (dw_app(x))^2  
  psi_app = zeros(lx)
  aa1 = -sgr/sinh(sgr)
  psi_app1(x) = aa1*cosh(sgr*(x - 1))
  psi_app2(x) = aa1*cosh(sgr*x)
  psi_app3(x) = cosh(sgr*x)*dw_app(x)
  psi_app4(x) = cosh(sgr*(x - 1))*dw_app(x)
  for id = 1:lx
    temp1, _ = quadgk(uu, 0, xx[id], rtol=reltol, atol=epsil)
    temp2, _ = quadgk(uu, xx[id], 1, rtol=reltol, atol=epsil)      
    u_app[id] = 0.5*(xx[id] - 1)*temp1 + 0.5*xx[id]*temp2
	  temp3, _ = quadgk(psi_app3, 0, xx[id], rtol=reltol, atol=epsil)
    temp4, _ = quadgk(psi_app4, xx[id], 1, rtol=reltol, atol=epsil)      
    psi_app[id] = psi_app1(xx[id])*temp3 + psi_app2(xx[id])*temp4
  end
  w_app1 = w_app.(xx)
  return w_app1, u_app, psi_app

end
