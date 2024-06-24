# Copyright (C) 2024 nikoloz katchakhidze
# Created: 2024-03-21
using QuadGK

function right_hand_fun(k02, E, h, sg, nu, n, reltol, abstol, al)
  
  k(x) = exp(-pi*x)
  sgr = sqrt(sg)
  a1 = E*h*al*(2/pi + pi*al*(1 - exp(2*pi))/8)/(1 - nu^2)
  f0(x) = 2*pi^2*al*exp(pi*x)*cos(pi*x)
  f1(x) = a1*(k(x) + f0(x))
  a2 = k02*E*h/(2*(1 + nu))
  a3 = 2*pi^2*sg*al/(((pi + sgr)^2 + pi^2)*((pi - sgr)^2 + pi^2))  
  f2(x) = exp(pi*x)*(2*pi^2*sin(pi*x) - sg*cos(pi*x))
  a4 = sg/sinh(sgr)
  f3(x) = a4*(exp(pi)*sinh(sgr*x) + sinh(sgr*(x - 1)))  
  f(x) = f1(x) - a2*(f0(x) + a3*(f2(x) - f3(x)))  
  fi = zeros(n)
  ki = zeros(n)
  for id = 1:n
      ff(x) = 4*f(x)*sin(id*pi*x)
      fk(x) = k(x)*sin(id*pi*x)
      fi[id], _ = quadgk(ff, 0, 1, rtol=reltol, atol=abstol)
      ki[id], _ = quadgk(fk, 0, 1, rtol=reltol, atol=abstol)      
  end
  return fi, ki

end
