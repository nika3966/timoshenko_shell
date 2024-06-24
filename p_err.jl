using LinearAlgebra

function p_err(xx, pnw, pnu, pnpsi, w_app, u_app, psi_app)
  
  deltaw2 = (pnw - w_app).^2
  deltau2 = (pnu - u_app).^2
  deltapsi2 = (pnpsi - psi_app).^2  
  lx = length(xx)
  hh = 1/lx
  err_w = sqrt(hh*sum(deltaw2))
  err_u = sqrt(hh*sum(deltau2))
  err_psi = sqrt(hh*sum(deltapsi2))  
  return err_w, err_u, err_psi

end
