# Copyright (C) 2024 nikoloz katchakhidze
# Created: 2024-03-21
include("pexact_sol.jl")
include("app_sol_calc.jl")
include("p_err.jl")
include("plot_res.jl")
using PyPlot

function output_results_py(xx, wn_m, n, reltol, epsil, sg, fid, it_count, al)
  
  (w_exx, u_exx, psi_exx) = pexact_sol(xx, n, reltol, epsil, sg, al)
  (w_app, u_app, psi_app) = app_sol_calc(xx, wn_m, n, reltol, epsil, sg)
  (err_w, err_u, err_psi) = p_err(xx, w_exx, u_exx, psi_exx, w_app, u_app, psi_app)     
  d = [1, 26, 51, 76, 101]   
  nodes = getindex(xx, d)  
  u_ex = getindex(u_exx, d)
  w_ex = getindex(w_exx, d)
  psi_ex = getindex(psi_exx, d)
  u_out = getindex(u_app, d)
  psi_out = getindex(psi_app, d)
  w_out = getindex(w_app, d)
  F = open(fid,"a")
  write(F, "Number of Galerkin coefficients: n = " * string(n) * "\n")
  write(F, "iteration_count is " * string(it_count) * "\n")
  write(F, "Algorithm error for w(x) is " * string(err_w) * "\n")
  write(F, "Algorithm error for u(x) is " * string(err_u) * "\n")
  write(F, "Algorithm error for psi(x) is " * string(err_psi) * "\n")
  write(F, "\n")
  write(F, "The nodes")
  write(F, "\n");
  write(F, "  " * string(nodes))
  write(F, "\n\n")
  write(F, "The values of exact solution w_ex")
  write(F, "\n")
  write(F, "  " * string(w_ex))
  write(F, "\n")
  write(F, "The values of approximate solution w_app")
  write(F, "\n")
  write(F, "  " * string(w_out))
  write(F, "\n\n")
  write(F, "The values of exact solution u_ex")
  write(F, "\n")
  write(F, "  " * string(u_ex))
  write(F, '\n')
  write(F, "The values of approximate solution u_app")
  write(F, "\n")
  write(F, "  " * string(u_out))
  write(F, "\n\n")
  write(F, "The values of exact solution psi_ex")
  write(F, "\n");
  write(F, "  " * string(psi_ex))
  write(F, "\n")
  write(F, "The values of approximate solution psi_app")
  write(F, "\n")
  write(F, "  " * string(psi_out))
  write(F, "\n")
  close(F)
  
  fs = (6,5)
  ylab = L"\psi"
  uap = L"u_{10,138}(x)"
  wap = L"w_{10,138}(x)"
  psi_ap = L"\psi_{10,138}(x)"
  psi_ex = L"\psi(x)"
  plot_res(fs, xx, w_exx, w_app, "w(x)", wap, "w")
  plot_res(fs, xx, u_exx, u_app, "u(x)", uap, "u")
  plot_res(fs, xx, psi_exx, psi_app, psi_ex, psi_ap, ylab)

  
  
end



