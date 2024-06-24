include("iterations.jl")
#include("iterations1.jl")
include("output_results_py.jl")

# initial parameters
n = 20 # number of summands in galerkin sum

nu = 0.3
k02 = 2/3
E = 2.1e11
h = 0.005  # 0.01, 0.02
sg = 6*(1 - nu)*k02/(h^2)
max_iter = 500
reltol = 1e-8
epsil = 1e-8
xx = collect(0:0.01:1)
al = 0.1
# approximate solution
(wn_m, iterations_count) = iterations(max_iter, k02, nu, n, E, h, reltol, epsil,al)
# output_results
fid = "results" * string(n) * ".txt"
F = open(fid,"w");
write(F,"Initial parameters:")
write(F,"\n\n")
write(F,"k_0^2 = " * string(k02) * "\n")
write(F,"nu = " * string(nu) * "\n")
write(F,"E = " * string(E) * "\n")
write(F,"\n")
write(F,"The solution:")
write(F,"\n\n")
close(F)

output_results_py(xx, wn_m, n, reltol, epsil, sg, fid, iterations_count,al)