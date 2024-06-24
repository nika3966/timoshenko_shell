include("iterations.jl")
include("app_sol_calc.jl")
include("p_err.jl")
include("pexact_sol.jl")
include("log_plot.jl")
#using CSV
#using DataFrames
using Base.Threads


@time begin
    # initial parameters
    N = 50 # number of summands in galerkin sum

    nu = 0.3
    k02 = 2/3
    E = 2.1e11
    h = 0.005  # 0.01, 0.02
    sg = 6*(1 - nu)*k02/(h^2)
    max_iter = 200
    reltol = 1e-8
    abstol = 1e-8
    al = 0.1

    n_layers = Integer(N/5)
    
    err_w_n = zeros(n_layers)
    err_u_n = zeros(n_layers)
    err_psi_n = zeros(n_layers)
    xx = collect(0:0.01:1)

    Threads.@threads for idx = 5:5:N
        
        idx_n = Integer(idx/5)
        # approximate solution coefficients and iteration count on layers
        (wn_m, iterations_count) = iterations(max_iter, k02, nu, idx, E, h, reltol, abstol,al)
        # error calculation
        (w_app, u_app, psi_app) = app_sol_calc(xx, wn_m, idx, reltol, abstol, sg)
	(pnw, pnu, pnpsi)  = pexact_sol(xx, idx, reltol, abstol, sg,al)
        (err_w, err_u, err_psi) = p_err(xx, pnw, pnu, pnpsi, w_app, u_app, psi_app)
        err_w_n[idx_n] = err_w
        err_u_n[idx_n] = err_u
        err_psi_n[idx_n] = err_psi
        
    end
    nn = 5:5:N
    fs = (6.5, 5)    
    y_min = 1e-5
    y_max = 1e-1 
    x_min = 0.4e+1
    x_max = 0.56e+2
    ylab = L"e_{n,m}^{(1)}"
    log_plot(fs, nn, err_u_n, h, x_min, x_max, y_min, y_max, ylab)
    #figure	
    y_min = 1e-6
    y_max = 1e-2 
    ylab = L"e_{n,m}^{(2)}"	   
    log_plot(fs, nn, err_w_n, h, x_min, x_max, y_min, y_max, ylab)    
    #figure
    y_min = 1e-5
    y_max = 1e-1 
    ylab = L"e_{n,m}^{(3)}"       
    log_plot(fs, nn, err_psi_n, h, x_min, x_max, y_min, y_max, ylab)
    # df = DataFrame(Galerkin_coeffs = nn, Algorithm_error = err_u_n)
    # CSV.write("res/err_u.csv", df)
	# df1 = DataFrame(Galerkin_coeffs = nn, Algorithm_error = err_w_n)
    # CSV.write("res/err_w.csv", df1)
	# df2 = DataFrame(Galerkin_coeffs = nn, Algorithm_error = err_psi_n)
    # CSV.write("res/err_psi.csv", df2)
end
