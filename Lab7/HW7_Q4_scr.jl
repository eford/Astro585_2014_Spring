Pkg.add("HDF5")
using HDF5, JLD

@everywhere include("HW6_Q2_planet_populations_once.jl")  
@everywhere include("HW6_Q2_planet_populations.jl")
min_eta = 0.0
max_eta = 1.0
min_shape = 0.0001
max_shape = 1.0
min_log_scale = log10(0.3)
max_log_scale = log10(3.0)  
num_eta = 4
num_shape = 4
num_scale = 4
num_evals = 1
etas = linspace(min_eta,max_eta,num_eta)
scales = 10.0.^linspace(min_log_scale,max_log_scale,num_scale)
shapes = linspace(min_shape,max_shape,num_shape)
num_stars = 16000
eta = 0.2
shape = 0.1
scale = 1.0;

srand(42)
@time result = eval_model_on_grid_loops(etas,shapes,scales,num_stars; num_evals = num_evals, true_eta = eta, true_shape = shape, true_scale = scale);

@save "result.jld" result, min_eta, max_eta, min_shape, max_shape, min_log_scale, max_log_scale, num_eta, num_shape, num_scale, num_evals, etas, scales, shapes, num_stars, eta, shape, scale

