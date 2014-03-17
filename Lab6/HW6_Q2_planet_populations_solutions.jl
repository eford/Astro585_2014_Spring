
function evaluate_model_parallel(stats_obs, eta::Real, shape::Real, scale::Real, num_stars::Integer; minP=0.0, maxP = Inf, num_evals::Integer = 1 )
  @assert(0.0<=eta<=1.0)
  @assert(num_stars>0)
  @assert(minP>=0.0)
  @assert(maxP>minP)
  @assert(num_evals>=1)  
  sum_distance = @parallel (+)  for i in 1:num_evals
      data_sim = generate_transiting_planet_sample(eta, shape, scale, num_stars; minP=minP, maxP=maxP)
      stats_sim = compute_stats(data_sim)
      evaluate_model(stats_obs,stats_sim)
  end
  sum_distance / num_evals      
end


function eval_model_on_grid_parallel(etas::Array, shapes::Array, scales::Array, num_stars = 1600; num_evals = 1, 
                            true_eta = 0.2, true_shape = 0.1, true_scale = 1.0)
const solar_radius_in_AU = 0.00464913034
minP = solar_radius_in_AU^1.5
maxP = 4*days_in_year/3
data_obs = generate_transiting_planet_sample(true_eta, true_shape, true_scale, num_stars; minP=minP, maxP=4*days_in_year/3)
stats_obs = compute_stats(data_obs)

dist = Array(Float64,(length(etas), length(shapes), length(scales) ) )
for k in 1:length(scales)
    for j in 1:length(shapes)
        for i in 1:length(etas)
            dist[i,j,k] = evaluate_model_parallel(stats_obs, etas[i], shapes[j], scales[k], num_stars; minP=minP, maxP=maxP, num_evals=num_evals )
        end # for i
    end # for j
end # for k
dist
end # function

function eval_model_on_grid_map(etas::Array, shapes::Array, scales::Array, num_stars = 1600; num_evals = 1, 
                            true_eta = 0.2, true_shape = 0.1, true_scale = 1.0)
const solar_radius_in_AU = 0.00464913034
minP = (2.0*solar_radius_in_AU)^1.5
maxP = 4*days_in_year/3
data_obs = generate_transiting_planet_sample(true_eta, true_shape, true_scale, num_stars; minP=minP, maxP=4*days_in_year/3)
stats_obs = compute_stats(data_obs)

idx = [ (i,j,k) for i in 1:length(etas), j in 1:length(shapes), k in 1:length(scales) ]
map(tuple -> evaluate_model(stats_obs, etas[tuple[1]], shapes[tuple[2]], scales[tuple[3]], num_stars; minP=minP, maxP=maxP, num_evals=num_evals ), idx)

end # function


# returns distributed array
function eval_model_on_grid_map_distributed(etas::Array, shapes::Array, scales::Array, num_stars = 1600; num_evals = 1,  true_eta = 0.2, true_shape = 0.1, true_scale = 1.0)
const solar_radius_in_AU = 0.00464913034
minP = solar_radius_in_AU^1.5
maxP = 4*days_in_year/3
data_obs = generate_transiting_planet_sample(true_eta, true_shape, true_scale, num_stars; minP=minP, maxP=4*days_in_year/3)
stats_obs = compute_stats(data_obs)

idx = distribute([ (i,j,k) for i in 1:length(etas), j in 1:length(shapes), k in 1:length(scales) ])
map(tuple -> evaluate_model(stats_obs, etas[tuple[1]], shapes[tuple[2]], scales[tuple[3]], num_stars; minP=minP, maxP=maxP, num_evals=num_evals ), idx)

end # function

