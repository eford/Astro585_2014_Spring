#Kepler's Law
semimajor_axis(P::Real, M::Real, G::Real=4pi*pi) = (G*M/(4pi*pi)*(P/days_in_year)^2)^(1.0/3.0)
@vectorize_2arg Real semimajor_axis

# Simplistic functions to generate star properties
function generate_star_mass_sample(N::Integer)
  @assert(N>0)
  fill(1.0,N)
end

function generate_star_radius_sample(N::Integer)
  @assert(N>0)
  solar_radius_in_AU = 0.00464913034
  fill(1.0*solar_radius_in_AU,N)
end

# Simplistic function pick a fraction of targets to have a planet (can also be used with an eta array)
function generate_number_of_planets_sample(eta, N::Integer)
  @assert(N>=0)
  num_planets = zeros(Int64,N)                           # Set all to have no planets
  idx_w_planets = rand(N) .<= eta                        # single probability of having exactly one planet
  num_planets[idx_w_planets] = 1                         # Indicate which have a planet
  num_planets  
end

# Draw period from truncated inverse gamma function
function generate_planet_period(shape::Real, scale::Real; minP = 0.0, maxP = Inf )
  @assert(shape>0.0)
  @assert(scale>0.0)
  @assert(minP>=0.0)
  @assert(maxP>minP)
  period_distribution = InverseGamma(shape,scale)  
  period = -1.0
  while !(minP <= period <= maxP)
    period = rand(period_distribution)
  end
  period
end

function generate_planet_period_sample(shape::Real, scale::Real, N::Integer; minP = 0.0, maxP = Inf )
  periods = [ generate_planet_period(shape,scale;minP=minP,maxP=maxP) for i in 1:N ]
end

# Returns a tuple of lists of orbital periods in days for planets, star masses and star radii for planet host stars
# eta = Fraction of stars with a planet
# Period distribution ~ Inverse Gamma distribution
#   pdf(P) = scale^shape /Gamma(shape)  P^(-shape-1) exp(-scale/P) 
# num_stars = Number of target stars
function generate_planet_sample(eta::Real, shape::Real, scale::Real, num_stars::Integer; minP = 0.0, maxP = Inf )
  @assert(0.0<=eta<=1.0)
  @assert(minP>=0.0)
  @assert(maxP>minP)
  @assert(num_stars>0)
  # Generate star properties
  mstar_list_all = generate_star_mass_sample(num_stars)
  rstar_list_all = generate_star_radius_sample(num_stars)
  num_planets_all = generate_number_of_planets_sample(eta,num_stars)
  total_planets = sum(num_planets_all)                         # Total number of planets around all targets
    
  # Select stars w/ planets
  idx_w_planets = find(num_planets_all)
  mstar_list = mstar_list_all[idx_w_planets]
  rstar_list = rstar_list_all[idx_w_planets]
  
  # Generate planet properties
  period_list = generate_planet_period_sample(shape,scale,total_planets;minP=minP,maxP=maxP)
  (period_list, mstar_list, rstar_list)
end

# Returns list of orbital periods in days for transiting planets
# eta = Fraction of stars with a planet
# Period distribution ~ Inverse Gamma distribution
#   pdf(P) = scale^shape /Gamma(shape)  P^(-shape-1) exp(-scale/P) 
# num_stars = Number of target stars
function generate_transiting_planet_sample(eta::Real, shape::Real, scale::Real, num_stars::Integer; minP = 0.0, maxP = Inf )
  @assert(0.0<=eta<=1.0)
  @assert(num_stars>0)
  (period_list, mstar_list, rstar_list) = generate_planet_sample(eta,shape,scale,num_stars; minP=minP, maxP=maxP)
  @assert(length(period_list) == length(mstar_list) == length(rstar_list) )
  total_planets = length(period_list)
  
  transit_probs = min(1.0,(rstar_list./semimajor_axis(period_list,mstar_list)))  # Assumes circular orbits, R/a<<1, independent orientations
  num_transiting_planets = generate_number_of_planets_sample(transit_probs, total_planets)   # Number of transiting planets for each target w/ a transiting planet
  idx_w_transiting_planets = find(num_transiting_planets)
  period_list[idx_w_transiting_planets]   
end




# returns a list of stats about the samples from data_obs
function compute_stats(data_obs::Array)
  logP = data_obs   # for now, just one array 
  N = length(logP)
  logP = log10(data_obs)
  mean_logP = (N>=1) ? mean(logP) : 0.0
  stddev_logP = (N>=2) ? std(logP) : 0.0
  skewness_logP = (N>=3) ? skewness(logP) : 0.0
  return [N, mean_logP, stddev_logP, skewness_logP] 
end

# Calculates a "distance" between observations and model using stats_obs (precomputed from observed sample) and stats calculated from sample in stats_sim
function evaluate_model(stats_obs::Array, stats_sim::Array )
  @assert(typeof(stats_sim) == typeof(stats_obs))
  @assert(size(stats_sim) == size(stats_obs))
  idx_rel = [ 1, 2, 3 ]  # Use relative difference for N, mean, and stddev
  idx_abs = [ 4 ]  # Don't use relative difference for skewness, since near zero
  max_d_rel = maximum(abs((stats_sim[idx_rel].-stats_obs[idx_rel])./stats_obs[idx_rel]))
  max_d_abs = maximum(abs((stats_sim[idx_abs].-stats_obs[idx_abs])))
  max_d = max(max_d_rel,max_d_abs)
end

# Calculates a "distance" between observations and draw(s) from a model 
# Optionally calculates an average distance using num_evals different realizations from same model
function evaluate_model(stats_obs, eta::Real, shape::Real, scale::Real, num_stars::Integer; minP=0.0, maxP = Inf, num_evals::Integer = 1 )
  @assert(0.0<=eta<=1.0)
  @assert(num_stars>0)
  @assert(minP>=0.0)
  @assert(maxP>minP)
  @assert(num_evals>=1)  
  sum_distance = 0.0
  for i in 1:num_evals
      data_sim = generate_transiting_planet_sample(eta, shape, scale, num_stars; minP=minP, maxP=maxP)
      stats_sim = compute_stats(data_sim)
      sum_distance += evaluate_model(stats_obs,stats_sim)
  end
  sum_distance / num_evals      
end


function eval_model_on_grid_loops(etas::Array, shapes::Array, scales::Array, num_stars = 1600; num_evals = 1, 
                            true_eta = 0.2, true_shape = 0.1, true_scale = 1.0)
  const solar_radius_in_AU = 0.00464913034
  minP = (2.0*solar_radius_in_AU)^1.5
  maxP = 4*days_in_year/3
  data_obs = generate_transiting_planet_sample(true_eta, true_shape, true_scale, num_stars; minP=minP, maxP=4*days_in_year/3)
  stats_obs = compute_stats(data_obs)

  dist = Array(Float64,(length(etas), length(shapes), length(scales) ) )
  for k in 1:length(scales)
    for j in 1:length(shapes)
      for i in 1:length(etas)
        dist[i,j,k] = evaluate_model(stats_obs, etas[i], shapes[j], scales[k], num_stars; minP=minP, maxP=maxP, num_evals=num_evals )
      end # for i
    end # for j
  end # for k
  dist
end # function

