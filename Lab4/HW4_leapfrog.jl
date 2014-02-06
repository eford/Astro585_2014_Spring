# HW4_Q3_leapfrog.jl

function update_derivs_pos!(state::Vector{Float64}, derivs::Vector{Float64} )
 # Input: state = [x,y,vx,vy], a vector of two 2-d positions and velocities for a test particle
 # Output: The derivatives of the position are updated in the preallocated array derivs.
  @assert length(state) == 4
  @assert length(derivs) == 4
  v_x = state[3]
  v_y = state[4]
  derivs[1] = v_x
  derivs[2] = v_y
  return derivs;
end

function update_derivs_vel!(state::Vector{Float64}, derivs::Vector{Float64} )
  # Input: state = [x,y,vx,vy], a vector of two 2-d positions and velocities for a test particle
  # Output: The derivatives of the velocity are updated in the preallocated array derivs.  
  @assert length(state) == 4
  @assert length(derivs) == 4
  GM = 1.0
  r_x = state[1]
  r_y = state[2]
  r2 = r_x*r_x+r_y*r_y
  a = -GM/r2
  r = sqrt(r2)
  a_x = a * r_x/r
  a_y = a * r_y/r
  derivs[3] = a_x
  derivs[4] = a_y
  return derivs;
end

function update_derivs!(state::Vector{Float64}, derivs::Vector{Float64} )
  # Input: state = [x,y,vx,vy], a vector of two 2-d positions and velocities for a test particle
  # Output: The derivatives are updated in the preallocated array derivs.  
  update_derivs_vel!(state,derivs)
  update_derivs_pos!(state,derivs)
  return derivs
end
  
function advance_leapfrog!(state::Vector{Float64},derivs::Vector{Float64}, dt::Float64; derivs_current::Bool = false)
  # Input/Output: state = array of two 2-d positions and velocities for a test particle
  # Temporary space: The derivatives are updated in the preallocated array derivs.
  # Input: dt is the fixed time step 
  # Optional param: derivs_current: whether need to calculate derivatives at beginning
  @assert length(state) == length(derivs)
  
  if !derivs_current 
    update_derivs_pos!(state,derivs);
  end
  {state[i] = state[i] + 0.5*dt*derivs[i]  for i in 1:2}
  update_derivs_vel!(state,derivs);
  {state[i] = state[i] + dt*derivs[i]      for i in 3:4}
  update_derivs_pos!(state,derivs);
  {state[i] = state[i] + 0.5*dt*derivs[i]  for i in 1:2}    
end

# Input/Output: state = [x,y,vx,vy], an array of two 2-d positions and velocities for a test particle
# Input: dt is the fixed time step 
# Input: duration is the total 
function integrate_leapfrog!(state::Vector{Float64}, dt::Float64, duration::Float64; max_num_log::Integer = 100000)
  @assert(length(state)==4) 
  @assert(dt>0.0)
  @assert(duration>0.0)
  
  # Preallocate array to hold data log (including  initial state)
  nsteps = iceil(duration/dt);
  nskip = (nsteps<max_num_log) ? 1 : iceil(nsteps/(max_num_log-1))
  num_log = iceil(nsteps/nskip)+1   
  log = Array(Float64,(num_log,length(state)));

  # Pre-allocate and pre-compute derivaties
  derivs = Array(Float64,4);  
  update_derivs!(state,derivs);

  # Log initial state 
  log_pos = 1
  log[log_pos,:] = deepcopy(state) 

  n = 0
  t = 0.0
  while t<duration
    # ensure don't integrate for more than duration
    dt_tmp = (t+dt<=duration) ? dt : duration-t;

	# advance system by one time step
    advance_leapfrog!(state,derivs,dt_tmp, derivs_current=true)
    t = t + dt_tmp
    n = n + 1

    if (n%nskip==0) # Log data
	   log_pos += 1
   	   @assert( log_pos<=length(log) )
	   @assert( length(log[log_pos,:])==length(state) )
	   log[log_pos,:] = deepcopy(state) 
	end
  end
  return log
end

function calc_error_leapfrog_old(dur::Float64, dt::Float64 = 2pi/200.0)
  state = [1.,0.,0.,1.];  
  integrate_leapfrog!(state,dt,dur*2pi);
  dist = state[1]^2+state[2]^2
  phase = atan2(state[2],state[1])
  offset = sum((state[1:2].-[1.0,0.0]).^2)
  return (dist-1.0,phase,offset)
end

function calc_end_distance_leapfrog(dur::Integer, dt::Float64 = 2pi/200.0, state::Vector{Float64} = [1., 0., 0., 1.] )
  #state = [1.,0.,0.,1.];  
  dist_init = sqrt(state[1]^2+state[2]^2)
  phase_init = atan2(state[2],state[1])
  integrate_leapfrog!(state,dt,dur*2pi);
  # Calculate three metrics of the accuracy of the integration  
  dist = sqrt(state[1]^2+state[2]^2)
  phase = atan2(state[2],state[1])
  offset = sqrt(sum((state[1:2].-[1.0,0.0]).^2))
  return (dist-dist_init,phase-phase_init,offset)
end

using Base.Test
function test_leapfrog(dur::Integer, dt::Float64 = 2pi/200.0)  
  err = calc_end_distance_leapfrog(dur,dt)
  @test_approx_eq_eps(err[1], 0., 1e-6)
  @test_approx_eq_eps(err[2], 0., 1e-1)
  @test_approx_eq_eps(err[3], 0., 1e-1)
end

