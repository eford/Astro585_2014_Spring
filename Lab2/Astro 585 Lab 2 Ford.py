# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

normal_pdf(z, y, sigma) = exp(-0.5*((y-z)/sigma)^2)/(sqrt(2*pi)*sigma);
function log_likelihood(y::Array, sigma::Array, z::Array)
   n = length(y);
   sum  = zero(y[1]);
   for i in 1:n
    sum = sum + log(normal_pdf(y[i],z[i],sigma[i]));
   end;
   return sum;
end

# <codecell>

function calc_time_log_likelihood(Nobs::Int, Mtimes::Int = 1)
  @assert (Nobs>=1);
  srand(42);
  z = zeros(Nobs);
  sigma = 2. * ones(Nobs);
  y = z + sigma .* randn(Nobs);
  total = 0.;
  for i in 1:Mtimes
    total = total + @elapsed log_likelihood(y,sigma,z);
  end
  return total;
end

# <codecell>

calc_time_log_likelihood(100,1)

# <codecell>

calc_time_log_likelihood(100,1)

# <codecell>

using PyPlot
n_list = [ 2^i for i=1:10 ]
elapsed_list = map(calc_time_log_likelihood,n_list)
plot(log10(n_list), log10(elapsed_list), color="red", linewidth=2, marker="+", markersize=12);
xlabel("log N"); ylabel("log (Time/s)"); 

# <codecell>

calc_time_log_likelihood(100)

# <codecell>

function log_likelihood(y::Array, sigma::Array, z::Array)
   n = length(y);
   assert(length(sigma)==n);    
   assert(length(z)==n);
   sum = zero(y[1]);
   for i in 1:n
    sum = sum + log(normal_pdf(y[i],z[i],sigma[i]));
   end;
   return sum;
end
calc_time_log_likelihood(100,100);
calc_time_log_likelihood(100,100)

# <codecell>

function log_likelihood(y::Array{Float64}, sigma::Array{Float64}, z::Array{Float64})
   n = length(y);
   assert(length(sigma)==n);    
   assert(length(z)==n);
  sum = zero(y[1]);
   for i in 1:n
    sum = sum + log(normal_pdf(y[i],z[i],sigma[i]));
   end;
  return sum;
end
calc_time_log_likelihood(100,100);
calc_time_log_likelihood(100,100)

# <codecell>

function log_likelihood_generic(y::Array{Float64}, sigma::Array{Float64}, z::Array{Float64}, pdf::Function)
   n = length(y);
   assert(length(sigma)==n);    
   assert(length(z)==n);
  sum = zero(y[1]);
   for i in 1:n
    sum = sum + log(pdf(y[i],z[i],sigma[i]));
   end;
  return sum;
end
function calc_time_log_likelihood(Nobs::Int, pdf::Function, Mtimes::Int = 1)
  @assert (Nobs>=1);
  srand(42);
  z = zeros(Nobs);
  sigma = 2. * ones(Nobs);
  y = z + sigma .* randn(Nobs);
  total = 0.;
  total = sum( [@elapsed log_likelihood_generic(y,sigma,z,pdf) for i in 1:Mtimes] );
end
calc_time_log_likelihood(100,normal_pdf,100);
calc_time_log_likelihood(100,normal_pdf,100)

# <codecell>

function log_likelihood(y::Array{Float64}, sigma::Array{Float64}, z::Array{Float64})
   n = length(y);
   assert(length(sigma)==n);    
   assert(length(z)==n);
   sum = zero(y[1]);
   for i in 1:n
    sum = sum -0.5*((y[i]-z[i])/sigma[i])^2 - log(sigma[i]);
   end;
  sum = sum - 0.5*n*log(2*pi);
  return sum;
end
calc_time_log_likelihood(100,100);
calc_time_log_likelihood(100,100)

# <codecell>

function log_likelihood(y::Array{Float64}, sigma::Array{Float64}, z::Array{Float64})
   n = length(y);
   assert(length(sigma)==n);    
   assert(length(z)==n);
   ll = sum(((y.-z)./sigma).^2 .- log(sigma));
   ll = - 0.5*(ll+n*log(2*pi));
   return ll;
end
calc_time_log_likelihood(100,100);
calc_time_log_likelihood(100,100)

# <codecell>

Pkg.add("Devectorize")
using Devectorize;

# <codecell>

function log_likelihood(y::Array{Float64}, sigma::Array{Float64}, z::Array{Float64})
   n = length(y);
   assert(length(sigma)==n);    
   assert(length(z)==n);
   @devec ll = sum(((y.-z)./sigma).^2 .- log(sigma));
   ll = - 0.5*(ll+n*log(2*pi));
   return ll;
end
calc_time_log_likelihood(100,100);
calc_time_log_likelihood(100,100)

# <codecell>

N = 1000000;
println("rand:  ", 1./(@elapsed x = rand(N)));
println("randn: ", 1./(@elapsed y = randn(N)));
println(".*:    ", 1./(@elapsed x.*x));
println(".*:    ", 1./(@elapsed x.*y));
println(".+:    ",1./(@elapsed x.+y));
println("sqrt:  ", 1./(@elapsed sqrt(x)));
println("sin:  ", 1./(@elapsed sqrt(x)));
println("log:   ", 1./(@elapsed log(x)));

# <codecell>

function update_derivs!(state::Array{Float64}, derivs::Array{Float64} )
  @assert length(state) == 4
  @assert length(derivs) == 4
  GM = 1.0
  r_x = state[1]
  r_y = state[2]
  v_x = state[3]
  v_y = state[4]
  r2 = r_x*r_x+r_y*r_y
  a = -GM/r2
  r = sqrt(r2)
  a_x = a * r_x/r
  a_y = a * r_y/r
  derivs[1] = v_x
  derivs[2] = v_y
  derivs[3] = a_x
  derivs[4] = a_y
  return derivs;
end
function advance_euler!(state::Array{Float64,1},derivs::Array{Float64,1}, dt::Float64)
  update_derivs!(state,derivs);
  {state[i] = state[i] + dt*derivs[i]  for i in 1:4}
end
advance!(state::Array{Float64,1},derivs::Array{Float64,1}, dt::Float64) = advance_euler!(state,derivs,dt)

# <codecell>

function integrate!(state::Array{Float64,1},dt::Float64, duration::Float64)
  nsteps = iceil(duration/dt);
  log = Array(Float64,(nsteps,length(state)));
  n = 0
  t = 0.0
  derivs = Array(Float64,4);
  while t<duration
    dt_tmp = (t+dt<=duration) ? dt : duration-t;
    advance!(state,derivs,dt_tmp);
    t = t + dt_tmp
    n = n + 1
    #@assert(n<=nsteps)
    if n<=nsteps log[n,:] = deepcopy(state) end
  end
  return log
end

# <codecell>

state = [1.,0.,0.,1.];  
#euler_log = integrate_euler!(state,2pi/200,3.0*2pi);
@time euler_log = integrate!(state,2pi/200,3.0*2pi);

# <codecell>

using PyPlot;
plot(euler_log[:,1],euler_log[:,2], color="red", linewidth=2, linestyle="-")
xlabel("x"); ylabel("y");  title("Trajectory w/ Euler Integrator");

# <codecell>

function update_derivs_pos!(state::Array{Float64}, derivs::Array{Float64} )
  @assert length(state) == 4
  @assert length(derivs) == 4
  v_x = state[3]
  v_y = state[4]
  derivs[1] = v_x
  derivs[2] = v_y
  return derivs;
end
function update_derivs_vel!(state::Array{Float64}, derivs::Array{Float64} )
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
function advance_leapfrog!(state::Array{Float64,1},derivs::Array{Float64,1}, dt::Float64)
  update_derivs_pos!(state,derivs);
  {state[i] = state[i] + 0.5*dt*derivs[i]  for i in 1:2}
  update_derivs_vel!(state,derivs);
  {state[i] = state[i] + dt*derivs[i]      for i in 3:4}
  update_derivs_pos!(state,derivs);
  {state[i] = state[i] + 0.5*dt*derivs[i]  for i in 1:2}    
end

# <codecell>

state = [1.,0.,0.,1.];  
advance!(state::Array{Float64,1},derivs::Array{Float64,1}, dt::Float64) = advance_leapfrog!(state,derivs,dt)
# Need to evaluate integrate again, so it knows to use the new advance is 

# <codecell>

@time leapfrog_log  = integrate!(state,2pi/200,3.0*2pi);

# <codecell>

using PyPlot;
plot(leapfrog_log[:,1],leapfrog_log[:,2], color="red", linewidth=2, linestyle="-")
xlabel("x"); ylabel("y"); title("Trajectory w/ Leapfrog Integrator");

# <codecell>

time10000 = @elapsed integrate!(state,2pi/200.,10000*2pi);
println(time10000, "s -> ", 4.5e9/10000*time10000/(24.*60.*60.), "d");

# <codecell>

function int_end_distance(dur::Float64, dt::Float64 = 2pi/200.0)
  state = [1.,0.,0.,1.];  
  integrate!(state,dt,dur*2pi);
  dist = state[1]^2+state[2]^2
  return dist-1.0;
end
int_end_distance(3.,2pi/100.)/int_end_distance(3.,2pi/200.)

# <codecell>

inlist= [ 2.0^i for i=1:10 ]
outlist = map(int_end_distance,inlist)
int_end_dist2(dur) = int_end_distance(dur,2pi/100.0)
outlist2 = map(int_end_dist2,inlist)

# <codecell>

using PyPlot
plot(log10(inlist),log10(abs(outlist)), log10(inlist),log10(abs(outlist2)))
xlabel("log Duration"); ylabel("log | |x|_final - 1 |");  title("Accuracy vs Integration Duration");

# <codecell>

log_likelihood_reduce_map(y::Array, sigma::Array, z::Array) = reduce(+,map(log,map(normal_pdf,y,z,sigma)));
@elapsed { log_likelihood_reduce_map(y,sigma,z) for i in 1:10000 } 

# <codecell>

likelihood_mapreduce(y::Array, sigma::Array, z::Array) = mapreduce(normal_pdf,*,y,z,sigma);

# <codecell>

@elapsed { likelihood_mapreduce(y,sigma,z) for i in 1:10000 } 

# <codecell>

using NumericExtensions

# <codecell>

function likelihood_nofunc(y::Array, sigma::Array, z::Array)
   n = length(y);
   assert(length(sigma)==n);    
   assert(length(z)==n);
   prod = one(y[1]);
   for i in 1:n
    prod = prod * exp(-0.5*((y[i]-z[i])/sigma[i])^2)/(sqrt(2*pi)*sigma[i]);
   end;
   return prod;
end

# <codecell>

@elapsed { likelihood_nofunc(y,sigma,z) for i in 1:10000 } 

# <codecell>

N=3; M=4;
datalog = Array(Float64,(N,M))

# <codecell>

datalog[1:end,2]

# <codecell>

datalog = zeros(N,M)

# <codecell>

push!(datalog[:,1],2.)

# <codecell>


