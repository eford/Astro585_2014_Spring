# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

srand(42);
N = 10000;
true_mean = 10000.;
y = true_mean+randn(N);
sample_mean = mean(y);
sample_var = var(y); 
sample_mean , sample_var

# <codecell>

y32bit = convert(Array{Float32,1},y);
y16bit = convert(Array{Float16,1},y);
sample_mean32bit = mean(y32bit);
sample_mean16bit = mean(y16bit);
sample_var32bit = var(y32bit);
sample_var16bit = var(y16bit);
(sample_mean32bit - sample_mean, sample_mean16bit - sample_mean)

# <codecell>

(sample_var32bit -sample_var, sample_var16bit -sample_var)

# <codecell>

srand(42);
N = 10000;
true_mean = 100000.;
y = true_mean+randn(N);
sample_mean = mean(y);
sample_var = var(y); 
sample_mean , sample_var

# <codecell>

y32bit = convert(Array{Float32,1},y);
y16bit = convert(Array{Float16,1},y);
sample_mean32bit = mean(y32bit);
sample_mean16bit = mean(y16bit);
sample_var32bit = var(y32bit);
sample_var16bit = var(y16bit);
sample_mean32bit - sample_mean, sample_mean16bit - sample_mean

# <codecell>

sample_var32bit -sample_var, sample_var16bit -sample_var

# <codecell>

function var_1pass(y::Array)
   n = length(y);
   sum = zero(y[1]);
   sumsq = zero(y[1]);
   for i in 1:n
      sum = sum + y[i];
	  sumsq = sumsq + y[i]*y[i];
   end
   variance = (sumsq - (sum*sum)/n)/(n-1);
end

# <codecell>

function var_2pass(y::Array)
   n = length(y);
   sum1 = zero(y[1]);
   for i in 1:n
      sum1 = sum1 + y[i];
   end   
   mean = sum1 / n;
   sum2 = zero(y[1]);
   for i in 1:n
	  sum2 = sum2 + (y[i]-mean)*(y[i]-mean);
   end   
   variance = sum2/(n-1);
end

# <codecell>

srand(42);
N = 100000;
true_mean = 100000.;
y = true_mean+randn(N);
var_1pass(y)-var(y), var_2pass(y)-var(y) 

# <codecell>

function var_online(y::Array)
  n = length(y);
  sum1 = zero(y[1]);
  mean = zero(y[1]);
  M2 = zero(y[1]);
  for i in 1:n
	  diff_by_i = (y[i]-mean)/i;
	  mean = mean +diff_by_i;
	  M2 = M2 + (i-1)*diff_by_i*diff_by_i+(y[i]-mean)*(y[i]-mean); 
  end;  
  variance = M2/(n-1);
end
var_online(y)-var(y)

# <codecell>

	srand(42);
	Nobs = 100;
	z = zeros(Nobs);
	sigma = 2. * ones(Nobs);
	y = z + sigma .* randn(Nobs);

# <codecell>

normal_pdf(z, y, sigma) =	exp(-0.5.*((y.-z)./sigma).^2)./(sqrt(2*pi).*sigma);

# <codecell>

function likelihood(y::Array, sigma::Array, z::Array)
   n = length(y);
   assert(length(sigma)==n);    
   assert(length(z)==n);
   prod = one(y[1]);
   for i in 1:n
      prod = prod * normal_pdf(y[i],z[i],sigma[i]);
   end;
   return prod;
end

# <codecell>

@time likelihood(y,sigma,z)

# <codecell>

likelihood_prod(y::Array, sigma::Array, z::Array) = prod(normal_pdf(y,z,sigma));

# <codecell>

likelihood_prod(y,sigma,z) 

# <codecell>

likelihood_reduce_map(y::Array, sigma::Array, z::Array) = reduce(*,map(normal_pdf,y,z,sigma));

# <codecell>

likelihood_reduce_map(y,sigma,z) 

# <codecell>

likelihood_mapreduce(y::Array, sigma::Array, z::Array) = mapreduce(i->normal_pdf(y[i],z[i],sigma[i]),*,1:length(y));

# <codecell>

likelihood_mapreduce(y,sigma,z) 

# <codecell>

srand(42);
Nobs = 600;
z = zeros(Nobs);
sigma = 2. * ones(Nobs);
y = z + sigma .* randn(Nobs);

# <codecell>

likelihood(y,sigma,z) 

# <codecell>

function log_likelihood(y::Array, sigma::Array, z::Array)
   n = length(y);
   assert(length(sigma)==n);
   assert(length(z)==n);
   ll= zero(y[1]);
   for i in 1:n
      ll = ll + log(normal_pdf(z[i],y[i],sigma[i]));
   end;
   return ll;
end

# <codecell>

log_likelihood(y,sigma,z)

# <codecell>

log_normal_pdf(z, y, sigma) = -0.5*( ((y-z)/sigma)^2 + log(2*pi*sigma*sigma));

# <codecell>

function log_likelihood_2(y::Array, sigma::Array, z::Array)
   n = length(y);
   assert(length(sigma)==n);
   assert(length(z)==n);
   ll= zero(y[1]);
   for i in 1:n
      ll = ll + log_normal_pdf(z[i],y[i],sigma[i]);
   end;
   return ll;
end

# <codecell>

log_likelihood_2(y,sigma,z)

# <codecell>

function round_down_to_power_of_ten(x::Real)
  if(!isfinite(x))
    if(isnan(x)) return nan(x) end
    if(isinf(x)) return x end
  end
  assert(x>0.);
  z = 1.0;
  if x >= 1.0;
    while z*10.<=x
      z = z * 10.0;
    end
  else
    while z > x # correct 
     z = z / 10.0;
    end
  end
 return z;
end

# <codecell>

round_down_to_power_of_ten(0.1)

# <codecell>

round_down_to_power_of_ten(-0.1)

# <codecell>

round_down_to_power_of_ten(NaN)

# <codecell>

round_down_to_power_of_ten(x::Real) = 10.0^(floor(log10(x)));

# <codecell>


