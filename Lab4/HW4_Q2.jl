# integrand that we'll use for testing
std_norm_pdf(x) = exp(-0.5*x.*x)./sqrt(2pi) # standard normal pdf

# Calculate \int_a^b dx func(x) using n function evaluations
# Approximates integral as uniformly spaced rectangles
# Avoids evaluating func at a or b, in case of singularities
# Uses a for loop to evaluate func
function int_loop(func::Function, a::Real, b::Real, n::Integer = 1000000)
  @assert(n>2)
  integral = 0.
  for i in 1:n
    x = a+i*(b-a)/(n+1)
    eval = func(x)
    integral += eval
  end
  integral *= (b-a)/n;
end

# Calculate \int_a^b dx func(x) using n function evaluations
# Approximates integral as uniformly spaced rectangles
# Avoids evaluating func at a or b, in case of singularities
# Uses vector notation to evaluate func
function int_vec(func::Function, a::Real,b::Real, n::Integer = 1000000)
  @assert(n>2)
  idx = [1:n]
  x = a+idx.*(b-a)./(n+1) 
  integral = sum(func(x)) * (b-a)/n 
end

# Calculate \int_a^b dx func(x) using n function evaluations
# Approximates integral as uniformly spaced rectangles
# Avoids evaluating func at a or b, in case of singularities
# Uses vector notation to evaluate func and 
#     a comprehension to construct function arguements
function int_vec2(func::Function, a::Real,b::Real, n::Integer = 1000000)
  @assert(n>2)
  x = [ a+i*(b-a)/(n+1) for i in 1:n ]
  integral = sum(func(x)) * (b-a)/n 
end

# Calculate \int_a^b dx func(x) using n function evaluations
# Approximates integral as uniformly spaced rectangles
# Avoids evaluating func at a or b, in case of singularities
# Uses map to evaluate func 
function int_map(func::Function, a::Real,b::Real, n::Integer = 1000000)
  @assert(n>2)
  x = [ a+i*(b-a)/(n+1) for i in 1:n ]
  integral = sum(map(func,x)) * (b-a)/n 
end

# Calculate \int_a^b dx func(x) using n function evaluations
# Approximates integral as uniformly spaced rectangles
# Avoids evaluating func at a or b, in case of singularities
# Uses pmap to evaluate func 
function int_pmap(func::Function, a::Real,b::Real, n::Integer = 1000000)
  @assert(n>2)
  x = [ a+i*(b-a)/(n+1) for i in 1:n ]
  integral = sum(pmap(func,x)) * (b-a)/n 
end

# Calculate \int_a^b dx func(x) using n function evaluations
# Approximates integral as uniformly spaced rectangles
# Avoids evaluating func at a or b, in case of singularities
# Uses mapreduce to evaluate sum of func 
function int_mapreduce(func::Function, a::Real,b::Real, n::Integer = 1000000)
  @assert(n>2)
  x = [ a+i*(b-a)/(n+1) for i in 1:n ]
  integral = mapreduce(func,+,x)
  integral *= (b-a)/n 
end

# Calculate \int_a^b dx func(x) using n function evaluations
# Approximates integral as uniformly spaced rectangles
# Avoids evaluating func at a or b, in case of singularities
# Uses vector notation to evaluate the sum of func and
#     the Devectorize package to minimize temporaries
using Devectorize
function int_devec_std_normal(a::Real,b::Real, n::Integer = 1000000)
  @assert(n>2)
  x = [ a+i*(b-a)/n for i in 1:n ]
  @devec integral = sum(exp(-0.5.*x.*x)./sqrt(2.*pi)) .* ((b-a)./n) # Write it out, rather than calling func so devectorize will do it's magic
end


# Calculate \int_a^b dx standard_normal_pdf(x) using n function evaluations
# Approximates integral as uniformly spaced rectangles
# Avoids evaluating func at a or b, in case of singularities
# Uses a hard wired for loop
function int_loop2_std_normal(a::Real,b::Real, n::Integer = 1000000)
  @assert(n>2)
  integral = 0.
  invsqrttwopi = 1.0/sqrt(2.*pi);
  x = [ a+i*(b-a)/n for i in 1:n ]
  for z in x
    integral += exp(-0.5*z*z)*invsqrttwopi
  end
  integral = integral *(b-a)/(n);
end

# Calculate \int_a^b dx standard_normal_pdf(x) using n function evaluations
# Approximates integral as uniformly spaced rectangles
# Avoids evaluating func at a or b, in case of singularities
# Uses a hard wired for loop and a temporary array of function arguements
function int_loop3_std_normal(a::Real,b::Real, n::Integer = 1000000)
  @assert(n>2)
  integral = 0.
  invsqrttwopi = 1.0/sqrt(2.*pi)
  for i in 1:n
    x = a+i*(b-a)/n 
    integral += exp(-0.5*x*x)*invsqrttwopi
  end
  integral = integral* (b-a)/(n);
end

# Calculate \int_a^b dx standard_normal_pdf(x) using n function evaluations
# Approximates integral as uniformly spaced rectangles
# Avoids evaluating func at a or b, in case of singularities
# Uses a hard wired for loop and a temporary array of function arguements
function int_loop4_std_normal(a::Real,b::Real, n::Integer = 1000000)
  @assert(n>2)
  integral = 0.
  invsqrttwopi = 1.0/sqrt(2.*pi)
  integral = @parallel (+) for i in 1:n
    x = a+i*(b-a)/n 
    exp(-0.5*x*x)*invsqrttwopi
  end
  integral = integral* (b-a)/(n);
end



# Tests for accuracy
using Base.Test
function test_int_std_normal(f::Function, N::Integer = 1000000, tol::Real = 0.0001)
  for i in 0:5
    param = float64(i);
	@test_approx_eq_eps(f(std_norm_pdf,-param,param,N),erf(param/sqrt(2.0)), tol )
  end 
end

function test_int_hardwired_std_normal(f::Function, N::Integer = 1000000, tol::Real = 0.0001)
  for i in 0:5
    param = float64(i);  @test_approx_eq_eps(f(-param,param,N),erf(param/sqrt(2.0)), tol )
  end 
end

function test_all()
  test_int_hardwired_std_normal(int_loop4_std_normal)
  test_int_std_normal(int_loop)
  test_int_std_normal(int_vec)
  test_int_std_normal(int_vec2)
  test_int_std_normal(int_map)
  test_int_std_normal(int_pmap)
  test_int_std_normal(int_mapreduce)
  test_int_hardwired_std_normal(int_devec_std_normal)
  test_int_hardwired_std_normal(int_loop2_std_normal)
  test_int_hardwired_std_normal(int_loop3_std_normal)
end 

# Benchmark each implementation of integration routine
function benchmark_all(Nevals::Integer = 1000000, Nrepeats::Integer = 10 )
  println(@elapsed { int_loop(std_norm_pdf,-1.0,1.0,Nevals) for i in 1:Nrepeats } )
  println(@elapsed { int_vec(std_norm_pdf,-1.0,1.0,Nevals) for i in 1:Nrepeats } )
  println(@elapsed { int_vec2(std_norm_pdf,-1.0,1.0,Nevals) for i in 1:Nrepeats } )
  println(@elapsed { int_map(std_norm_pdf,-1.0,1.0,Nevals) for i in 1:Nrepeats } )
  println(@elapsed { int_pmap(std_norm_pdf,-1.0,1.0,Nevals) for i in 1:Nrepeats } )
  println(@elapsed { int_mapreduce(std_norm_pdf,-1.0,1.0,Nevals) for i in 1:Nrepeats } )
  println(@elapsed { int_devec_std_normal(-1.0,1.0,Nevals) for i in 1:Nrepeats } )
  println(@elapsed { int_loop3_std_normal(-1.0,1.0,Nevals) for i in 1:Nrepeats } )
  println(@elapsed { int_loop4_std_normal(-1.0,1.0,Nevals) for i in 1:Nrepeats } )
end

