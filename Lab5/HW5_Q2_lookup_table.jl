# abstract base type for methods to approximate a function via a lookup table
abstract lookup_table_type

# derived type for lookup_table using linear interpolation
type lookup_table_linear_type <:  lookup_table_type
  y::Vector{Float64}
  dy::Vector{Float64}
  a::Real	
  b::Real
  n::Integer
end

# derived type for lookup_table using quadratic interpolation
type lookup_table_quadratic_type <:  lookup_table_type
  y::Vector{Float64}
  dy::Vector{Float64}
  d2y::Vector{Float64}                  
  a::Real
  b::Real
  n::Integer
end

# check if the required common elements (implied "base") of a lookup_table are in a valid state
function is_valid_base(table::lookup_table_type)
  if(table.n<=0) return false end
  if(table.b<=table.a) return false end
  return true
end

# check if table is in a valid state
function is_valid(table::lookup_table_linear_type)
  if !is_valid_base(table) return false end
  if(length(table.y) != length(table.dy)) return false end
  if(length(table.y)<2) return false end
  return true
end

# check if table is in a valid state
function is_valid(table::lookup_table_quadratic_type)
  if !is_valid_base(table) return false end
  if(length(table.y) != length(table.dy)) return false end
  if(length(table.d2y) != length(table.dy)) return false end
  if(length(table.y)<2) return false end
  return true
end

# check if x is in range for interpolation 
isinrange(t::lookup_table_type, x::Real) = (t.a<= x <=t.b)

# Construct a lookup_table using function f and it's derivative g, over range [a,b] with n evaluations
function make_table_linear(f::Function,g::Function,a::Real,b::Real,n::Integer = 10)
  @assert(a<b)
  @assert(n>=2)
  x = linspace(a,b,n)
  y = similar(x)
  dy = similar(x)
  for i in 1:length(x)
      y[i] = f(x[i])
      dy[i] = g(x[i])
      # TODO: Might consider adding code to check that g is close to numerical derivative of f  
  end
  return lookup_table_linear_type(y,dy,a,b,n)
end

using Calculus
# Construct a lookup_table using function f and it's derivative computed numerically, over range [a,b] with n evaluations
function make_table_linear(f::Function,a::Real,b::Real,n::Integer = 10)
  @assert(a<b)
  @assert(n>=2)
  x = linspace(a,b,n)
  y = similar(x)
  dy = similar(x)
  for i in 1:length(x)
      y[i] = f(x[i])
      dy[i] = derivative(f, x[i])
  end
  return lookup_table_linear_type(y,dy,a,b,n)
end

# Construct a lookup_table using function f and it's derivative computed numerically, over range [a,b] with n evaluations
function make_table_quadratic(f::Function,a::Real,b::Real,n::Integer = 10)
  @assert(a<b)
  @assert(n>=2)
  x = linspace(a,b,n)
  y = similar(x)
  dy = similar(x)
  d2y = similar(x)
  for i in 1:length(x)
      y[i] = f(x[i])
      dy[i] = derivative(f, x[i])
      d2y[i] = second_derivative(f, x[i])
  end
  return lookup_table_quadratic_type(y,dy,d2y,a,b,n)
end

# approximate f(x) using linear lookup table
function lookup( table::lookup_table_linear_type, x::Real )
   @assert(is_valid(table))
   @assert(isinrange(table,x))
   # compute index of function evaluation closest to x
   idx = iround(1+(x-table.a)/(table.b-table.a)*(table.n-1))
   # calculate location of function evaluation closest to x
   xo = table.a+(table.b-table.a)*(idx-1)/(table.n-1)
   dx = x - xo
   # perform linear interpolation
   table.y[idx] + table.dy[idx] * dx
end

# approximate f(x) using quadratic lookup table
# Warning:  This uses the function and its derivatives at one point, so it does not result in a continuous function
function lookup( table::lookup_table_quadratic_type, x::Real )
   @assert(is_valid(table))
   @assert(isinrange(table,x))
   # compute index of function evaluation closest to x
   idx = iround(1+(x-table.a)/(table.b-table.a)*(table.n-1))
   # calculate location of function evaluation closest to x
   xo = table.a+(table.b-table.a)*(idx-1)/(table.n-1)
   dx = x - xo
   # perform quadratic interpolation
   table.y[idx] + dx * ( table.dy[idx] + dx* 0.5*table.d2y[idx])
end

# approximate f(x) for an array of x's using table and write output into output array
function lookup!( table::lookup_table_type, x::Array{Float64}, output::Array{Float64} )
   @assert(size(output) == size(x) )
   for i in 1:length(x)
       output[i] = lookup(table,x[i])
   end
   output
end

# approximate f(x) using table using table 
function lookup( table::lookup_table_type, x::Array{Float64} )
   output = similar(x)
   for i in 1:length(x)
       output[i] = lookup(table,x[i])
   end
   output
end

