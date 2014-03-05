
@everywhere normal_pdf(x) = exp(-0.5*x.*x)./sqrt(2.*pi)
@everywhere using Base.Test

@everywhere function test_int_normal_pdf(func::Function; n::Integer = 1000000, eps::Real = 1.0e-6)
  limits = 1:5
  for limit in limits
    @test_approx_eq_eps func(-limit,limit) erf(limit/sqrt(2.0)) eps
  end
end

@everywhere function int_normal_pdf_loop(a::Real, b::Real, n::Integer = 1000000)
  @assert(n>2)
  dx = (b-a)/n  
  x = (a+0.5*dx):dx:(b-0.5*dx)  
  integral = @parallel (+) for i in 1:n
    normal_pdf(x[i])
  end
  integral *= (b-a)/n;
end

function int_normal_pdf_map(a::Real, b::Real, n::Integer = 1000000)
  @assert(n>2)
  x = distribute([ a+i*(b-a)/(n+1) for i in 1:n ])
  integral = (b-a)/n * reduce(+,map(fetch,{ (@spawnat p sum(normal_pdf(localpart(x)))) for p in procs(x) }))
end

function time_int_normal_pdf_functions(n::Integer = 1000000)
  # make sure functions compiled
  int_normal_pdf_loop(-1.0,1.0,100)
  int_normal_pdf_map(-1.0,1.0,100)
  timeloop = @elapsed  int_normal_pdf_loop(-1.0,1.0,n)
  timemap = @elapsed  int_normal_pdf_map(-1.0,1.0,n)
  (timeloop, timemap)
end

