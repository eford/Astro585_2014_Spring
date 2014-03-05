@everywhere include("HW4_Q2.jl")

# Run Tests
test_all()

println("Benchmarking each implementation...")
benchmark_all(10000000,10)

if(0>1)
addprocs(2)
@everywhere include("HW_Q2.jl")
println("Benchmarking w/ parallel turn on for each implementation...")
benchmark_all(10000000,10)
rmprocs(2)
end

# Profile
if(0>1)
println("\n\nProfiling int_vec2")
Profile.clear()
@profile { int_vec2(std_norm_pdf,-1.0,1.0,1000000) for i in 1:10 };
Profile.print()

println("Profiling int_devec_std_normal")
Profile.clear()
@profile { int_devec_std_normal(-1.0,1.0,1000000) for i in 1:10 };
Profile.print()
end