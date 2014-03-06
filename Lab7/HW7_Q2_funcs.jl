function int_normal_gpu(a::Real, b::Real, n::Integer = 1000000)

@assert(typeof(a)==typeof(b)==Float64)

# generate random arrays 
dx = (b-a)/n
x_cpu = linspace( convert(Float64,a+0.5*dx), convert(Float64,b-0.5*dx), n) 

println("Timing the upload of inputs from CPU to GPU and memory allocation on GPU.")
tic()
# load input array onto GPU
x_gpu = CuArray( x_cpu )

# create an array on GPU to store results
y_gpu = CuArray(Float64, (n))
toc()

# choose the block and grid sizes based on the problem size
block_size = choose_block_size(n)
grid_size = choose_grid_size(n,block_size)

stream1 = CUDA.null_stream()

println("Timing the gpu calculations.")
tic()
# run the kernel 
# syntax: launch(kernel, grid_size, block_size, arguments)
# here, grid_size and block_size could be an integer or a tuple of integers
launch(normal_pdf_kernel, grid_size, block_size, (x_gpu, y_gpu, n), stream=stream1)
synchronize(stream1)
toc()

println("Timing the download of results from GPU to CPU")
tic()
# download the results from GPU
y_cpu = to_host(y_gpu)   # c is a Julia array on CPU (host)
toc()

# release GPU memory
free(x_gpu)
free(y_gpu)

println("Timing the summation being performed on the CPU")
tic()
result = sum(y_cpu) * (b-a)/n
toc()

return result
end



