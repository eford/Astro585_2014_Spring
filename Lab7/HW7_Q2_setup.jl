# Pkg.add("CUDA")
using CUDA

function init_cuda(gpuid::Integer = 0)
  # select a CUDA device
  dev = CuDevice(gpuid)

  # create a context (like a process in CPU) on the selected device
  ctx = create_context(dev)

  # Technically we should query GPU to get these parameters, but these should work with GPUs on sagan and tesla
  const global max_grid_dim_size = 65535
  const global max_block_dim_size = 512

  return ctx
end

function load_functions()
  # load the PTX module (each module can contain multiple kernel functions)
  md = CuModule("normal_pdf_gpu.ptx")

  # retrieve the kernel functions from the module
  global normal_pdf_kernel = CuFunction(md, "normal_pdf_gpu")

  return md
end
 
function unload_functions(md::CuModule)
  # finalize: unload module and destroy context
  unload(md)
end

function close_cuda(ctx::CuContext)
  destroy(ctx)
end

#  Could optimize the chocie of block sizes better
function choose_block_size(n::Integer)
  min(n,max_block_dim_size)
end

#  Could optimize the chocie of grid sizes better
function choose_grid_size(n::Integer, block_size::Integer)
  @assert(block_size <= max_block_dim_size)

  grid_size = iceil(n/block_size)
  if grid_size > max_grid_dim_size
    grid_size_x = iceil(sqrt(grid_size))
    @assert(grid_size_x <= max_grid_dim_size)
    grid_size = (grid_size_x, grid_size_x)
  end
  return grid_size
end

