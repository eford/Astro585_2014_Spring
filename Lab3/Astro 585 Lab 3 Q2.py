# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

N = 1024
x = rand(N);

# <codecell>

dt_sprintf_stdout = @elapsed for i in 1:length(x)
  str = @sprintf("%f",x[i])
end
dt_sprintf_stdout

# <codecell>

x = rand(1024)
dt_printf_stdout = @elapsed for i in 1:length(x)
  str = @printf("%f ",x[i])
end
dt_printf_stdout

# <codecell>

dt_print_loop_stdout = @elapsed for i in 1:length(x)
  println(x[i])
  end
dt_print_loop_stdout

# <codecell>

dt_print_stdout = @elapsed println(x)
dt_print_stdout

# <codecell>

dt_printf_loop_file = @elapsed open("outfile", "w") do f
    { @printf(f, "%f ", x[i]) for i in 1:length(x) }
    end
dt_printf_loop_file

# <codecell>

N = 1024*1024
x = rand(N);

# <codecell>

dt_printf_loop_file = @elapsed open("outfile", "w") do f
    { @printf(f, "%f ", x[i]) for i in 1:length(x) }
    end
dt_printf_loop_file

# <codecell>

dt_print_loop_file = @elapsed open("outfile", "w") do f
{ println(f,x[i]) for i in 1:length(x) }
    end
dt_print_loop_file

# <codecell>

function dt_writedlm(N::Int) 
  x = rand(N); 
  try rm("outfile_$N.csv") catch end
  @elapsed writedlm("outfile_$N.csv",x)
end
function dt_readdlm(N::Int) 
  dt_writedlm(N)
  @elapsed x = readdlm("outfile_$N.csv")
end
println(1024,"  ", dt_writedlm(1024))
println(1024*1024,"  ",dt_readdlm(1024*1024))

# <codecell>

function dt_print_file(N::Int) 
  x = rand(N); 
  try rm("outfile_$N.txt") catch end
  @elapsed open("outfile_$N.txt", "w") do f  print(f,x)  end 
end
println(1024,"  ", dt_print_file(1024))
println(1024*1024,"  ",dt_print_file(1024*1024))

# <codecell>

function dt_write_loop(N::Int)
  x = rand(N); 
  try rm("outfile_$N") catch end
  @elapsed open("outfile_$N", "w") do f
    { write(f,x[i])  for i in 1:length(x) }
end
end
println(1024,"  ", dt_write_loop(1024))
println(1024*1024,"  ",dt_write_loop(1024*1024))

# <codecell>

function dt_serialize(N::Int)
  x = rand(N); 
  try rm("outfile_$N.bin") catch end
  @elapsed open("outfile_$N.bin", "w") do f
      serialize(f,x)
  end
end
function dt_deserialize(N::Int)
  dt_serialize(N)
  @elapsed open("outfile_$N.bin", "r") do f
    x = deserialize(f)
  end
end
println(1024,"  ", dt_serialize(1024))
println(1024,"  ", dt_deserialize(1024))
println(1024*1024,"  ",dt_serialize(1024*1024))
println(1024*1024,"  ",dt_deserialize(1024*1024))

# <codecell>

function dt_write(N::Int)
  x = rand(N);
  try rm("outfile_$N.bin2") catch end
  @elapsed open("outfile_$N.bin2", "w") do f
    write(f,x) 
  end
end
function dt_read(N::Int)
  x = rand(N);
  @elapsed open("outfile_$N.bin2", "r") do f
  { x[i] = read(f,typeof(x[i])) for i in 1:N }
  end
end
println(1024,"  ", dt_write(1024))
println(1024,"  ",dt_read(1024))
println(1024*1024,"  ", dt_write(1024*1024))
println(1024*1024,"  ",dt_read(1024*1024))

# <codecell>

using HDF5, JLD

# <codecell>

function dt_hdf_save(N::Int)
  x = rand(N);
  try rm("output_$N.hdf5") catch end
  @elapsed @save "output_$N.hdf5" x
end
function dt_hdf_load(N::Int)
  dt_hdf_save(N)
  @elapsed @load "output_$N.hdf5" x
end

println(1024,"  ", dt_hdf_save(1024))
println(1024,"  ", dt_hdf_load(1024))
println(1024*1024,"  ",dt_hdf_save(1024*1024))
println(1024*1024,"  ",dt_hdf_load(1024*1024))

# <codecell>

function dt_hdf_h5write(N::Int)
  x = rand(N);
  try rm("output2_$N.hdf5") catch end
  @elapsed h5write("output2_$N.hdf5", "mydata/x", x)
end
println(1024,"  ", dt_hdf_h5write(1024))
println(1024*1024,"  ",dt_hdf_h5write(1024*1024))

# <codecell>

Nlog2 = 25
inlist= [ 2^i for i=1:Nlog2 ]
loginlist= [ log10(2.0^i) for i=1:Nlog2 ]
outlist1 = log10(convert(Array{Float64},map(dt_writedlm,inlist)))
outlist2 = log10(convert(Array{Float64},map(dt_serialize,inlist)))
outlist3 = log10(convert(Array{Float64},map(dt_hdf_save,inlist)))
outlist6 = log10(convert(Array{Float64},map(dt_readdlm,inlist)))
outlist5 = log10(convert(Array{Float64},map(dt_deserialize,inlist)))
outlist4 = log10(convert(Array{Float64},map(dt_hdf_load,inlist)))

# <codecell>

using PyPlot
plot(loginlist,outlist1,color="red",label="writedlm");
plot(loginlist,outlist2,color="blue", label="serialize");
plot(loginlist,outlist3,color="green", label="hdf_save");
plot(loginlist,outlist6,color="magenta",label="readdlm");
plot(loginlist,outlist5,color="cyan", label="deserialize");
plot(loginlist,outlist4,color="black", label="hdf_load");
xlabel("log Data Size"); ylabel("log (Time/sec)");  
title("I/O Performance"); 
legend(loc=2)

# <codecell>


