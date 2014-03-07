# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>
include("pmap_darray.jl")

addprocs(2)

# <codecell>

@everywhere function func(i, j; optparam= 0.0) 
  i*100+j+optparam
end

@everywhere function func(idx, input::DArray, output::DArray; optparam = 0.0) 
  #run(`sleep 3`)
  i,j = input[idx]
  output[idx] = func(i,j,optparam=optparam) 
end

# <codecell>

n = 4
m = 5
arr_local = [ (i,j) for i in 1:n, j in 1:m ]
arr_comp = [ func(i,j,optparam=0.5) for i in 1:n, j in 1:m ]
arr_dist = distribute(arr_local)
arr_comp

# <codecell>

pmap_darray(arr_dist,x->func(x[1],x[2],optparam=0.5))

# <codecell>

n = 8
arr_local = [ i for i in 1:n ]
arr_comp = [ x->x^2 for i in 1:n ]
arr_dist = distribute(arr_local)
pmap_darray(arr_dist,x->x^2)

# <codecell>


