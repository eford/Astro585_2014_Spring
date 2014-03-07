# <codecell>

function pmap_darray_1d!(arr_dist::DArray, func::Function, arr_output::Array)
@assert(size(arr_dist) == size(arr_output) )
proclist = procs(arr_dist)
refs = [@spawnat proclist[p] map(func, localpart(arr_dist)) for p in 1:length(proclist) ]
for p in 1:length(proclist)
  arr_sub = fetch(refs[p])
  subidx = arr_dist.indexes[p]
  for i in 1:length(subidx[1])
    arr_output[subidx[1][i]] = arr_sub[i]
  end
end
arr_output
end

# <codecell>

function pmap_darray_2d!(arr_dist::DArray, func::Function, arr_output::Array)
@assert(size(arr_dist) == size(arr_output) )
proclist = procs(arr_dist)
refs = [@spawnat proclist[p] map(func, localpart(arr_dist)) for p in 1:length(proclist) ]
for p in 1:length(proclist)
  arr_sub = fetch(refs[p])
  subidx = arr_dist.indexes[p]
  for i in 1:length(subidx[1]), j in 1:length(subidx[2])
    arr_output[subidx[1][i],subidx[2][j]] = arr_sub[i,j]
  end
end
arr_output
end

# <codecell>

function pmap_darray_3d!(arr_dist::DArray, func::Function, arr_output::Array)
@assert(size(arr_dist) == size(arr_output) )
proclist = procs(arr_dist)
refs = [@spawnat proclist[p] map(func, localpart(arr_dist)) for p in 1:length(proclist) ]
for p in 1:length(proclist)
  arr_sub = fetch(refs[p])
  subidx = arr_dist.indexes[p]
  for i in 1:length(subidx[1]), j in 1:length(subidx[2]), k in 1:length(subidx[3])
    arr_output[subidx[1][i],subidx[2][j],subidx[3][k]] = arr_sub[i,j,k]
  end
end
arr_output
end

# <codecell>

function pmap_darray!(arr_dist::DArray, func::Function, arr_output::Array)
  numd = length(size(arr_dist))
  @assert(1 <= numd <= 3)
  if numd == 1    return pmap_darray_1d!(arr_dist,func,arr_output)  end
  if numd == 2    return pmap_darray_2d!(arr_dist,func,arr_output)  end
  if numd == 3    return pmap_darray_3d!(arr_dist,func,arr_output)  end
end

# <codecell>

function pmap_darray(arr_dist::DArray, func::Function)
  arr_output = cell(size(arr_dist))
  pmap_darray!(arr_dist,func,arr_output)
  arr_output
end

