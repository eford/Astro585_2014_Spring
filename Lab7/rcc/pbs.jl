function check_pbs_tasknum_one()
  pbs_tasknum = int(ENV["PBS_TASKNUM"])
  if (pbs_tasknum!=1) error("Somehow Julia started with a tasknum !=1") end
end

function read_pbs_nodefile(nodefilename::String = ENV["PBS_NODEFILE"]) 
  nodefile = open(nodefilename)
  local nodelist 
  try
    nodelist = readlines(nodefile)
  finally
    close(nodefile)
  end
  @assert(length(nodelist)>=1)  # ensure read at least one nodename 
  # trim CR/LF at end of each node name
  { nodelist[i] = chomp(nodelist[i]) for i in 1:length(nodelist) }
  return nodelist
end

function addprocs_pbs( ; nodefilename::String = ENV["PBS_NODEFILE"], all_procs_work::Integer = -1 )
  nodelist = read_pbs_nodefile(nodefilename)
  # Automatically decide whether to set the first node as a worker if not specified
  if all_procs_work == -1
     if length(nodelist) <= 16
        all_procs_work = 1
     else
        all_procs_work = 0
     end
  end
  first_proc = all_procs_work > 0 ? 1 : 2
  proclist = try
    addprocs(nodelist[first_proc:end])
  catch
    warn("addprocs threw an exception")
    workers()
  end
end

@everywhere function print_node_status()
  println("# myid= ",myid())
  println("# nprocs= ",nprocs())
  println("# nworkers= ",nworkers())
  flush_cstdio()
end

@everywhere function srand_all_procs(seed::Integer; same::Bool=false, debug::Bool = false)
  srand(seed)
  worker_seeds = same ? fill(seed, length(workers())) : rand(1:typemax(Uint32),length(workers()) )
  @sync for (proc,seed) in zip(workers(),worker_seeds)
     if debug println("Running srand($seed) on proc ",proc) end
     @spawnat(proc, srand(seed) )
  end
  if debug 
     @sync for proc in workers()
        println("rand() from proc returns ",@fetchfrom(proc, rand() ) )
     end
  end
end

