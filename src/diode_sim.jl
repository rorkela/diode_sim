module diode_sim
using Plots
using LinearAlgebra
using UnicodePlots
include("definations.jl")
include("poisson.jl")
include("continuity.jl")
println("Load zyada hora")
function run_sim(dx,dt)
    params=Parameters(dx,dt)
    vbound1=0.0
    vboundn=0.0
    # Initialize state
    s=State(v_init(vbound1,vboundn), n_init(), p_init());
    for i = 1:100000
       sprev=deepcopy(s)
       for i=1:20
           poisson!(s, vbound1, vboundn,params)
           continuity!(s,sprev,params)
       end
    if i%100==0 println("Time step: ", i)
    end
    end
    plotstate(s,params)
end
end # module
