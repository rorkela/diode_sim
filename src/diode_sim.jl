module diode_sim
using Plots
using LinearAlgebra
using UnicodePlots
include("definations.jl")
include("poisson.jl")
include("continuity.jl")
println("Load zyada hora")

vbound1=-5.0
vboundn=0.0
# Initialize state
s=State(v_init(), n_init(), p_init());
plotstate(s)
for i = 1:1000
    sprev=deepcopy(s)
    for i=1:20
        poisson!(s, vbound1, vboundn)
        continuity!(s,sprev)
    end
    plotstate(s)
end
end # module
