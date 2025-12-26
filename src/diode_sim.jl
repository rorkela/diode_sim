module diode_sim
using Plots
using LinearAlgebra
using UnicodePlots
include("definations.jl")
include("poisson.jl")
include("continuity.jl")
println("Load zyada hora")

vbound1=5.0
vboundn=0.0
# Initialize state
s=State(v_init(), n_init(), p_init());
plotstate(s)
poisson!(s, vbound1, vboundn)
continuity!(s)
plotstate(s)
end # module
