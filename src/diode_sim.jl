module diode_sim
using Plots
using LinearAlgebra
using UnicodePlots
include("definations.jl")
include("io.jl")
include("poisson.jl")
include("continuity.jl")
println("Load zyada hora")

function run_sim(params::Parameters,vbound1::Float64,vboundn::Float64)
    # Initialize state
    s=State(v_init(params,vbound1,vboundn), n_init(params), p_init(params));
    for i = 1:params.t_steps
       sprev=deepcopy(s)
       for i=1:20
           poisson!(params,s, vbound1, vboundn)
           continuity!(params,s,sprev)
       end
    # if i%100==0 println("Time step: ", i) end
    end
    return s
end

function iv_curve(params::Parameters)
    v=[params.v_start + (i-1)*(params.v_end - params.v_start)/(params.v_steps - 1) for i = 1:params.v_steps]
    currents=Vector{Float64}(undef,params.v_steps)
    Threads.@threads for i in 1:params.v_steps
        vapp=v[i]
        s=run_sim(params,vapp,0.0)
        jn=current_den_n(params,s,s.v./(k_B*T/q))
        jp=current_den_p(params,s,s.v./(k_B*T/q))
        jsum=jn .+ jp
        currents[i]=jsum[end]
        println("Applied Voltage: ", vapp, " V, Current Density: ", jsum[end], " A/m^2")
    end
    plot(v, currents, title="IV Curve", xlabel="Applied Voltage (V)", ylabel="Current Density (A/m^2)", legend=false)
    gui()
    readline()
    return v,currents
end
end # module
