
const q = 1.602176634e-19
const k_B = 1.380649e-23
const T = 300.0
const ep_0 = 8.854187817e-12
const ep_r_Si = 11.7
const n = 1000
const dx = 1e-8
const dt = 1e-12
const mu_n = 0.135
const mu_p = 0.045
const N_A_ptype = 1e18
const N_D_ntype = 1e18
const n_i=1.1e16

mutable struct State
    v::Vector{Float64}
    n::Vector{Float64}
    p::Vector{Float64}
end

function Na(index)
    if index < n/2
        return N_A_ptype
    else
        return 0.0
    end
end

function Nd(index)
    if index >= n/2
        return N_D_ntype
    else
        return 0.0
    end
end

function v_init()
    return [vbound1*(n-i)/(n-1)+vboundn*(i-1)/(n-1) for i = 1:n]
end
function n_init()
    [i<n/2 ? n_i^2/N_A_ptype : N_D_ntype for i = 1:n]
end

function p_init()
    [i<n/2 ? N_A_ptype : n_i^2/N_D_ntype for i = 1:n]
end

function permittivity(index)
    return ep_r_Si*ep_0
end

function plotstate(s::State)
    gr()
    p1=plot(
        s.v,
        title = "Electrostatic Potential",
        xlabel = "Position (nm)",
        ylabel = "Potential (V)",
    )
    p2=plot(
        s.n,
        title = "Electron Concentration",
        xlabel = "Position (nm)",
        ylabel = "Concentration (m^-3)",
        yscale = :log10,
    )
    p3=plot(
        s.p,
        title = "Hole Concentration",
        xlabel = "Position (nm)",
        ylabel = "Concentration (m^-3)",
        yscale = :log10,
    )
    plot(p1, p2, p3, layout = (3, 1))
    gui()
    readline()
end
