
const q = 1.602176634e-19
const k_B = 1.380649e-23
const T = 300.0
const ep_0 = 8.854187817e-12
const ep_r_Si = 11.7
#const n = 128
#const dx = 1e-6
#const dt = 1e-20
const mu_n = 0.135
const mu_p = 0.045
#const N_A_ptype = 2.1e18
#const N_D_ntype = 2.1e18
const n_i=1.1e16
#const C_Rr=1.1e-8
#const Gr=1e16
mutable struct State
    v::Vector{Float64}
    n::Vector{Float64}
    p::Vector{Float64}
end
@kwdef mutable struct Parameters
    dx::Float64 = 1e-6
    dt::Float64 = 1e-20
    n::Int64 = 128
    N_A::Float64 = 2.1e18
    N_D::Float64 = 2.1e18
    C_Rr::Float64 = 1.1e-8
    Gr::Float64 = 1e16
    t_steps::Int64 = 1000
    v_start::Float64 = -2
    v_end::Float64 = 2
    v_steps::Int64 = 20

end
function Na(params::Parameters,index::Int64)
    if index < params.n/2
        return params.N_A
    else
        return 0.0
    end
end

function Nd(params::Parameters,index::Int64)
    if index >= params.n/2
        return params.N_D
    else
        return 0.0
    end
end

function v_init(params::Parameters,vbound1::Float64,vboundn::Float64)
    return [vbound1*(params.n-i)/(params.n-1)+vboundn*(i-1)/(params.n-1) for i = 1:params.n]
end
function n_init(params::Parameters)
    [i<params.n/2 ? n_i^2/params.N_A : params.N_D for i = 1:params.n]
end

function p_init(params::Parameters)
    [i<params.n/2 ? params.N_A : n_i^2/params.N_D for i = 1:params.n]
end

function permittivity(index)
    return ep_r_Si*ep_0
end

function charge_den(params,s::State,i)
    return q*(Nd(params,i)-Na(params,i)-s.n[i]+s.p[i])
end

function current_den_n(params::Parameters,s::State,v::Vector{Float64})
    c=k_B*T*mu_n/q
    Jn=Vector{Float64}(undef, params.n-1)
    for i = 1:params.n-1
            Jn[i]=c*(s.n[i+1]*B(v[i+1]-v[i])-s.n[i]*B(v[i]-v[i+1]))
    end
    return Jn
end

function current_den_p(params::Parameters,s::State,v::Vector{Float64})
    c=k_B*T*mu_p/q
    Jp=Vector{Float64}(undef, params.n-1)
    for i = 1:params.n-1
            Jp[i]=-(k_B*T*mu_p/q)*(s.p[i+1]*B(-v[i+1]+v[i])-s.p[i]*B(-v[i]+v[i+1]))
    end
    return Jp
end

