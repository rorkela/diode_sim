
const q = 1.602176634e-19
const k_B = 1.380649e-23
const T = 300.0
const ep_0 = 8.854187817e-12
const ep_r_Si = 11.7
const n = 128
#const dx = 1e-6
#const dt = 1e-20
const mu_n = 0.135
const mu_p = 0.045
const N_A_ptype = 2.1e18
const N_D_ntype = 2.1e18
const n_i=1.1e16
const C_Rr=1.1e-8
const Gr=1e16
mutable struct State
    v::Vector{Float64}
    n::Vector{Float64}
    p::Vector{Float64}
end

mutable struct Parameters
    dx::Float64
    dt::Float64
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

function v_init(vbound1,vboundn)
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

function charge_den(s::State,i)
    return q*(Nd(i)-Na(i)-s.n[i]+s.p[i])
end

function current_den_n(s::State,v::Vector{Float64},dx::Float64)
    c=k_B*T*mu_n/q
    Jn=Vector{Float64}(undef, n-1)
    for i = 1:n-1
            Jn[i]=c*(s.n[i+1]*B(v[i+1]-v[i])-s.n[i]*B(v[i]-v[i+1]))
    end
    return Jn
end

function current_den_p(s::State,v::Vector{Float64},dx::Float64)
    c=k_B*T*mu_p/q
    Jp=Vector{Float64}(undef, n-1)
    for i = 1:n-1
            Jp[i]=-(k_B*T*mu_p/q)*(s.p[i+1]*B(-v[i+1]+v[i])-s.p[i]*B(-v[i]+v[i+1]))
    end
    return Jp
end

function plotstate(s::State,params::Parameters)
    dx=params.dx
    dt=params.dt
    charge=[charge_den(s,i) for i=1:n]
    ef=(s.v[1:end-1]-s.v[2:end])/dx;
    jn=current_den_n(s,s.v./(k_B*T/q),dx)
    jp=current_den_p(s,s.v./(k_B*T/q),dx)
    jsum=jn .+ jp
    gr()
    default(
        titlefontsize = 9,
        guidefontsize = 8,
        tickfontsize  = 7,
        legendfontsize = 7,
        linewidth = 2,
    )
    p1=plot(s.v,title = "Electrostatic Potential",ylabel = "V",legend = false)
    p2=plot(ef,title="Electric Field",ylabel="E",legend=false)
    p3=plot(s.n,yscale=:log10,label="n",title="Carrier Conc.")
    plot!(p3,s.p,yscale=:log10,label="p")
    p4=plot(jn,linestyle=:dash,title="Current Density",ylabel="n",legend=false)
    plot!(p4,linestyle=:dash,jp,label="p")
    plot!(p4,linestyle=:solid,jsum,label="Total")
    p5=plot(charge,title="Charge Density",legend=false)
    plot(p1, p2, p3,p4,p5, layout = (3, 2))
    gui()
    readline()
end
