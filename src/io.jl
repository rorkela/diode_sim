function plotstate(params::Parameters,s::State)
    dx=params.dx
    dt=params.dt
    charge=[charge_den(params,s,i) for i=1:params.n]
    ef=(s.v[1:end-1]-s.v[2:end])/dx;
    jn=current_den_n(params,s,s.v./(k_B*T/q))
    jp=current_den_p(params,s,s.v./(k_B*T/q))
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
