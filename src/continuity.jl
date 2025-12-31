norm(v) = v/(k_B*T/q)
B(z) = (z==0.0) ? (1.0) : (z/(exp(z)-1.0))

function jacobi_n(params::Parameters,s::State,sprev::State,vnorm::Vector{Float64},vprevnorm::Vector{Float64})
    mu=mu_n
    dx=params.dx
    dt=params.dt
    con=k_B*T*mu/(2.0*q*dx*dx)
    jl=Vector{Float64}(undef, params.n)
    jm=Vector{Float64}(undef, params.n)
    ju=Vector{Float64}(undef, params.n)
    for i = 1:params.n
        if i==1 || i==params.n
            jl[i]=0.0
            jm[i]=1.0
            ju[i]=0.0
        else
            jl[i]= -con*B(vnorm[i-1]-vnorm[i])
            jm[i]= params.C_Rr*s.p[i]/2+con*B(-vnorm[i-1]+vnorm[i])+con*B(vnorm[i]-vnorm[i+1])+1/dt
            ju[i]= -con*B(-vnorm[i]+vnorm[i+1])
        end
    end
    
    return Tridiagonal(jl[2:end], jm, ju[1:(end-1)])
end

function jacobi_p(params::Parameters,s::State,sprev::State,vnorm::Vector{Float64},vprevnorm::Vector{Float64})
    mu=mu_p
    dx=params.dx
    dt=params.dt
    con=k_B*T*mu/(2.0*q*dx*dx)
    jl=Vector{Float64}(undef, params.n)
    jm=Vector{Float64}(undef, params.n)
    ju=Vector{Float64}(undef, params.n)
    for i = 1:params.n
        if i==1 || i==params.n
            jl[i]=0.0
            jm[i]=1.0
            ju[i]=0.0
        else
            jl[i]= -con*B(-vnorm[i-1]+vnorm[i])
            jm[i]= params.C_Rr*s.n[i]/2+con*B(vnorm[i-1]-vnorm[i])+con*B(-vnorm[i]+vnorm[i+1])+1/dt
            ju[i]= -con*B(vnorm[i]-vnorm[i+1])
        end
    end

    return Tridiagonal(jl[2:end], jm, ju[1:(end-1)])
end

function residual_n(params::Parameters,s::State,sprev::State,vnorm::Vector{Float64},vprevnorm::Vector{Float64})
    R=Vector{Float64}(undef, params.n)
    dx=params.dx
    dt=params.dt
    cur=current_den_n(params,s,vnorm)
    curp=current_den_n(params,sprev,vprevnorm)
    for i = 1:params.n
        if i==1
            R[i]=s.n[i]-n_i^2/params.N_A
        elseif i==params.n
            R[i]=s.n[i]-params.N_D
        else
            R[i]=(s.n[i]-sprev.n[i])/dt - ((cur[i] - cur[i-1] + curp[i] - curp[i-1])/(2*q*dx)) - (params.Gr - params.C_Rr * (s.n[i]*s.p[i]/2 + sprev.n[i]*sprev.p[i]/2 - n_i*n_i));
        end
    end
    return R
end

function residual_p(params::Parameters,s::State,sprev::State,vnorm::Vector{Float64},vprevnorm::Vector{Float64})
    R=Vector{Float64}(undef, params.n)
    dx=params.dx
    dt=params.dt
    cur=current_den_p(params,s,vnorm)
    curp=current_den_p(params,sprev,vprevnorm)
    for i = 1:params.n
        if i==1
            R[i]=s.p[i]-params.N_A
        elseif i==params.n
            R[i]=s.p[i]-n_i^2/params.N_D
        else
            R[i]=(s.p[i]-sprev.p[i])/dt - (-(cur[i] - cur[i-1] + curp[i] - curp[i-1])/(2*q*dx)) - (params.Gr - params.C_Rr * (s.n[i]*s.p[i]/2 + sprev.n[i]*sprev.p[i]/2 - n_i*n_i));
        end
    end
    return R
    
end

function continuity!(params::Parameters,s::State,sprev::State)
    vnorm=s.v ./(k_B*T/q)
    dx=params.dx
    dt=params.dt
    vprevnorm=sprev.v ./(k_B*T/q)
    maxiter=10
    for i in 1:maxiter
    update =jacobi_p(params,s,sprev,vnorm,vprevnorm)\(-1.0 .* residual_p(params,s,sprev,vnorm,vprevnorm))
    update[1]=0.0
    update[end]=0.0
    s.p = s.p .+ 0.2 .* update
    end
    for i in 1:maxiter
    update =jacobi_n(params,s,sprev,vnorm,vprevnorm)\(-1.0 .* residual_n(params,s,sprev,vnorm,vprevnorm))
    update[1]=0.0
    update[end]=0.0
    s.n = s.n .+ 0.2 .* update
    end
end
