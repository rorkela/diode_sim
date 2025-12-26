norm(v) = v/(k_B*T/q)
B(z) = (z==0.0) ? (1.0) : (z/(exp(z)-1.0))
function current_den_n(s::State,v::Vector{Float64})
    c=k_B*T*mu_n/q
    Jn=Vector{Float64}(undef, n)
    for i = 1:n
        if i==n
            Jn[i]=c*(N_D_ntype-s.n[n])
        else
            Jn[i]=c*(s.n[i+1]*B(v[i+1]-v[i])-s.n[i]*B(v[i]-v[i+1]))
        end
    end
    return Jn
end

function current_den_p(s::State,v::Vector{Float64})
    c=k_B*T*mu_p/q
    Jp=Vector{Float64}(undef, n)
    for i = 1:n
        if i==n
            Jp[i]=-c*(n_i^2/N_D_ntype-s.p[n])
        else
            Jp[i]=-(k_B*T*mu_p/q)*(s.p[i+1]*B(-v[i+1]+v[i])-s.p[i]*B(-v[i]+v[i+1]))
        end
    end
    return Jp
end

function jacobi_n(s,sprev,vnorm,vprevnorm)
    mu=mu_n
    con=k_B*T*mu/(2.0*q*dx*dx)
    jl=Vector{Float64}(undef, n)
    jm=Vector{Float64}(undef, n)
    ju=Vector{Float64}(undef, n)
    for i = 1:n
        if i==1 || i==n
            jl[i]=0.0
            jm[i]=1.0
            ju[i]=0.0
        else
            jl[i]= -con*B(vnorm[i-1]-vnorm[i])
            jm[i]= C_Rr*s.p[i]/2+con*B(-vnorm[i-1]+vnorm[i])+con*B(vnorm[i]-vnorm[i+1])+1/dt
            ju[i]= -con*B(-vnorm[i]+vnorm[i+1])
        end
    end

    return Tridiagonal(jl[2:end], jm, ju[1:(end-1)])
end

function jacobi_p(s,sprev,vnorm,vprevnorm)
    mu=mu_p
    con=k_B*T*mu/(2.0*q*dx*dx)
    jl=Vector{Float64}(undef, n)
    jm=Vector{Float64}(undef, n)
    ju=Vector{Float64}(undef, n)
    for i = 1:n
        if i==1 || i==n
            jl[i]=0.0
            jm[i]=1.0
            ju[i]=0.0
        else
            jl[i]= -con*B(-vnorm[i-1]+vnorm[i])
            jm[i]= C_Rr*s.n[i]/2+con*B(vnorm[i-1]-vnorm[i])+con*B(-vnorm[i]+vnorm[i+1])+1/dt
            ju[i]= -con*B(vnorm[i]-vnorm[i+1])
        end
    end

    return Tridiagonal(jl[2:end], jm, ju[1:(end-1)])
end
function residual_n(s::State,sprev::State,vnorm,vprevnorm)
    R=Vector{Float64}(undef, n)
    cur=current_den_n(s,vnorm)
    curp=current_den_p(s,vprevnorm)
    for i = 1:n
        if i==1
            R[i]=s.n[i]-n_i^2/N_A_ptype
        elseif i==n
            R[i]=s.n[i]-N_D_ntype
        else
            R[i]=(s.n[i]-sprev.n[i])/dt - ((cur[i] - cur[i-1] + curp[i] - curp[i-1])/(2*q*dx)) - (Gr - C_Rr * (s.n[i]*s.p[i]/2 + sprev.n[i]*sprev.p[i]/2 - n_i*n_i));
        end
    end
    return R
end

function residual_p(s::State,sprev::State,vnorm,vprevnorm)
    R=Vector{Float64}(undef, n)
    cur=current_den_n(s,vprevnorm)
    curp=current_den_p(s,vnorm)
    for i = 1:n
        if i==1
            R[i]=s.p[i]-N_A_ptype
        elseif i==n
            R[i]=s.p[i]-n_i^2/N_D_ntype
        else
            R[i]=(s.p[i]-sprev.p[i])/dt - (-(curp[i] - curp[i-1] + cur[i] - cur[i-1])/(2*q*dx)) - (Gr - C_Rr * (s.n[i]*s.p[i]/2 + sprev.n[i]*sprev.p[i]/2 - n_i*n_i));
        end
    end
    return R
    
end

function continuity!(s::State,sprev::State)
    vnorm=s.v ./(k_B*T/q)
    vprevnorm=sprev.v ./(k_B*T/q)
    maxiter=10
    for i in 1:maxiter
    update =jacobi_n(s,sprev,vnorm,vprevnorm)\(-1.0 .* residual_n(s,sprev,vnorm,vprevnorm))
    s.n = s.n .+ 0.2 .* update
    end
    for i in 1:maxiter
    update =jacobi_p(s,sprev,vnorm,vprevnorm)\(-1.0 .* residual_p(s,sprev,vnorm,vprevnorm))
    s.p = s.p .+ 0.2 .* update
    end
end
