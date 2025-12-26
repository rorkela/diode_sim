norm(v) = v/(k_B*T/q)
B(z) = (z==0.0)?(1.0):(z/(exp(z)-1.0))
function current_den_n(s::State,v::Vector{Float64})
    c=k_B*T*mu_n/q
    Jn=Vector{Float64}(undef, n)
    for i = 1:n
        if i==1
            Jn[i]=c*(s.n[i+1] - n_teq)
        elseif i==n
            Jn[i]=c*(n_teq-s.n[n])
        else
            Jn[i]=c*(s.n[i+1]*B(v[i+1]-v[i])-s.n[i]*B(v[i]-v[i+1]))
        end
    end
    return Jn
end

function current_den_p(s::State,v::Vector{Float64})
    c=k_B*T*mu_n/q
    Jp=Vector{Float64}(undef, n)
    for i = 1:n
        if i==1
            Jn[i]=-c*(s.p[i+1] - p_teq)
        elseif i==n
            Jn[i]=-c*(p_teq-s.p[n])
        else
            Jp[i]=-(k_B*T*mu_p/q)*(s.p[i+1]*B(-v[i+1]+v[i])-s.p[i]*B(-v[i]+v[i+1]))
        end
    end
    return Jp
end

function jacobi_n(s,sprev,vnorm,vprevnorm)
    mu=mu_n
    

end

end
function continuity!(s::State,sprev::State)
    vnorm=[s.v ./(k_B*T/q)]
    vprevnorm=[sprev.v ./(k_B*T/q)]
    maxiter=20
    for i in 1:maxiter
    update =jacobi_n(s,sprev,vnorm,vprevnorm)\(-1.0 .* residual_n())
    s.n = s.n .+ 0.2 .* update
    end
    for i in 1:maxiter
    update =jacobi_p(s,sprev,vnorm,vprevnorm)\(-1.0 .* residual_p())
    s.p = s.p .+ 0.2 .* update
    end
end
