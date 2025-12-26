
function poisson!(s::State, bound1::Float64, boundn::Float64)
    Jl=Vector{Float64}(undef, n)
    Jm=Vector{Float64}(undef, n)
    Ju=Vector{Float64}(undef, n)
    R=Vector{Float64}(undef, n)
    for iter = 1:1
        for i = 1:n
            if i==1 || i==n
                Jl[i]=0.0
                Jm[i]=1.0
                Ju[i]=0.0
                R[i]=s.v[i] - ((i==1) ? bound1 : boundn)
            else
                Jl[i]=(permittivity(i)+permittivity(i-1))/(2.0*dx^2)
                Jm[i]=-(
                    permittivity(i+1)+permittivity(i-1)+2.0*permittivity(i)
                )/(2*dx^2)-q^2*(s.p[i]+s.n[i])/k_B/T
                Ju[i]=(permittivity(i)+permittivity(i+1))/(2.0*dx^2)
                R[i]=(
                    1/(2dx^2)
                )*(
                    (
                        permittivity(i)+permittivity(i-1)
                    )*s.v[i-1]+(
                        permittivity(i+1)+permittivity(i)
                    )*s.v[i+1]-(
                        2.0*permittivity(i)+permittivity(i+1)+permittivity(i-1)
                    )*s.v[i]
                )+q*(Nd(i)-Na(i)-s.n[i]+s.p[i])
            end
        end
        R=-1.0 .* R
        update=(Tridiagonal(Jl[2:end], Jm, Ju[1:(end-1)])\R);
        s.v[2:(end-1)] .+= 1.0 .* update[2:(end-1)]
        s.v[1] = bound1
        s.v[end] = boundn
    end


end
