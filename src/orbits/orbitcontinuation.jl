function OrbitContinuation(f!, hp_ini::HopfPoint{U}, params, pmin::U, pmax::U, 
                           Δs::U, maxsteps::T, tol::U, ite::T,
                           integtol::U, integorder::T; 
                           θ::U = 0.0, integmaxsteps::T = 500, parse_eqs::Bool = true) where {U<:Real, T<:Integer}

    ###
    n = length(hp_ini.x) + 1
    ###
    x = Array{U,2}(undef, maxsteps, n)

    ###
    xT = TaylorN.(1:n, order = 1)

    ###
    x0 = [copy(hp_ini.x) ; U(2.0*π / abs(imag(hp_ini.value)))]
    Φ0 = zeros(U, n)
    Φ0[1:n-2] = sin(θ) * real(hp_ini.vector) + cos(θ) * imag(hp_ini.vector) 
    dx0 = [ abs(imag(hp_ini.value)) * (
        cos(θ) * real(hp_ini.vector) -
        sin(θ) * imag(hp_ini.vector)
    ) ; 0.0; 0.0]

    normalize!(Φ0)

    ###
    x1 = copy(x0)
    Φ1 = copy(Φ0)

    ###
    cache = TaylorIntegration.init_cache(0.0, xT, integorder, f!, params; parse_eqs)

    ###
    J = zeros(U, n, n)
    F = zeros(U, n)
    v = zeros(U, n)
    v[n] = one(U)

    ###
    x[1, :] .= x0

    ###
    i = 2

    while i <= maxsteps && (pmin <= x1[n-1] <= pmax)

        for k in 1:n
            x1[k] = x0[k] + Δs * Φ0[k]
            xT[k][0][1] = x1[k]
            for j in 1:n
                if k == j
                    xT[k][1][j] = 1.0
                else
                    xT[k][1][j] = 0.0
                end
            end
        end

        taylorinteg_optim!(f!, xT, 0.0, 1.0, integtol, cache, params; maxsteps = integmaxsteps)

        OrbitJacobian!(J, xT, dx0, Φ0, n)
        OrbitSystem!(F, xT, x1, x0, dx0, Φ0, Δs, n)

        j = 1

        while j <= ite && norm(F) > tol

            if abs(det(J)) == 0.0
                j = ite + 1
                break
            end

            x1 .-= J \ F 

            for k in 1:n
                xT[k][0][1] = x1[k]
                for l in 1:n
                    if k == l
                        xT[k][1][l] = 1.0
                    else
                        xT[k][1][l] = 0.0
                    end
                end
            end
    
            taylorinteg_optim!(f!, xT, 0.0, 1.0, integtol, cache, params; maxsteps = integmaxsteps)

            OrbitJacobian!(J, xT, dx0, Φ0, n)
            OrbitSystem!(F, xT, x1, x0, dx0, Φ0, Δs, n)

            j += 1

        end

        if j > ite
            @warn("System norm did not converge to tolerance. Stopping continuation.")
            break
        end

        x[i, :] .= x1

        Φ1 .= J \ v

        normalize!(Φ1)

        x0 .= x1
        Φ0 .= Φ1

        f!(dx0, x0, params, 0.0)

        i += 1

        print(" \r i = $i : ||F|| = $(norm(F))")

    end

    return x[1:i-1, :]

end