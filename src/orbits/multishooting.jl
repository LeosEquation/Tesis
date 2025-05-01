function OrbitContinuation(f!, hp_ini::HopfPoint{U}, params, pmin::U, pmax::U, 
                           Δs::U, maxsteps::T, tol::U, ite::T, Npoints::T,
                           integtol::U, integorder::T;
                           integmaxsteps::T = 500, parse_eqs::Bool = true) where {U<:Real, T<:Integer}

    ###
    n = length(hp_ini.x) + 1
    m = Npoints - 1
    N = m*(n-2) + 2

    Δτ = one(U) / m
    τ = LinRange(zero(U), one(U), Npoints)

    # return τ, Δτ

    ###
    x = Array{U, 2}(undef, maxsteps, N)

    ###
    xT = TaylorN.([j for i in 1:m, j in 1:n], order = 1)

    ###

    x0 = zeros(U, N)
    x00 = zeros(U, n)
    ω0 = abs(imag(hp_ini.value))

    comp3 = 0
    for comp1 in 1:m
        for comp2 in 1:n-2
            comp3 += 1
            x0[comp3] = hp_ini.x[comp2]
        end
    end
    
    x0[N-1] = hp_ini.x[n-1]
    x0[N] = 2π / ω0

    Φ0 = zeros(U, N)

    comp3 = 0
    for comp1 in 1:m 
        for comp2 in 1:n-2
            comp3 += 1
            Φ0[comp3] = cos(2π*τ[comp1]) * real(hp_ini.vector[comp2]) + sin(2π*τ[comp1]) * imag(hp_ini.vector[comp2])
        end
    end

    normalize!(Φ0)

    dx0 = [ ω0 * imag(hp_ini.vector) ; 0.0; 0.0]

    ###
    x1 = copy(x0)
    Φ1 = copy(Φ0)
    ω1 = copy(ω0)

    ###
    # caches = [TaylorIntegration.init_cache(zero(U), xT[i, :], integorder, f!, params; parse_eqs) for i in 1:m]

    cache = TaylorIntegration.init_cache(zero(U), xT[1, :], integorder, f!, params; parse_eqs)

    ###
    J = zeros(U, N, N)
    F = zeros(U, N)
    v = zeros(U, N)
    v[N] = one(U)

    ###
    x[1, :] .= x0

    ###
    i = 2

    while i <= maxsteps && (pmin <= x1[n-1] <= pmax)

        # println(" i = $i :")

        for comp1 in 1:N
            x1[comp1] = x0[comp1] + Δs * Φ0[comp1]
        end

        comp3 = 0
        for comp1 in 1:m 
            for comp2 in 1:n-2
                comp3 += 1
                xT[comp1, comp2][0][1] = x1[comp3]
                for comp4 in 1:n
                    if comp2 == comp4
                        xT[comp1, comp2][1][comp4] = 1.0
                    else
                        xT[comp1, comp2][1][comp4] = 0.0
                    end
                end
            end
            xT[comp1, n-1][0][1] = x1[N-1]
            xT[comp1, n][0][1] = x1[N]
            for comp4 in 1:n
                if n-1 == comp4
                    xT[comp1, n-1][1][comp4] = 1.0
                elseif n == comp4
                    xT[comp1, n][1][comp4] = 1.0
                else
                    xT[comp1, n-1][1][comp4] = 0.0
                    xT[comp1, n][1][comp4] = 0.0
                end
            end
        end

        # return xT

        for comp1 in 1:m
            xTview = @view xT[comp1, :]
            taylorinteg_optim!(f!, xTview, 0.0, Δτ, integtol, cache, params; maxsteps = integmaxsteps)
        end

        # return xT, x1

        OrbitJacobian!(J, xT, dx0, Φ0, n, m, N)
        OrbitSystem!(F, xT, x1, x0, dx0, Φ0, Δs, n, m, N)

        # return J, F

        j = 1

        # println("\t j = $j : norm = $(norm(F))")

        while j <= ite && norm(F) > tol

            if abs(det(J)) == 0.0
                j = ite + 1
                break
            end

            δ = J \ F

            # @show δ

            x1 .-= δ

            comp3 = 0
            for comp1 in 1:m 
                for comp2 in 1:n-2
                    comp3 += 1
                    xT[comp1, comp2][0][1] = x1[comp3]
                    for comp4 in 1:n
                        if comp2 == comp4
                            xT[comp1, comp2][1][comp4] = 1.0
                        else
                            xT[comp1, comp2][1][comp4] = 0.0
                        end
                    end
                end
                xT[comp1, n-1][0][1] = x1[N-1]
                xT[comp1, n][0][1] = x1[N]
                for comp4 in 1:n
                    if n-1 == comp4
                        xT[comp1, n-1][1][comp4] = 1.0
                    elseif n == comp4
                        xT[comp1, n][1][comp4] = 1.0
                    else
                        xT[comp1, n-1][1][comp4] = 0.0
                        xT[comp1, n][1][comp4] = 0.0
                    end
                end
            end

            for comp1 in 1:m
                xTview = @view xT[comp1, :]
                taylorinteg_optim!(f!, xTview, 0.0, Δτ, integtol, cache, params; maxsteps = integmaxsteps)
            end

            OrbitJacobian!(J, xT, dx0, Φ0, n, m, N)
            OrbitSystem!(F, xT, x1, x0, dx0, Φ0, Δs, n, m, N)

            j += 1

            # @show F

            # println("\t j = $j : norm = $(norm(F))")

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

        for comp1 in 1:n-2
            x00[comp1] = x0[comp1]
        end

        x00[n-1] = x0[N-1]
        x00[n] = x0[N]

        f!(dx0, x00, params, 0.0)

        i += 1

        print(" \r i = $i : ||F|| = $(norm(F))")

    end

    return x[1:i-1, :]

end