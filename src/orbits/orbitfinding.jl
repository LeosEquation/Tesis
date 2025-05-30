
function HOrbitFinding!(f!, H, E, bc!, bcU!, x_ini::Array{U, 1}, params,
                       tol::U, ite::T, Npoints::T,
                       integtol::U, integorder::T;
                       integmaxsteps::T = 500, parse_eqs::Bool = true) where {U<:Real, T<:Integer}

    @info "This function is running with $(Threads.nthreads()) threads."

    ###
    n = length(x_ini)
    m = Npoints
    N = m*(n-1) + 1

    Δτ = one(U) / (m - 1)
    τ = LinRange(zero(U), one(U), Npoints)

    ###
    xT = TaylorN.([j for i in 1:m-1, j in 1:n], order = 1)
    Mvec = [zeros(U, n-1, n-1) for i in 1:m-1]
    Haux = TaylorN(1, order = 1)
    xaux = TaylorN.(1:n, order = 1)

    ###

    x = zeros(U, N)

    comp3 = 0
    for comp1 in 1:m
        for comp2 in 1:n-1
            comp3 += 1
            x[comp3] = x_ini[comp2]
        end
    end
    
    x[N] = x_ini[n]

    for comp1 in 1:n
        xaux[comp1][0][1] = x_ini[comp1]
    end

    Haux .= H(xaux, params, 0.0)

    ###

    x_temp = copy(x_ini)

    cache_ini = TaylorIntegration.init_cache(zero(U), x_ini, integorder, f!, params; parse_eqs)

    comp3 = n - 1

    for comp1 in 2:m-1

        taylorinteg_optim!(f!, bcU!, x_temp, 0.0, Δτ, integtol, cache_ini, params; maxsteps = integmaxsteps)

        for comp2 in 1:n-1

            comp3 += 1

            x[comp3] = x_temp[comp2]

        end

    end

    ###
    J = zeros(U, N, N)
    F = zeros(U, N)

    comp3 = 0
    for comp1 in 1:m-1
        for comp2 in 1:n-1
            comp3 += 1
            xT[comp1, comp2][0][1] = x[comp3]
            for comp4 in 1:n
                if comp2 == comp4
                    xT[comp1, comp2][1][comp4] = 1.0
                else
                    xT[comp1, comp2][1][comp4] = 0.0
                end
            end
        end
        xT[comp1, n][0][1] = x[N]
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

    ###

    caches = [
                TaylorIntegration.init_cache(zero(U), xT[1, :], integorder, f!, params; parse_eqs) 
                for _ in 1:Threads.nthreads()
             ]

    ###

    Threads.@threads for comp1 in 1:m-1
        tid = Threads.threadid()
        xTview = @view xT[comp1, :]
        taylorinteg_optim!(f!, bc!, xTview, 0.0, Δτ, integtol, caches[tid], params; maxsteps = integmaxsteps)
    end

    HOrbitJacobian!(J, xT, Haux, n, m, N)
    HOrbitSystem!(F, xT, Haux[0][1], E, x, n, m, N)

    # @show F

    j = 1

    while j <= ite && norm(F) > tol

        println(" Newton iteration = $j \t \t ||F|| = $(norm(F)) \t \t  ")

        if abs(det(J)) == 0.0
            break
        end

        x .-= J \ F

        for comp1 in 1:n-1
            xaux[comp1][0][1] = x[comp1]
        end
        xaux[n][0][1] = x[N]

        Haux .= H(xaux, params, 0.0)
        
        comp3 = 0
        for comp1 in 1:m-1
            for comp2 in 1:n-1
                comp3 += 1
                xT[comp1, comp2][0][1] = x[comp3]
                for comp4 in 1:n
                    if comp2 == comp4
                        xT[comp1, comp2][1][comp4] = 1.0
                    else
                        xT[comp1, comp2][1][comp4] = 0.0
                    end
                end
            end
            xT[comp1, n][0][1] = x[N]
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


        Threads.@threads for comp1 in 1:m-1
            tid = Threads.threadid()
            xTview = @view xT[comp1, :]
            taylorinteg_optim!(f!, bc!, xTview, 0.0, Δτ, integtol, caches[tid], params; maxsteps = integmaxsteps)
        end        

        HOrbitJacobian!(J, xT, Haux, n, m, N)
        HOrbitSystem!(F, xT, Haux[0][1], E, x, n, m, N)

        j += 1

    end

    if norm(F) > tol
        @warn("System norm did not converge to tolerance. Stopping continuation.")
    end

    for comp1 in 1:n-1
        x_ini[comp1] = x[comp1]
    end
    x_ini[n] = x[N]

    for comp1 in 1:m-1
        for comp2 in 1:n-1
            for comp3 in 1:n-1
                Mvec[comp1][comp2, comp3] = xT[comp1, comp2][1][comp3]
            end
        end
    end
    
    pS = pschur!(Mvec, :L)

end