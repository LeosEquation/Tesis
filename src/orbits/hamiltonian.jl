
function HOrbitContinuation(f!, bc!, hp_ini::HopfPoint{U}, params,
                            Δs::U, maxsteps::T, tol::U, ite::T, Npoints::T,
                            integtol::U, integorder::T;
                            integmaxsteps::T = 500, parse_eqs::Bool = true) where {U<:Real, T<:Integer}

    @info "This function is running with $(Threads.nthreads()) threads."

    ###
    n = length(hp_ini.x) + 1
    m = Npoints
    N = m*(n-2) + 2

    Δτ = one(U) / (m - 1)
    τ = LinRange(zero(U), one(U), Npoints)

    # return τ, Δτ

    ###
    x = Array{U, 3}(undef, maxsteps, m, n)
    μ = Array{Complex{U}, 2}(undef, maxsteps, n-2)
    stbl = Array{Bool, 1}(undef, maxsteps)
    PDtest = Array{U, 1}(undef, maxsteps)
    LPCtest = Array{U, 1}(undef, maxsteps)

    ###
    xT = TaylorN.([j for i in 1:m-1, j in 1:n], order = 1)
    Mvec = [zeros(U, n-2, n-2) for i in 1:m-1]
    # In = LinearAlgebra.I(n-2)

    ###

    x0 = zeros(U, N)
    PDtest0 = 0.0
    LPCtest0 = 0.0
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

    # @show x0

    Φ0 = zeros(U, N)

    comp3 = 0
    for comp1 in 1:m 
        for comp2 in 1:n-2
            comp3 += 1
            Φ0[comp3] = cos(2π*τ[comp1]) * real(hp_ini.vector[comp2]) + sin(2π*τ[comp1]) * imag(hp_ini.vector[comp2])
        end
    end

    normalize!(Φ0)

    dx0 = [ ω0 * imag(hp_ini.vector) ; 0.0 ; 0.0 ]

    ###
    x1 = copy(x0)
    PDtest1 = 0.0
    LPCtest1 = 0.0
    Φ1 = copy(Φ0)
    # ω1 = copy(ω0)

    ###
    J = zeros(U, N, N)
    F = zeros(U, N)
    v = zeros(U, N)
    v[N] = one(U)

    ###

    comp3 = 0
    for comp1 in 1:m 
        for comp2 in 1:n-2
            comp3 += 1
            x[1, comp1, comp2] = x0[comp3]
        end
        x[1, comp1, n-1] = x0[N-1]
        x[1, comp1, n] = x0[N]
    end

    # @show x1

    comp3 = 0
    for comp1 in 1:m-1
        for comp2 in 1:n-2
            comp3 += 1
            xT[comp1, comp2][0][1] = x0[comp3]
            for comp4 in 1:n
                if comp2 == comp4
                    xT[comp1, comp2][1][comp4] = 1.0
                else
                    xT[comp1, comp2][1][comp4] = 0.0
                end
            end
        end
        xT[comp1, n-1][0][1] = x0[N-1]
        xT[comp1, n][0][1] = x0[N]
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

    caches = [TaylorIntegration.init_cache(zero(U), xT[1, :], integorder, f!, params; parse_eqs) for _ in 1:Threads.nthreads()]

    ###

    Threads.@threads for comp1 in 1:m-1
        tid = Threads.threadid()
        xTview = @view xT[comp1, :]
        taylorinteg_optim!(f!, bc!, xTview, 0.0, Δτ, integtol, caches[tid], params; maxsteps = integmaxsteps)
    end

    for comp1 in 1:m-1
        for comp2 in 1:n-2
            for comp3 in 1:n-2
                Mvec[comp1][comp2, comp3] = xT[comp1, comp2][1][comp3]
            end
        end
    end

    pS = pschur!(Mvec, :L)

    PDtest0 = prod(real(pS.values) .+ 1.0)
    PDtest[1] = PDtest0
    # LPCtest0 = 0.0

    μ[1, :] .= pS.values

    if all( abs.(μ[1, :]) .<= 1.0) && all( abs.( abs.(μ[1, :]) .- 1.0 ) .> tol)
        stbl[1] = false 
    else 
        stbl[1] = true 
    end

    ###
    i = 2

    while i <= maxsteps

        for comp1 in 1:N
            x1[comp1] = x0[comp1] + Δs * Φ0[comp1]
        end

        # @show x1

        # if g(x1, params, 0.0)

        comp3 = 0
        for comp1 in 1:m-1
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

        Threads.@threads for comp1 in 1:m-1
            tid = Threads.threadid()
            xTview = @view xT[comp1, :]
            taylorinteg_optim!(f!, bc!, xTview, 0.0, Δτ, integtol, caches[tid], params; maxsteps = integmaxsteps)
        end

        OrbitJacobian!(J, xT, dx0, Φ0, n, m, N)
        OrbitSystem!(F, xT, x1, x0, dx0, Φ0, Δs, n, m, N)

        j = 1

        while j <= ite && norm(F) > tol

            if abs(det(J)) == 0.0
                break
            end

            x1 .-= J \ F

            comp3 = 0
            for comp1 in 1:m-1
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

            Threads.@threads for comp1 in 1:m-1
                tid = Threads.threadid()
                xTview = @view xT[comp1, :]
                taylorinteg_optim!(f!, bc!, xTview, 0.0, Δτ, integtol, caches[tid], params; maxsteps = integmaxsteps)
            end        

            OrbitJacobian!(J, xT, dx0, Φ0, n, m, N)
            OrbitSystem!(F, xT, x1, x0, dx0, Φ0, Δs, n, m, N)

            j += 1

        end

        if norm(F) > tol
            @warn("System norm did not converge to tolerance. Stopping continuation.")
            break
        end

        comp3 = 0

        for comp1 in 1:m
            for comp2 in 1:n-2
                comp3 += 1
                x[i, comp1, comp2] = x1[comp3]
            end
            # x1[N-1] = 0.0
            x[i, comp1, n-1] = x1[N-1]
            x[i, comp1, n] = x1[N]
        end

        for comp1 in 1:m-1
            for comp2 in 1:n-2
                for comp3 in 1:n-2
                    Mvec[comp1][comp2, comp3] = xT[comp1, comp2][1][comp3]
                end
            end
        end
    
        pS = pschur!(Mvec, :L)
    
        PDtest1 = prod(real(pS.values) .+ 1.0)
        PDtest[i] = PDtest1
        

        μ[i, :] .= pS.values

        if all( abs.(μ[i, :]) .<= 1.0) && all( abs.( abs.(μ[i, :]) .- 1.0 ) .> tol)
            stbl[i] = false 
        else 
            stbl[i] = true 
        end

        if sign(PDtest1) != sign(PDtest0) && abs(PDtest1) > tol && abs(PDtest0) > tol
            println(" PD Bifurcation change is at index $i ")
        end

        Φ1 .= J \ v

        normalize!(Φ1)

        LPCtest1 = Φ1[n]
        LPCtest[i] = LPCtest1

        
        if sign(LPCtest1) != sign(LPCtest0) && abs(LPCtest1) > tol && abs(LPCtest0) > tol
            println(" LPC Bifurcation change is at index $i ")
        end

        x0 .= x1
        Φ0 .= Φ1
        PDtest0 = PDtest1
        LPCtest0 = LPCtest1

        for comp1 in 1:n-2
            x00[comp1] = x0[comp1]
        end

        x00[n-1] = x0[N-1]
        x00[n] = x0[N]

        f!(dx0, x00, params, 0.0)

        i += 1

        print(" \r i = $i \t \t ||F|| = $(norm(F)) \t \t  ")

    end

    println(" ")

    return PeriodicBranch(τ, x[1:i-1, :, :], μ[1:i-1, :], stbl[1:i-1])

end



















function HOrbitContinuation(f!, bcU!, bc!, x_ini::Array{U, 1}, params,
                            Δs::U, maxsteps::T, tol::U, ite::T, Npoints::T,
                            integtol::U, integorder::T;
                            integmaxsteps::T = 500, parse_eqs::Bool = true) where {U<:Real, T<:Integer}

    @info "This function is running with $(Threads.nthreads()) threads."

    ###
    n = length(x_ini)
    m = Npoints
    N = m*(n-2) + 2

    Δτ = one(U) / (m - 1)
    τ = LinRange(zero(U), one(U), Npoints)

    # return τ, Δτ

    ###
    x = Array{U, 3}(undef, maxsteps, m, n)
    μ = Array{Complex{U}, 2}(undef, maxsteps, n-2)
    stbl = Array{Bool, 1}(undef, maxsteps)
    PDtest = Array{U, 1}(undef, maxsteps)
    LPCtest = Array{U, 1}(undef, maxsteps)

    ###
    xT = TaylorN.([j for i in 1:m-1, j in 1:n], order = 1)
    Mvec = [zeros(U, n-2, n-2) for i in 1:m-1]
    # In = LinearAlgebra.I(n-2)

    ###

    x0 = zeros(U, N)
    PDtest0 = 0.0
    LPCtest0 = 0.0
    x00 = zeros(U, n)
    ω0 = 2π / x_ini[n]

    comp3 = 0
    for comp1 in 1:m
        for comp2 in 1:n-2
            comp3 += 1
            x0[comp3] = x_ini[comp2]
        end
    end
    
    x0[N-1] = x_ini[n-1]
    x0[N] = x_ini[n]

    Φ0 = zeros(U, N)

    dx0 = zeros(n)

    f!(dx0, x_ini, params, 0.0)

    ###
    x1 = copy(x0)
    PDtest1 = 0.0
    LPCtest1 = 0.0
    Φ1 = copy(Φ0)

    ###
    J = zeros(U, N, N)
    F = zeros(U, N)
    v = zeros(U, N)
    v[N] = one(U)

    ###

    x_temp = copy(x_ini)

    cache_ini = TaylorIntegration.init_cache(zero(U), x_ini, integorder, f!, params; parse_eqs)

    comp3 = n - 2

    for comp1 in 2:m-1

        taylorinteg_optim!(f!, bcU!, x_temp, 0.0, Δτ, integtol, cache_ini, params; maxsteps = integmaxsteps)

        for comp2 in 1:n-2

            comp3 += 1

            x0[comp3] = x_temp[comp2]
            x1[comp3] = x_temp[comp2]

        end

    end

    comp3 = 0

    for comp1 in 1:m
        for comp2 in 1:n-2
            comp3 += 1
            x[1, comp1, comp2] = x0[comp3]
        end
        # x1[N-1] = 0.0
        x[1, comp1, n-1] = x0[N-1]
        x[1, comp1, n] = x0[N]
    end

    comp3 = 0
    for comp1 in 1:m-1
        for comp2 in 1:n-2
            comp3 += 1
            xT[comp1, comp2][0][1] = x0[comp3]
            for comp4 in 1:n
                if comp2 == comp4
                    xT[comp1, comp2][1][comp4] = 1.0
                else
                    xT[comp1, comp2][1][comp4] = 0.0
                end
            end
        end
        xT[comp1, n-1][0][1] = x0[N-1]
        xT[comp1, n][0][1] = x0[N]
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

    caches = [TaylorIntegration.init_cache(zero(U), xT[1, :], integorder, f!, params; parse_eqs) for _ in 1:Threads.nthreads()]

    ###

    Threads.@threads for comp1 in 1:m-1
        tid = Threads.threadid()
        xTview = @view xT[comp1, :]
        taylorinteg_optim!(f!, bc!, xTview, 0.0, Δτ, integtol, caches[tid], params; maxsteps = integmaxsteps)
    end

    OrbitJacobian!(J, xT, dx0, Φ0, n, m, N)
    OrbitSystem!(F, xT, x1, x0, dx0, Φ0, Δs, n, m, N)

    # @show F

    # @show cond(J)

    NS = nullspace(J)

    if size(NS, 2) == 0
        error("Initial point does not have a valid branch direction. Cannot continue.")
    end

    if size(NS, 2) > 1
        @show NS
        error("Initial point is a codimension-2 bifurcation. Cannot continue.")
    end

    Φ0 .= NS
    # @show Φ0

    for comp1 in 1:m-1
        for comp2 in 1:n-2
            for comp3 in 1:n-2
                Mvec[comp1][comp2, comp3] = xT[comp1, comp2][1][comp3]
            end
        end
    end

    pS = pschur!(Mvec, :L)

    PDtest0 = prod(real(pS.values) .+ 1.0)
    PDtest[1] = PDtest0

    μ[1, :] .= pS.values

    if all( abs.(μ[1, :]) .<= 1.0) && all( abs.( abs.(μ[1, :]) .- 1.0 ) .> tol)
        stbl[1] = false 
    else 
        stbl[1] = true 
    end

    ###
    i = 2

    while i <= maxsteps

        for comp1 in 1:N
            x1[comp1] = x0[comp1] + Δs * Φ0[comp1]
        end

        # @show x1

        # if g(x1, params, 0.0)

        comp3 = 0
        for comp1 in 1:m-1
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

        Threads.@threads for comp1 in 1:m-1
            tid = Threads.threadid()
            xTview = @view xT[comp1, :]
            taylorinteg_optim!(f!, bc!, xTview, 0.0, Δτ, integtol, caches[tid], params; maxsteps = integmaxsteps)
        end

        OrbitJacobian!(J, xT, dx0, Φ0, n, m, N)
        OrbitSystem!(F, xT, x1, x0, dx0, Φ0, Δs, n, m, N)

        # @show cond(J)

        j = 1

        while j <= ite && norm(F) > tol

            if abs(det(J)) == 0.0
                break
            end

            x1 .-= J \ F

            # @show x1

            comp3 = 0
            for comp1 in 1:m-1
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

            # @show evaluate(xT)

            Threads.@threads for comp1 in 1:m-1
                tid = Threads.threadid()
                xTview = @view xT[comp1, :]
                taylorinteg_optim!(f!, bc!, xTview, 0.0, Δτ, integtol, caches[tid], params; maxsteps = integmaxsteps)
            end        

            OrbitJacobian!(J, xT, dx0, Φ0, n, m, N)
            OrbitSystem!(F, xT, x1, x0, dx0, Φ0, Δs, n, m, N)

            # @show cond(J)

            j += 1

        end

        if norm(F) > tol
            @warn("System norm did not converge to tolerance. Stopping continuation.")
            break
        end

        comp3 = 0

        for comp1 in 1:m
            for comp2 in 1:n-2
                comp3 += 1
                x[i, comp1, comp2] = x1[comp3]
            end
            # x1[N-1] = 0.0
            x[i, comp1, n-1] = x1[N-1]
            x[i, comp1, n] = x1[N]
        end

        for comp1 in 1:m-1
            for comp2 in 1:n-2
                for comp3 in 1:n-2
                    Mvec[comp1][comp2, comp3] = xT[comp1, comp2][1][comp3]
                end
            end
        end
    
        pS = pschur!(Mvec, :L)
    
        PDtest1 = prod(real(pS.values) .+ 1.0)
        PDtest[i] = PDtest1
        

        μ[i, :] .= pS.values

        if all( abs.(μ[i, :]) .<= 1.0) && all( abs.( abs.(μ[i, :]) .- 1.0 ) .> tol)
            stbl[i] = false 
        else 
            stbl[i] = true 
        end

        if sign(PDtest1) != sign(PDtest0) && abs(PDtest1) > tol && abs(PDtest0) > tol
            println(" PD Bifurcation change is at index $i ")
        end

        Φ1 .= J \ v

        normalize!(Φ1)

        LPCtest1 = Φ1[n]
        LPCtest[i] = LPCtest1

        
        # if sign(LPCtest1) != sign(LPCtest0) && abs(LPCtest1) > tol && abs(LPCtest0) > tol
        #     println(" LPC Bifurcation change is at index $i ")
        # end

        x0 .= x1
        Φ0 .= Φ1
        PDtest0 = PDtest1
        LPCtest0 = LPCtest1

        for comp1 in 1:n-2
            x00[comp1] = x0[comp1]
        end

        x00[n-1] = x0[N-1]
        x00[n] = x0[N]

        f!(dx0, x00, params, 0.0)

        i += 1

        print(" \r i = $i \t \t ||F|| = $(norm(F)) \t \t  ")

    end

    println(" ")

    return PeriodicBranch(τ, x[1:i-1, :, :], μ[1:i-1, :], stbl[1:i-1])

end