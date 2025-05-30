

function LPContinuation(f!, g, lp_ini::LimitPoint{U}, 
                        params,
                        Δs::U, maxsteps::T, tol::U, ite::T) where {U<:Real, T<:Integer}

    #-

    n = length(lp_ini.x)
    ordtup = [ntuple(k -> count(==(k), (i, j)), TaylorSeries.get_numvars()) for i in 1:n, j in 1:n]

    x = Array{U, 2}(undef, maxsteps, n)


    ###

    # variable_names = [string("δx", TaylorSeries.subscriptify(i)) for i in 1:n]

    # TaylorSeries.set_variables(U, variable_names, order = 2)

    δx = TaylorN.(1:n, order = 2)

    xaux = copy(δx)
    dx = zero(δx)

    Jeval = zeros(U, n - 2, n)
    dxeval = zeros(U, n - 2)

    ###

    # Inicializando valores previos

    q0 = zeros(U, 2*n - 2)
    q0[1:n] .= lp_ini.x
    q0[n+1:2*n-2] .= lp_ini.nullvec

    ###

    # Inicializando valores posteriores

    q1 = copy(q0)

    ###

    # dx = Array{TaylorN{U},1}(undef, n)

    ###

    # Inicializando variables para la continuación

    J = zeros(U, 2*n - 2, 2*n - 2)
    F = zeros(U, 2*n - 2)
    Φ = zeros(U, 2*n - 2)
    v = zeros(U, 2*n - 2)
    v[2*n - 2] = one(U)

    

    ###

    # Iniciamos el primer paso

    x[1, :] .= q0[1:n]

    for j in 1:n
        xaux[j][0][1] = q0[j]
    end

    f!(dx, xaux, params, zero(U))

    for comp1 in 1:n-2
        for comp2 in 1:n
            Jeval[comp1,comp2] = dx[comp1][1][comp2]
        end
        dxeval[comp1] = dx[comp1][0][1]
    end
    
    LPJacobian!(J, Jeval, dx, q0, Φ, ordtup, n)
    LPSystem!(F, Jeval, dxeval, q1, q0, Φ, Δs, n)

    # @show F

    NS = nullspace(J)

    if size(NS, 2) == 0
        # @show(J)
        error("The jacobian in the point has not nullspace. Exiting.")
    end

    if size(NS, 2) > 1
        display(NS)
        @show J 
        @show F
        println( "corank J = $(size(Jeval, 2) - rank(Jeval))" )
        error("The jacobian in the point is 2 or more corank. Exiting.")
    end

    Φ .= NS

    #

    i = 2

    while i <= maxsteps ## && pmin[1] <= q1[n-1] <= pmax[1] && pmin[2] <= q1[n] <= pmax[2]

        # print("$i.")

        for j in 1:2*n-2
            q1[j] = q0[j] + Φ[j] * Δs
        end

        if g(q1, params, 0.0)
            break
        end

        for j in 1:n
            xaux[j][0][1] = q1[j]
        end

        f!(dx, xaux, params, 0.0)

        for comp1 in 1:n-2
            for comp2 in 1:n
                Jeval[comp1,comp2] = dx[comp1][1][comp2]
            end
            dxeval[comp1] = dx[comp1][0][1]
        end

        # TaylorSeries.evaluate!(dx, xzero, dxeval)
        # TaylorSeries.jacobian!(Jeval, dx)

        LPJacobian!(J, Jeval, dx, q1, Φ, ordtup, n)
        LPSystem!(F, Jeval, dxeval, q1, q0, Φ, Δs, n)

        j = 1

        while j <= ite && norm(F) > tol

            if abs(det(J)) == 0.0
                break
            end

            q1 .-= J \ F
        
            if g(q1, params, 0.0)
                break
            end

            for k in 1:n
                xaux[k][0][1] = q1[k]
            end
    
            f!(dx, xaux, params, 0.0)
    
            for comp1 in 1:n-2
                for comp2 in 1:n
                    Jeval[comp1,comp2] = dx[comp1][1][comp2]
                end
                dxeval[comp1] = dx[comp1][0][1]
            end
    
            LPJacobian!(J, Jeval, dx, q1, Φ, ordtup, n)
            LPSystem!(F, Jeval, dxeval, q1, q0, Φ, Δs, n)

            j += 1

        end

        if j > ite
            @warn("Tolerance branch was exceded. Exiting. ")
            break
        end

        x[i, :] .= q1[1:n]

        Φ .= J \ v

        normalize!(Φ)

        q0 .= q1

        i += 1

    end

    return x[1:i-1, :]

end