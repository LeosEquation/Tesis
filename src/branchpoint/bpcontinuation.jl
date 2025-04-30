

function BPContinuation(f!, bp_ini::BranchPoint{U}, params, pmin::Array{U, 1}, pmax::Array{U, 1},
                        Δs::U, maxsteps::T, tol::U, ite::T) where {U<:Real, T<:Integer}

#-

    n = length(bp_ini.x) # Dimensión del sistema
    ordtup = [ntuple(k -> count(==(k), (i, j)), TaylorSeries.get_numvars()) for i in 1:n, j in 1:n]

    x = Array{Float64,2}(undef, maxsteps, n) # Rama de equilibrio

    # variable_names = [string("δx", TaylorSeries.subscriptify(i)) for i in 1:n]

    # TaylorSeries.set_variables(U, variable_names, order = 2)

    δx = TaylorN.(1:n, order = 2)

    xaux = copy(δx)
    dx = zero(δx)

    q0 = zeros(U, 2*n - 1)
    q0[1:n] .= bp_ini.x
    q0[n+1:2*n-2] .= bp_ini.vector

    q1 = zeros(U, 2*n - 1)
    q1[1:n] .= bp_ini.x
    q1[n+1:2*n-2] .= bp_ini.vector

    dx = Array{TaylorN{U},1}(undef, n)

    J = zeros(U, 2*n - 1, 2*n - 1)
    F = zeros(U, 2*n - 1)
    Φ = zeros(U, 2*n - 1)
    v = zeros(U, 2*n - 1)
    v[2*n - 1] = one(U)

    Jeval = zeros(U, n, n)
    dxeval = zeros(U, n)
    xzero = zeros(U, n) 
    
    ###

    for j in 1:n
        x[1, j] = q0[j]
    end

    for j in 1:n
        xaux[j][0][1] = q0[j]
    end

    f!(dx, xaux, params, zero(U))

    TaylorSeries.jacobian!(Jeval, dx)
    TaylorSeries.evaluate!(dx, xzero, dxeval)

    BPJacobian!(J, Jeval, q0, dx, Φ, ordtup, n)
    BPSystem!(F, dxeval, Jeval, q1, q0, Φ, Δs, n)

    # @show J
    # @show F

    NS = nullspace(J)

    if size(NS, 2) == 0
        error("No hay dirección a seguir")
    end

    if size(NS, 2) > 1
        error("Hay más de una dirección inicial a seguir")
    end

    Φ .= NS

    ###

    i = 2

    while i <= maxsteps && pmin[1] <= q1[n-1] <= pmax[1] && pmin[2] <= q1[n] <= pmax[2]

        for j in 1:2*n-1
            q1[j] = q0[j] + Φ[j] * Δs
        end

        for j in 1:n
            xaux[j][0][1] = q1[j]
        end

        f!(dx, xaux, params, zero(U))

        TaylorSeries.jacobian!(Jeval, dx)
        TaylorSeries.evaluate!(dx, xzero, dxeval)

        BPJacobian!(J, Jeval, q1, dx, Φ, ordtup, n)
        BPSystem!(F, dxeval, Jeval, q1, q0, Φ, Δs, n)

        j = 1

        while j <= ite && norm(F) > tol

            if abs(det(J)) == 0.0
                break
            end

            q1 .-= J \ F

            for k in 1:n
                xaux[k][0][1] = q1[k]
            end

            f!(dx, xaux, params, zero(U))

            TaylorSeries.jacobian!(Jeval, dx)
            TaylorSeries.evaluate!(dx, xzero, dxeval)

            BPJacobian!(J, Jeval, q1, dx, Φ, ordtup, n)
            BPSystem!(F, dxeval, Jeval, q1, q0, Φ, Δs, n)

            j += 1

        end

        if norm(F) >= tol
            @warn("La norma del sistema no convergió. Abortando.")
            break
        end

        for j in 1:n
            x[i, j] = q1[j]
        end

        Φ .= J \ v

        normalize!(Φ)

        q0 .= q1

        i += 1

    end

    return x[1:i-1, :]

end