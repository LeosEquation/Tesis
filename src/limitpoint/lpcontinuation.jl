

function LPContinuation(f!, x_ini::Array{U, 1}, v_ini::Array{U, 1}, 
                        params, Δs::U, maxsteps::T, tol::U, ite::T) where {U<:Real, T<:Integer}

    #-

    n = length(x_ini)

    x = Array{U, 2}(undef, maxsteps, n)
    # vectors = Array{U,2}(undef, maxsteps, n-2)

    ###

    variable_names = [string("δx", TaylorSeries.subscriptify(i)) for i in 1:n]

    TaylorSeries.set_variables(U, variable_names, order = 2)

    δx = TaylorN.(1:n, order = 2)

    xaux = copy(δx)
    dx = zero(δx)

    ###

    # Inicializando valores previos

    q0 = zeros(U, 2*n - 2)
    q0[1:n] .= x_ini
    q0[n+1:2*n-2] .= v_ini

    ###

    # Inicializando valores posteriores

    q1 = zeros(U, 2*n - 2)
    q1[1:n] .= x_ini
    q1[n+1:2*n-2] .= v_ini

    ###

    dx = Array{TaylorN{U},1}(undef, n)

    ###

    # Inicializando variables para la continuación

    J = zeros(U, 2*n - 2, 2*n - 2)
    F = zeros(U, 2*n - 2)
    Φ = zeros(U, 2*n - 2)
    v = zeros(U, 2*n - 2)
    v[2*n - 2] = one(U)

    Jeval = zeros(U, n, n)
    dxeval = zeros(U, n)
    xzero = zeros(U, n)

    ###

    # Iniciamos el primer paso

    x[1, :] .= q0[1:n]

    for j in 1:n
        xaux[j][0][1] = q0[j]
    end

    f!(dx, xaux, params, zero(U))

    TaylorSeries.jacobian!(Jeval, dx)
    
    LPJacobian!(J, Jeval, dx, q0, Φ, n)

    NS = nullspace(J)

    if size(NS, 2) == 0
        @show(J)
        error("No hay dirección a seguir")
    end

    if size(NS, 2) > 1
        display(NS)
        error("El espacio nulo del jacobiano es de dimensión 2")
    end

    Φ .= NS

    #

    i = 2

    while i <= maxsteps

        for j in 1:2*n-2
            q1[j] = q0[j] + Φ[j] * Δs
        end

        for j in 1:n
            xaux[j][0][1] = q1[j]
        end

        f!(dx, xaux, params, 0.0)

        TaylorSeries.evaluate!(dx, xzero, dxeval)
        TaylorSeries.jacobian!(Jeval, dx)

        LPJacobian!(J, Jeval, dx, q1, Φ, n)
        LPSystem!(F, Jeval, dxeval, q1, q0, Φ, Δs, n)

        # @show F

        j = 1

        while j <= ite && norm(F) > tol

            # println(" ite = $j, ||F|| = $(norm(F))")

            if abs(det(J)) == 0.0
                break
            end

            q1 .-= J \ F
    
            for k in 1:n
                xaux[k][0][1] = q1[k]
            end
    
            f!(dx, xaux, params, 0.0)
    
            TaylorSeries.evaluate!(dx, xzero, dxeval)
            TaylorSeries.jacobian!(Jeval, dx)
    
            LPJacobian!(J, Jeval, dx, q1, Φ, n)
            LPSystem!(F, Jeval, dxeval, q1, q0, Φ, Δs, n)

            j += 1

        end

        if norm(F) >= tol
            @warn("La norma del sistema no convergió. Abortando.")
            break
        end

        x[i, :] .= q1[1:n]

        Φ .= J \ v

        normalize!(Φ)

        # print("\r i = $i : norm = $(norm(F)) \t")

        q0 .= q1

        i += 1

    end

    return x[1:i-1, :]

end