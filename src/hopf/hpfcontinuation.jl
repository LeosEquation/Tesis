function HPContinuation(f!, hp_ini::HopfPoint{U}, params, pmin::Array{U, 1}, pmax::Array{U, 1},
                        Δs::U, maxsteps::T, tol::U, ite::T) where {U<:Real, T<:Number}

#-

    n = length(hp_ini.x) # Dimensión del sistema

    x = Array{U,2}(undef, maxsteps, n)
    ω = Array{U,1}(undef, maxsteps)

###

# Inicializando TaylorN

    set_variables("δx", numvars = n, order = 2)

    δx = TaylorN.(1:n, order = 2)

    xaux = copy(δx)
    dx = zero(δx)

    ###

    # Inicializando valores previos

    q0 = zeros(U, 2*n - 1)
    q0[1:n] .= copy(hp_ini.x)
    q0[n+1:2*n-2] .= real(hp_ini.vector) / norm(real(hp_ini.vector))
    q0[2*n-1] = imag(hp_ini.value)^2

    ###

    A = zeros(U, n-1, n-2)

    q1 = copy(q0)
    w = zeros(U, n-2)

    ###

    # Inicializando variables para la continuación

    J = zeros(U, 2*n - 1, 2*n - 1)
    F = zeros(U, 2*n - 1)
    Φ = zeros(U, 2*n - 1)
    v = zeros(U, 2*n - 1)
    v[2*n - 1] = 1.0

    xzero = zero(hp_ini.x)

    dxeval = zeros(U, n)
    Jeval = zeros(U, n, n)

    ###

    # Iniciamos el primer paso

    for j in 1:n
        x[1, j] = q0[j]
    end

    for j in 1:n
        xaux[j][0][1] = q0[j]
    end

    f!(dx, xaux, params, zero(U))

    TaylorSeries.evaluate!(dx, xzero, dxeval)
    TaylorSeries.jacobian!(Jeval, dx)

    for i in 1:n-2
        for j in 1:n-2
            A[i,j] = sum(Jeval[j,k] * Jeval[k,i] for k in 1:n-2)
        end
        A[i, i] += q0[2*n-1]
        A[n-1, i] = q0[n+i]
    end

    _, S, Vt = svd(A)
    w .= Vt[:, end]    # Último vector columna
    normalize!(w)

    # @show S[end]

    # w .= nullspace(A)

    HopfJacobian!(J, Jeval, dx, q0, w, Φ, n)
    # HopfSystem!(F, dxeval, Jeval, q0, q0, Φ, w, Δs, n)

    # @show F

    NS = nullspace(J)

    if size(NS, 2) == 0
    @show(J)
    error("No hay dirección a seguir")
    end

    if size(NS, 2) > 1
    display(NS)
    error("Hay más de una dirección inicial a seguir")
    end

    Φ .= NS

    # @show Φ

    for j in 1:n
        x[1, j] = q0[j]
    end

    ω[1] = sqrt(q0[2*n-1])

#

    i = 2

    while i <= maxsteps && pmin[1] <= q1[n-1] <= pmax[1] && pmin[2] <= q1[n] <= pmax[2]

        for j in 1:2*n-1
            q1[j] = q0[j] + Φ[j] * Δs
        end

        for j in 1:n
            xaux[j][0][1] = q1[j]
        end
    
        f!(dx, xaux, params, zero(U))
    
        TaylorSeries.evaluate!(dx, xzero, dxeval)
        TaylorSeries.jacobian!(Jeval, dx)
    
        HopfJacobian!(J, Jeval, dx, q1, w, Φ, n)
        HopfSystem!(F, dxeval, Jeval, q1, q0, w, Φ, Δs, n)

        j = 1

        # @show cond(J)

        # println("ite = $j, ||F|| = $(norm(F))")

        while j <= ite && norm(F) > tol

            if abs(det(J)) == 0.0
                break
            end

            q1 .-= J \ F
            # println("‖δq‖ = ", norm(δq))
            # q1 .-= δq

            for k in 1:n
                xaux[k][0][1] = q1[k]
            end
        
            f!(dx, xaux, params, zero(U))
        
            TaylorSeries.evaluate!(dx, xzero, dxeval)
            TaylorSeries.jacobian!(Jeval, dx)
        
            HopfJacobian!(J, Jeval, dx, q1, w, Φ, n)
            HopfSystem!(F, dxeval, Jeval, q1, q0, w, Φ, Δs, n)

            j += 1

            # @show cond(J)

            # println("ite = $j, ||F|| = $(norm(F))")

        end

        if norm(F) >= tol
            @warn("La norma del sistema no convergió. Abortando.")
            break
        end

        if q1[2*n-1] <= tol
            @warn("La parte imaginaria se anuló. Abortando")
            break
        end

        for j in 1:n
            x[i, j] = q1[j]
        end

        ω[i] = sqrt(q1[2*n-1])

        for j in 1:n-2
            for k in 1:n-2
                A[j,k] = sum(Jeval[k,l] * Jeval[l,j] for l in 1:n-2)
            end
            A[j, j] += q1[2*n-1]
            A[n-1, j] = q1[n+j]
        end

        _, S, Vt = svd(A)
        w .= Vt[:, end]    # Último vector columna
        normalize!(w)
        # @show S[end]

        # w .= nullspace(A)

        Φ .= J \ v

        normalize!(Φ)

        # print("\r ite = $i : norm = $(norm(F)) \t")

        q0 .= q1

        i += 1

    end

    return x[1:i-1,:], ω[1:i-1]

end