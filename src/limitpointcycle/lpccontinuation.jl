

function LPCContinuation(f!, x_ini::Vector{Float64}, v_ini::Vector{Float64}, params, Δs::Float64, maxsteps::Int64, tol::Float64, ite::Int64)

    #-

    n = length(v_ini) # Dimensión del sistema

    x = Array{Float64,2}(undef, maxsteps, n+2) # Rama de equilibrio

    ###

    # Inicializando TaylorN

    δx = set_variables("δx", numvars = n+2, order = 2)

    ###

    # Inicializando valores previos

    x0 = zeros(Float64, 2*n + 2)
    x0[1:n+2] .= copy(x_ini)
    x0[n+3:2*n+2] .= copy(v_ini)

    ###

    # Inicializando valores posteriores

    x1 = zeros(Float64, 2*n + 2)
    x1[1:n+2] .= copy(x_ini)
    x1[n+3:2*n+2] .= copy(v_ini)

    ###

    # Inicializando variables para bifurcaciones

    dx = Array{TaylorN{Float64},1}(undef, n+2)
    Fx = Array{TaylorN{Float64},2}(undef, n, n)

    ###

    # Inicializando valores reutilizables

    λ = zeros(ComplexF64, n)

    ###

    # Inicializando variables para la continuación

    J = zeros(Float64, 2*n + 2, 2*n + 2)
    F = zeros(Float64, 2*n + 2)
    Φ = zeros(Float64, 2*n + 2)
    v = zeros(Float64, 2*n + 2)
    v[2*n + 2] = 1.0

    ###

    # Iniciamos el primer paso

    x[1, :] .= x0[1:n+2]

    f!(dx, x0[1:n+2] + δx, params, 0.0)

    for j in 1:n
        Fx[:, j] .= differentiate.(dx[1:n], j)
    end

    # return dx
    # LPSystem!(F, Fx, dx, x0, x1, Φ, Δs, n)
    LPJacobian!(J, Fx, dx, x0, x1, Φ, n)

    # return J

    # @show F

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
        # println("$i")

        x1 .= x0 .+ (Φ * Δs)

        f!(dx, x1[1:n+2] + δx, params, 0.0)
        
        for j in 1:n
            Fx[:, j] .= differentiate.(dx[1:n], j)
        end

        LPJacobian!(J, Fx, dx, x0, x1, Φ, n)
        LPSystem!(F, Fx, dx, x0, x1, Φ, Δs, n)

        j = 1

        while j <= ite && norm(F) > tol

            # println("$j : norm = $(norm(F))")

            if abs(det(J)) == 0.0
                break
            end

            x1 .-= J \ F
    
            f!(dx, x1[1:n+2] + δx, params, 0.0)
            for j in 1:n
                Fx[:, j] .= differentiate.(dx[1:n], j)
            end

            LPJacobian!(J, Fx, dx, x0, x1, Φ, n)
            LPSystem!(F, Fx, dx, x0, x1, Φ, Δs, n)

            j += 1

        end

        if norm(F) >= tol
            @warn("La norma del sistema no convergió. Abortando.")
            break
        end

        x[i, :] .= x1[1:n+2]

        # @show J

        Φ .= J \ v

        normalize!(Φ)

        print("\r i = $i : norm = $(norm(F)) \t")

        x0 .= x1

        i += 1

    end

    return x[1:i-1, :]

end