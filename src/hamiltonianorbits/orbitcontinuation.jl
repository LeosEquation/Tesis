

function PeriodicJacobianH!(J::Matrix{Float64}, H, sol::Vector{TaylorN{Float64}}, xₛ::Vector{Float64})

    @views begin

        for j in 1:5

            J[1, j] = evaluate(differentiate(sol[1], j))
            J[2, j] = evaluate(differentiate(sol[2], j))
            J[3, j] = evaluate(differentiate(sol[6], j))
            J[4, j] = evaluate(differentiate(H, j))
            J[5, j] = xₛ[j]

        end

        J[1, 1] -= 1.0
        J[2, 2] -= 1.0

    end

end



function PeriodicSystemH!(F::Vector{Float64}, H, E, sol::Vector{TaylorN{Float64}}, 
                          x1::Vector{Float64}, x0::Vector{Float64}, 
                          xₛ::Vector{Float64}, 
                          Δs::Float64)

    @views begin

        F[1] = evaluate(sol[1]) - x1[1]

        F[2] = evaluate(sol[2]) - x1[2]

        F[3] = evaluate(sol[6]) - x1[6]

        F[4] = evaluate(H) - E

        F[5] = dot(x1[1:5] - x0[1:5], xₛ) - Δs

    end

end


function PeriodicH(H, E, f!, h!, x_ini, params, Δs, maxsteps, tol, ite, integorder, integtol, integmaxsteps)

    #-

    x = Array{Float64,2}(undef, maxsteps, 6)

    μ = Array{ComplexF64, 2}(undef, maxsteps, 4)

    S = set_variables("s", numvars = 6, order = 1)

    x0 = copy(x_ini)
    dx0 = zero(x0)
    M0 = zeros(Float64, 4, 4)
    μ0 = zeros(ComplexF64, 4)
    v0 = zeros(ComplexF64, 4, 4)

    x1 = copy(x0)
    dx1 = zero(x1)
    M1 = zeros(Float64, 4, 4)
    μ1 = zeros(ComplexF64, 4)
    v1 = zeros(ComplexF64, 4, 4)

    F_ = zeros(Float64, 5)
    J_ = zeros(Float64, 5, 5)

    xₛ = zeros(Float64, 5)

    δv = zeros(Float64, 5)

    dir = [zeros(Float64, 4); 1.0]

    sol = zeros(TaylorN{Float64}, 4)

    #

    x[1, :] .= x0

    sol = taylorinteg_wrap(f!, h!, x1 + S, 0.0, 1.0, integorder, integtol, params; maxsteps = integmaxsteps).x[end, :]

    # f!(dx1, evaluate(sol), params, 1.0)
    # f!(dx0, x1, params, 0.0)

    H_TN = H(x1 + S, params)

    PeriodicJacobianH!(J_, H_TN, sol, xₛ)

    NS = nullspace(J_)

    if size(NS, 2) == 0
        @show(J_)
        error("No hay dirección a seguir")
    end

    if size(NS, 2) > 1
        @show(NS)
        error("Hay más de una dirección inicial a seguir")
    end

    xₛ .= NS

    PeriodicJacobianH!(J_, H_TN, sol, xₛ)
    PeriodicSystemH!(F_, H_TN, E, sol, x1, x0, xₛ, Δs)

    for i in 1:3
        for j in 1:3
            M0[i, j] = evaluate(TaylorSeries.differentiate(sol[i], j))
        end
        M0[4, i] = evaluate(TaylorSeries.differentiate(sol[6], i))
        M0[i, 4] = evaluate(TaylorSeries.differentiate(sol[i], 6))
    end
    M0[4, 4] = evaluate(TaylorSeries.differentiate(sol[6], 6))

    μ0 .= eigvals(M0)

    v0 .= eigvecs(M0)

    μ[1, :] .= μ0

    #

    i = 2

    while i <= maxsteps

        # println(" i = $i")

        x1[1:5] = x0[1:5] + (xₛ * Δs)

        h!(x1, params) 

        sol = taylorinteg_wrap(f!, h!, x1 + S, 0.0, 1.0, integorder, integtol, params; maxsteps = integmaxsteps).x[end, :]

        H_TN .= H(x1 + S, params)

        PeriodicJacobianH!(J_, H_TN, sol, xₛ)
        PeriodicSystemH!(F_, H_TN, E, sol, x1, x0, xₛ, Δs)

        j = 1

        while j <= ite && norm(F_) > tol

            # println(" j = $j : norm = $(norm(F_))")

            if abs(det(J_)) == 0.0
                print("\n")
                @show(det(J_))
                @warn("El jacobiano del sistema se volvió singular. Abortando.")
                return x[1:i-1, :], μ[1:i-1, :]
            end

            δv .= J_ \ F_

            x1[1:5] .-= δv

            if h!(x1, params)
                print("\n")
                @warn("Se ha salido del rango permitido. Abortando.")
                return x[1:i-1, :], μ[1:i-1, :]
            end

            sol = taylorinteg_wrap(f!, h!, x1 + S, 0.0, 1.0, integorder, integtol, params; maxsteps = integmaxsteps).x[end, :]

            H_TN .= H(x1 + S, params)

            PeriodicJacobianH!(J_, H_TN, sol, xₛ)
            PeriodicSystemH!(F_, H_TN, E, sol, x1, x0, xₛ, Δs)

            j += 1

        end

        if norm(F_) >= tol
            print("\n")
            @show(F_)
            @warn("El sistema periódico superó la tolerancia, abortando")
            return x[1:i-1, :], μ[1:i-1, :]
        end

        if abs(det(J_)) == 0.0
            print("\n")
            @show(det(J_))
            @warn("El jacobiano del sistema se volvió singular. Abortando.")
            return x[1:i-1, :], μ[1:i-1, :]
        end

        xₛ .= J_ \ dir

        xₛ .= xₛ / norm(xₛ)

        for k in 1:3
            for l in 1:3
                M1[k, l] = evaluate(TaylorSeries.differentiate(sol[k], l))
            end
            M1[4, k] = evaluate(TaylorSeries.differentiate(sol[6], k))
            M1[k, 4] = evaluate(TaylorSeries.differentiate(sol[k], 6))
        end
        M1[4, 4] = evaluate(TaylorSeries.differentiate(sol[6], 6))
    
        μ1 .= eigvals(M1)

        v1 .= eigvecs(M1)

        x[i, :] .= x1

        μ[i, :] .= μ1

        x0 .= x1
        μ0 .= μ1
        M0 .= M1

        print("\r i = $i : norm = $(norm(F_)) \t")

        i += 1

    end

    return x[1:i-1, :], μ[1:i-1, :]

end

