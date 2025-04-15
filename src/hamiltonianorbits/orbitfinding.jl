
function OrbitJacobianH!(J::Matrix{Float64}, sol::TaylorSolution, dx1::Vector{Float64}, dx0::Vector{Float64})

    @views begin

        J[1, 1] = evaluate(differentiate(sol.x[end, 1], 1)) - 1.0
        J[1, 2] = evaluate(differentiate(sol.x[end, 1], 2))
        J[1, 3] = evaluate(differentiate(sol.x[end, 1], 3))
        J[1, 4] = dx1[1]

        J[2, 1] = evaluate(differentiate(sol.x[end, 2], 1))
        J[2, 2] = evaluate(differentiate(sol.x[end, 2], 2)) - 1.0
        J[2, 3] = evaluate(differentiate(sol.x[end, 2], 3))
        J[2, 4] = dx1[2]

        J[3, 1] = evaluate(differentiate(sol.x[end, 4], 1))
        J[3, 2] = evaluate(differentiate(sol.x[end, 4], 2))
        J[3, 3] = evaluate(differentiate(sol.x[end, 4], 3))
        J[3, 4] = dx1[4]

        J[4, 1] =  dx0[2]
        J[4, 2] = -dx0[1]
        J[4, 3] =  dx0[4]
        # J[4, 5] = 0.0

    end

end


function OrbitSystemH!(F::Vector{Float64}, H, E, sol::TaylorSolution, x1::Vector{Float64}, params)

    @views begin

        F[1] = evaluate(sol.x[end, 1]) - x1[1]

        F[2] = evaluate(sol.x[end, 2]) - x1[2]

        F[3] = evaluate(sol.x[end, 4]) - x1[4]

        F[4] = H(x1, params) - E

    end

end



function FindPeriodicSolution(H, E::Float64, f!, h!, x_ini::Vector{Float64}, params, T_ini::Float64, 
                             integorder::Int64, integtol::Float64, integmaxst::Int64, 
                             tol::Float64, ite::Int64)
    
    # setprecision(BigFloat, 80)

    J = zeros(Float64, 4, 4)
    F = zeros(Float64, 4)
    x = copy(x_ini)
    T = copy(T_ini)
    dx0 = zeros(Float64, 4)
    dx1 = zeros(Float64, 4)

    S = set_variables(Float64, "s", numvars = 4, order = 1)

    sol = taylorinteg_wrap(f!, h!, x + S, 0.0, T, integorder, integtol, params; maxsteps = integmaxst)

    f!(dx1, evaluate(sol.x[end, :]), params, T)
    f!(dx0, x, params, 0.0)

    OrbitSystemH!(F, H, E, sol, x, params)
    OrbitJacobianH!(J, sol, dx1, dx0)

    δv = zeros(Float64, 4)

    j = 1

    while j <= ite && norm(F) > tol

        println(" j = $j : norm = $(norm(F))")

        if abs(det(J)) == 0.0
            print("\n")
            @show(det(J))
            @warn("El jacobiano del sistema se volvió singular. Abortando.")
            return x, T
        end

        δv .= J \ F

        x[1] -= δv[1]
        x[2] -= δv[2]
        x[3] -= δv[3]
        T    -= δv[4]

        if h!(x, params)
            print("\n")
            @warn("Se ha salido del rango permitido. Abortando.")
            return x, T
        end

        sol = taylorinteg_wrap(f!, h!, x + S, 0.0, T, integorder, integtol, params; maxsteps = integmaxst)

        f!(dx1, evaluate(sol.x[end, :]), params, T)
        f!(dx0, x, params, 0.0)

        OrbitSystemH!(F, H, E, sol, x, params)
        OrbitJacobianH!(J, sol, dx1, dx0)

        j += 1

    end

    if norm(F) >= tol
        print("\n")
        @show(F)
        @warn("El sistema periódico superó la tolerancia, abortando")
        return x, T
    end

    return x, T

end




function FindPeriodicSolution(H, E::Float64, f!, g, h!, x_ini::Vector{Float64}, params, Tmax::Float64, 
                              integorder::Int64, integtol::Float64, integmaxst::Int64, period_order::Int64, 
                              tol::Float64, ite::Int64)

    sol = taylorinteg_wrap(f!, g, h!, x_ini, 0.0, Tmax, integorder, integtol, params; maxsteps = integmaxst)

    S = zeros(Taylor1{Float64}, 5)
    S[3] = Taylor1(Float64, 1)

    f  = Taylor1(Float64, 1)
    df = Taylor1(Float64, 1)

    x_initial = copy(x_ini)

    i = findfirst(x -> x > tol, sol.tevents) - 1

    x_initial[1] = mean(sol.xevents[i+period_order:period_order:end, 1])
    x_initial[2] = mean(sol.xevents[i+period_order:period_order:end, 2])

    _NewtonRaphsonH!(H, E, x_initial, params, tol, ite, S, f, df)

    x_norm = norm(sol.xevents[i+period_order, :] - sol.x[1, :])

    print("\r j = 1 : norm(f) = $(x_norm)")

    j = 1

    while j <= ite && x_norm > tol

        sol = taylorinteg_wrap(f!, g, h!, x_initial, 0.0, Tmax, integorder, integtol, params; maxsteps = integmaxst)

        if length(sol.tevents) == 0
                @warn("No hay intersecciónes. Abortando.")
            break
        end

        i = findfirst(x -> x > tol, sol.tevents) - 1

        x_initial[1] = mean(sol.xevents[i+period_order:period_order:end, 1])
        x_initial[2] = mean(sol.xevents[i+period_order:period_order:end, 2])

        _NewtonRaphsonH!(H, E, x_initial, params, tol, ite, S, f, df)

        x_norm = norm(sol.xevents[i+period_order, :] - sol.x[1, :])

        j += 1

        print("\r j = $j : norm(f) = $(norm(x_initial - sol.x[1, :]))")

    end

    return sol, x_initial, sol.tevents[i+period_order]

end
