

function _NewtonRaphsonH!(H, E::Float64, x::Vector{Float64}, params, tol::Float64, max_iter::Int64, 
                         S::Vector{Taylor1{Float64}}, f::Taylor1{Float64}, df::Taylor1{Float64})
    
    iter_count = 1
    f .= H(x + S, params) - E
    df .= differentiate(f)
    
    while iter_count <= max_iter && abs(evaluate(f)) > tol
        @inbounds x[3] -= evaluate(f / df)
        f .= H(x + S, params) - E
        df .= differentiate(f)
        iter_count += 1
    end

end

function AdjustEnergy(H, E::Float64, Q1::Float64, P1::Float64, Q2lims, params, tol::Float64, ite::Int64; nrdivs::Int64 = 100)

    x_energy = Vector{Vector{Float64}}()

    S = zeros(Taylor1{Float64}, 4)
    S[3] = Taylor1(Float64, 1)

    f  = Taylor1(Float64, 1)
    df = Taylor1(Float64, 1)

    x = zeros(Float64, 4)
    @inbounds x[1] = Q1
    @inbounds x[2] = P1
    @inbounds x[4] = P2section

    x_temp = zeros(Float64, 4)
    @inbounds x_temp[1] = Q1
    @inbounds x_temp[2] = P1
    @inbounds x_temp[4] = P2section

    Q2min, Q2max = Q2lims(x, params)
    ΔQ2 = Q2max - Q2min

    δvar = ΔQ2 / nrdivs

    x[3] = δvar

    E_old = H(x, params)

    for _ in 1:(nrdivs - 2)

        x[3] += δvar

        E_new = H(x, params)

        if (E_old - E) * (E_new - E) <= 0.0

            x_temp .= x

            _NewtonRaphsonH!(H, E, x_temp, params, tol, ite, S, f, df)

            if abs(H(x_temp, params) - E) <= tol

                # @show(x_temp)

                push!(x_energy, copy(x_temp))

            end
        end

    E_old = E_new

    end

    return x_energy

end

function PoincareSection(H, E::Float64, nrand::Int64, Q1min::Float64, Q1max::Float64, Q2lims, 
                         P1min::Union{Float64, Irrational}, P1max::Union{Float64, Irrational}, P2section::Union{Float64, Irrational}, 
                         params, tol::Float64, ite::Int64,
                         f!, g, h!, tmax::Float64, integorder::Int64, integtol::Float64, maxsteps::Int64;
                         nrdivs::Int64 = 100, prev_message::String = "")
    

    # poincare_map =  Vector{Matrix{Float64}}(undef, nrand)

    poincare_map = Vector{Matrix{Float64}}()

    # M = Matrix{Float64}(undef, 0, 2)

    S = zeros(Taylor1{Float64}, 4)
    S[3] = Taylor1(Float64, 1)

    f  = Taylor1(Float64, 1)
    df = Taylor1(Float64, 1)

    x = Vector{Float64}(undef, 4)
    @inbounds x[4] = P2section

    x_temp = Vector{Float64}(undef, 4)
    @inbounds x_temp[4] = P2section

    ΔQ1 = Q1max - Q1min
    ΔP1 = P1max - P1min

    Q2min = 0.0
    Q2max = 0.0

    E_new = 0.0
    E_old = 0.0

    δvar = 0.0

    # rootfind = false

    # i = 1

    for i in 1:nrand

        # rootfind = false

        @inbounds x[1] = Q1min + rand() * ΔQ1

        @inbounds x[2] = P1min + rand() * ΔP1

        # @show(x[1:2])

        Q2min, Q2max = Q2lims(x, params)

        δvar = (Q2max - Q2min) / nrdivs

        @inbounds x[3] = Q2min + δvar

        # @show(x)

        E_old = H(x, params)
        
        for _ in 1:(nrdivs - 2)

            @inbounds x[3] += δvar

            E_new = H(x, params)
            
            if (E_old - E) * (E_new - E) <= 0.0

                x_temp .= x

                _NewtonRaphsonH!(H, E, x_temp, params, tol, ite, S, f, df)
                
                if abs(H(x_temp, params) - E) <= tol

                    sol = taylorinteg_wrap(f!, g, h!, x_temp, 0.0, tmax, integorder, integtol, params, maxsteps = maxsteps, show_warn = false)
                    
                    # @inbounds poincare_map[i] = sol.xevents[:, 1:2]

                    # @show sol.xevents[:, 3:4]
                    # @show H(sol.xevents[1, 1:4], params)
                    # @show H(sol.xevents[end, 1:4], params)
                    # @show x_temp[4]

                    push!(poincare_map, copy(sol.xevents[:, 1:2]))

                    # rootfind = true

                    sol = nothing

                    # break

                end

            end
            
            E_old = E_new

        end

        # if !(rootfind)
        #     @inbounds poincare_map[i] = M
        # end

        print("\r $prev_message Progreso = $(round(i*100/nrand, digits = 2)) % ")

    end

    return vcat(poincare_map...)

end

