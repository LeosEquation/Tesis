function PoincareSection!(
        plt,
        f!,
        bc!,
        g,
        H,
        rs!,
        nrfind,
        E::U,
        dim::T,
        params,
        tmax::U,
        abstol::U,
        order::T,
        x_idx::T,
        y_idx::T,
        nrand::T;
        parse_eqs::Bool= true,
        maxsteps::Int= 500,
        dense::Bool= true,
        ms::U= one(U),
        color = :black,
        eventorder::Int = 0,
        newtoniter::Int = 10,
        nrabstol::U = eps(U),
    ) where {U<:Real, T<:Int}

    ξ = zeros(dim)
    dξ = zeros(dim)

    rs!(ξ, params, 0.0)

    idx, Δξ, ξ_max, ξ_min = nrfind(ξ, params, 0.0)
    
    ξaux = [Taylor1(0.0, 1) for i in 1:dim]

    ξaux[idx][1] = 1.0

    Haux = Taylor1(0.0, 1) 

    cache = TaylorIntegration.init_cache(Val(dense), 0.0, ξ, maxsteps, order, f!, params; parse_eqs)

    npoints = 0

    for i in 1:nrand
        
        rs!(ξ, params, 0.0)

        idx, Δξ, ξ_max = nrfind(ξ, params, 0.0)
    
        ξ_old = Δξ 
    
        ξ_new = ξ_old + Δξ
    
        ξ[idx] = ξ_old
        
        nE = 0
        
        if ξ[idx] < ξ_max
            
            E_old = H(ξ, params, 0.0)
            
            while ξ_new < ξ_max
        
                ξ[idx] = ξ_new
                
                E_new = H(ξ, params, 0.0)
        
                if sign(E_old - E) != sign(E_new - E)

                    for j in 1:dim
                        ξaux[j][0] = ξ[j]
                    end

                    niter = 1

                    Haux .= H(ξaux, params, 0.0)

                    while niter < newtoniter && abs( Haux[0] - E ) > nrabstol
                        
                        ξaux[idx][0] -= ( Haux[0] - E ) / Haux[1]
                        Haux .= H(ξaux, params, 0.0)

                        niter += 1

                    end

                    if abs( Haux[0] - E ) > nrabstol
                        break
                    end

                    ξ[idx] = ξaux[idx][0]

                    f!(dξ, ξ, params, 0.0)

                    if dξ[idx] != 0
    
                        sol = TaylorIntegration.taylorinteg_wrap!(Val(dense), f!, bc!, g, ξ, 0.0, tmax, abstol, cache, params; maxsteps, eventorder, newtoniter, nrabstol)
                        scatter!(plt, sol.xevents[:, x_idx], sol.xevents[:, y_idx], ms = ms, color = color, label = "")
        
                        nE += 1

                        break

                    end
                        
                    # break
                    
                end
        
                ξ_old = ξ_new
        
                E_old = E_new
        
                ξ_new = ξ_old + Δξ
                
            end

        end

        if nE >= 1
            npoints += 1
        end

        if i % ceil(Int64, nrand / 100) == 0
            print("\r Progress : $( round( i / nrand , digits = 2 ) ) \t ")
        end
    end

    print("\n")

    println( " El área con soluciones es aproximadamente $(npoints / nrand) " )

end

function PPoincareSection!(
        plt,
        f!,
        bc!,
        g,
        H,
        rs!,
        nrfind,
        E::U,
        dim::T,
        params,
        tmax::U,
        abstol::U,
        order::T,
        x_idx::T,
        y_idx::T,
        nrand::T;
        parse_eqs::Bool= true,
        maxsteps::Int= 500,
        dense::Bool= true,
        ms::U= one(U),
        color = :black,
        eventorder::Int = 0,
        newtoniter::Int = 10,
        nrabstol::U = eps(U),
    ) where {U<:Real, T<:Int}

    @info "This function is running with $(Threads.nthreads()) threads."

    nT = Threads.nthreads()
    
    ξ = zeros(dim)

    idx0, _, _ = nrfind(ξ, params, 0.0)

    rs!(ξ, params, 0.0)

    caches = [
                TaylorIntegration.init_cache(Val(dense), 0.0, ξ, maxsteps, order, f!, params; parse_eqs)
                for _ in 1:nT
             ]

    ξs = [ zeros(dim) for i in 1:nT ]

    ξauxs = [ [Taylor1(0.0, 1) for i in 1:dim] for i in 1:nT ]

    for i in 1:nT

        ξauxs[i][idx0][1] = 1.0

    end

    Hauxs = [ Taylor1(0.0, 1) for i in 1:nT ]

    ξ_olds = [ 0.0 for i in 1:nT ]

    ξ_news = [ 0.0 for i in 1:nT ]

    E_olds = [ 0.0 for i in 1:nT ]

    E_news = [ 0.0 for i in 1:nT ]

    niters = [ 0 for i in 1:nT ]

    progress = [ 0 for i in 1:nT ]

    @Threads.threads for i in 1:nrand

        tid = Threads.threadid()

        progress[tid] += 1
        
        rs!(ξs[tid], params, 0.0)

        idx, Δξ, ξ_max = nrfind(ξs[tid], params, 0.0)
    
        ξ_olds[tid] = Δξ 
    
        ξ_news[tid] = ξ_olds[tid] + Δξ
    
        ξs[tid][idx] = ξ_olds[tid]
        
        if ξs[tid][idx] < ξ_max
            
            E_olds[tid] = H(ξs[tid], params, 0.0)
            
            while ξ_news[tid] < ξ_max
        
                ξs[tid][idx] = ξ_news[tid]
                
                E_news[tid] = H(ξs[tid], params, 0.0)
        
                if sign(E_olds[tid] - E) != sign(E_news[tid] - E)

                    for j in 1:dim
                        ξauxs[tid][j][0] = ξs[tid][j]
                    end

                    niters[tid] = 1

                    Hauxs[tid] .= H(ξauxs[tid], params, 0.0)

                    while niters[tid] < newtoniter && abs( Hauxs[tid][0] - E ) > nrabstol
                        
                        ξauxs[tid][idx][0] -= ( Hauxs[tid][0] - E ) / Hauxs[tid][1]
                        Hauxs[tid] .= H(ξauxs[tid], params, 0.0)

                        niters[tid] += 1
                        
                    end

                    if abs( Hauxs[tid][0] - E ) > nrabstol
                        break
                    end

                    ξs[tid][idx] = ξauxs[tid][idx][0]
    
                    sol = TaylorIntegration.taylorinteg_wrap!(Val(dense), f!, bc!, g, ξs[tid], 0.0, tmax, abstol, caches[tid], params; maxsteps, eventorder, newtoniter, nrabstol)

                    scatter!(plt, sol.xevents[:, x_idx], sol.xevents[:, y_idx], ms = ms, color = color, label = "")

                    # sol = nothing
                    
                end
        
                ξ_olds[tid] = ξ_news[tid]
        
                E_olds[tid] = E_news[tid]
        
                ξ_news[tid] = ξ_olds[tid] + Δξ
                
            end

        end

        print("\r Progress : $( round( sum(progress) / nrand , digits = 2 ) ) \t ")

    end

    # return vcat([vcat(i...) for i in section]...)

end