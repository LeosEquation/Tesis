function EnergyFinding(
        f!,
        H,
        seed!,
        nrfind,
        E::U,
        dim::T,
        params,;
        newtoniter::Int = 10,
        nrabstol::U = eps(U),
    ) where {U<:Real, T<:Int}

    ξ = zeros(dim)

    dξ = zeros(dim)

    seed!(ξ, params, 0.0)

    idx, Δξ, ξ_max, ξ_min = nrfind(ξ, params, 0.0)
    
    ξaux = [Taylor1(0.0, 1) for i in 1:dim]

    ξaux[idx][1] = 1.0

    Haux = Taylor1(0.0, 1) 

    ξ_old = Δξ 

    ξ_new = ξ_old + Δξ

    ξ[idx] = ξ_old
        
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
                    @warn "Energy tolerance was exceded : $(abs( Haux[0] - E ))"
                end

                ξ[idx] = ξaux[idx][0]

                f!(dξ, ξ, params, 0.0)

                if dξ[idx] != 0

                    return ξ

                end

            end
    
            ξ_old = ξ_new
    
            E_old = E_new
    
            ξ_new = ξ_old + Δξ
            
        end

    end

    @warn(" No se encontraron soluciones a esa energía. ")

end


function EnergyFinding!(
        f!,
        H,
        ξ::Array{U, 1},
        nrfind,
        E::U,
        params;
        newtoniter::Int = 10,
        nrabstol::U = eps(U),
    ) where {U<:Real}

    dim = length(ξ)

    dξ = zeros(dim)

    idx, Δξ, ξ_max, ξ_min = nrfind(ξ, params, 0.0)
    
    ξaux = [Taylor1(0.0, 1) for i in 1:dim]

    ξaux[idx][1] = 1.0

    Haux = Taylor1(0.0, 1) 

    ξ_old = ξ_min + Δξ 

    ξ_new = ξ_old + Δξ

    ξ[idx] = ξ_old
        
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
                    @warn "Energy tolerance was exceded : $(abs( Haux[0] - E ))"
                end

                ξ[idx] = ξaux[idx][0]

                f!(dξ, ξ, params, 0.0)

                if dξ[idx] != 0

                    return nothing

                end

            end
    
            ξ_old = ξ_new
    
            E_old = E_new
    
            ξ_new = ξ_old + Δξ
            
        end

    end

    @warn(" No se encontraron soluciones a esa energía. ")

end