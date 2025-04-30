
function EquilibriaContinuation(f!, bp_ini::BranchPoint{U}, params, pmin::U, pmax::U, 
                                Δs::U, maxsteps::T, tol::U, ite::T) where {U<:Real, T<:Integer}

    ###
    n = length(bp_ini.x)

    ordtup = [ntuple(k -> count(==(k), (i, j)), TaylorSeries.get_numvars()) for i in 1:n, j in 1:n]

    ###
    x = Array{U,2}(undef, maxsteps, n)
    λ = Array{Complex{U},2}(undef, maxsteps, n-1)
    lp = LimitPoint{U}[]
    bp = BranchPoint{U}[]
    hp = HopfPoint{U}[]
    stbl = Array{Bool,1}(undef, maxsteps)

    ###
    δx = TaylorN.(1:n, order = 1)
    δy = TaylorN.(1:n, order = 2)
    xaux = copy(δx)
    dx = zero(δx)
    yaux = copy(δy)
    dy = zero(δy)

    ###
    x0 = copy(bp_ini.x)
    lptest0 = 0.0
    htest0  = 0.0
    bptest0 = 0.0
    Φ0 = zero(bp_ini.dir)

    ###
    x1 = copy(bp_ini.x)
    lptest1 = 0.0
    htest1  = 0.0
    bptest1 = 0.0
    Φ1 = zero(bp_ini.dir)

    ###
    m = ((n - 1) * (n - 2)) ÷ 2
    C = zeros(U, m, m)
    Jeval = zeros(U, n-1, n)
    dxeval = zeros(U, n-1)

    ###
    λ_ = Array{Complex{U}, 1}(undef, n-1)

    ###
    y = copy(bp_ini.x)

    ###
    J = zeros(U, n, n)
    F = zeros(U, n)
    v = zeros(U, n)
    v[n] = one(U)

    ###
    JLP = zeros(U, 2*n-1, 2*n-1)
    FLP = zeros(U, 2*n-1)
    qLP = zeros(U, 2*n-1)
    vLP = zeros(U, n-1)

    ###
    JBP = zeros(U, 2*n, 2*n)
    FBP = zeros(U, 2*n)
    qBP = zeros(U, 2*n)
    ΦBP = zeros(U, n)
    wBP = zeros(U, n-1)
    h1 = zero(U)
    h0 = zero(U)

    ###
    JH = zeros(U, 3*n-1, 3*n-1)
    FH = zeros(U, 3*n-1)
    qH = zeros(U, 3*n-1)
    realvH = zeros(U, n-1)
    imagvH = zeros(U, n-1)
    realwH = zeros(U, n-1)
    imagwH = zeros(U, n-1)
    ωH = zero(U)

    ###
    x[1, :] .= x0

    for i in 1:n
        xaux[i][0][1] = x1[i]
    end

    f!(dx, xaux, params, 0.0)
    
    for j in 1:n
        for i in 1:n-1
            Jeval[i,j] = dx[i][1][j]
        end
    end

    NS = nullspace(Jeval)

    if size(NS, 2) != 2
        error("Initial point does not a valid branch point. Cannot continue.")
    end

    a, b = NS \ bp_ini.dir

    Φ0 .= - b * NS[:, 1] + a * NS[:, 2]

    EquilibriaJacobian!(J, dx, Φ0, n)
    EquilibriaSystem!(F, dx, x1, x0, Φ0, Δs, n)

    for j in 1:n
        for k in 1:n-1
            Jeval[k, j] = dx[k][1][j]
        end
    end

    λ_ .= eigvals(Jeval[:, 1:n-1])
    BiProduct!(2.0, Jeval, LinearAlgebra.I, C, n-1)

    lptest0 = Φ0[n]
    bptest0 = det(J)
    htest0  = det(C)

    for j in 1:n-1
        λ[1, j] = λ_[j]
    end

    if all(real(λ_) .<= 0.0)
        stbl[1] = true
    else
        stbl[1] = false
    end

    #
    i = 2

    while i <= maxsteps && (pmin <= x1[n] <= pmax)

        for k in 1:n
            x1[k] = x0[k] + Δs * Φ0[k]
            xaux[k][0][1] = x1[k]
        end

        f!(dx, xaux, params, 0.0)
        EquilibriaJacobian!(J, dx, Φ0, n)
        EquilibriaSystem!(F, dx, x1, x0, Φ0, Δs, n)

        j = 1

        while j <= ite && norm(F) > tol

            if abs(det(J)) == 0.0
                break
            end

            x1 .-= J \ F
    
            for k in 1:n
                xaux[k][0][1] = x1[k]
            end

            f!(dx, xaux, params, 0.0)
            EquilibriaJacobian!(J, dx, Φ0, n)
            EquilibriaSystem!(F, dx, x1, x0, Φ0, Δs, n)    

            j += 1

        end

        if j > ite
            @warn("System norm did not converge to tolerance. Stopping continuation.")
            break
        end

        x[i, :] .= x1

        Φ1 .= J \ v

        normalize!(Φ1)

        for j in 1:n
            for k in 1:n-1
                Jeval[k, j] = dx[k][1][j]
            end
        end

        λ_ .= eigvals(Jeval[:, 1:n-1])
        BiProduct!(2.0, Jeval, LinearAlgebra.I, C, n-1)

        lptest1 = Φ1[n]
        bptest1 = det(J)
        htest1  = det(C)

        for j in 1:n-1
            λ[i, j] = λ_[j]
        end
    
        if all(real(λ_) .<= 0.0)
            stbl[i] = true
        else
            stbl[i] = false
        end
    

        if sign(lptest0) != sign(lptest1) && sign(bptest0) == sign(bptest1)

            for j in 1:n
                y[j] = 0.5 * (x1[j] + x0[j])
                yaux[j][0][1] = y[j]
            end

            f!(dy, yaux, params, 0.0)

            for j in 1:n
                for k in 1:n-1
                    Jeval[k, j] = dy[k][1][j]
                end
            end

            Dxeigen = eigen(Jeval[:, 1:n-1])
            idx = argmin(abs.(real.(Dxeigen.values)))

            for j in 1:n-1
                vLP[j]  = real(Dxeigen.vectors[j, idx])
            end

            LPFinding!(f!, y, params, vLP,
                       qLP, JLP, FLP,
                       dxeval, Jeval,
                       dy, yaux, ite, tol, ordtup, n)

            Dxeigen = eigen(Jeval[:, 1:n-1])
            idx = argmin(abs.(real.(Dxeigen.values)))

            push!(lp, LimitPoint(copy(y), Dxeigen.vectors[:, idx], Dxeigen.values[idx]))

        end

        if sign(htest0) != sign(htest1) && count(!iszero, imag(λ_)) >= 2

            for j in 1:n
                y[j] = 0.5 * (x1[j] + x0[j])
                yaux[j][0][1] = y[j]
            end

            f!(dy, yaux, params, 0.0)

            for j in 1:n
                for k in 1:n-1
                    Jeval[k, j] = dy[k][1][j]
                end
            end

            Dxeigen = eigen(Jeval[:, 1:n-1])
            indices = findall(x -> imag(x) != 0, Dxeigen.values)
            closest_index = argmin(abs.(real.(Dxeigen.values[indices])))
            idx = indices[closest_index]

            ωH = imag(Dxeigen.values[idx])

            for j in 1:n-1
                realvH[j] = real(Dxeigen.vectors[j, idx])
                imagvH[j] = imag(Dxeigen.vectors[j, idx])
            end

            HopfFinding!(f!, y, params, realvH, imagvH, ωH,
                         realwH, imagwH,
                         qH, JH, FH, 
                         dxeval, Jeval, 
                         dy, yaux, ite, tol, ordtup, n)

            Dxeigen = eigen(Jeval[:, 1:n-1])
            indices = findall(x -> imag(x) != 0, Dxeigen.values)
            closest_index = argmin(abs.(real.(Dxeigen.values[indices])))
            idx = indices[closest_index]
        
            push!(hp, HopfPoint(copy(y), Dxeigen.vectors[:, idx], Dxeigen.values[idx]))

        end

        if sign(bptest0) != sign(bptest1)

            for j in 1:n
                y[j] = 0.5 * (x1[j] + x0[j])
                yaux[j][0][1] = y[j]
            end

            f!(dy, yaux, params, 0.0)

            for j in 1:n-1
                for k in 1:n-1
                    Jeval[j, k] = dy[k][1][j]
                end
            end

            Dxeigen = eigen(Jeval[:, 1:n-1])

            idx = argmin(abs.(real.(Dxeigen.values)))

            for j in 1:n-1
                wBP[j]  = real(Dxeigen.vectors[j, idx])
            end

            BPFinding!(f!, y, params, wBP, 
                       qBP, JBP, FBP, 
                       dxeval, Jeval, 
                       dy, yaux, ite, tol, ordtup, n)

            Dxeigen = eigen(Jeval[:, 1:n-1])
            
            idx = argmin(abs.(real.(Dxeigen.values)))

            h1 = sqrt(sum((x1[j] - y[j])^2 for j in 1:n))
            h0 = sqrt(sum((x0[j] - y[j])^2 for j in 1:n))

            for j in 1:n
                ΦBP[j] = (h0 * Φ0[j] + h1 * Φ1[j]) / Δs
            end

            push!(bp, BranchPoint(copy(y), Dxeigen.vectors[:, idx], Dxeigen.values[idx], copy(ΦBP)))

        end

        x0 .= x1
        Φ0 .= Φ1
        htest0 = htest1
        bptest0 = bptest1
        lptest0 = lptest1

        i += 1

    end

    return EquilibriumBranch(x[1:i-1, :], λ[1:i-1, :], lp, hp, bp, stbl[1:i-1])

end