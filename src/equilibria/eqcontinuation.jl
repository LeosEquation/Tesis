function EquilibriaContinuation(f!, x_ini::Array{U, 1}, Φ_ini::Array{U, 1}, params, 
                                Δs::U, maxsteps::T, tol::U, ite::T) where {U<:Real, T<:Integer}

    ###

    n = length(x_ini)
    xzero = zero(x_ini)

    ###

    x = Array{U,2}(undef, maxsteps, n)
    realev = Array{U,2}(undef, maxsteps, n-1)
    imagev = Array{U,2}(undef, maxsteps, n-1)
    lp = Array{U, 1}[]
    lpeval= Array{Complex{U}, 1}[]
    lpevec= Array{Complex{U}, 2}[]
    bp = Array{U, 1}[]
    bpdir  = Array{U, 1}[]
    bpeval= Array{Complex{U}, 1}[]
    bpevec= Array{Complex{U}, 2}[]
    h = Array{U, 1}[]
    heval= Array{Complex{U}, 1}[]
    hevec= Array{Complex{U}, 2}[]
    stbl = Array{Bool,1}(undef, maxsteps)

    ###

    variable_names = [string("δx", TaylorSeries.subscriptify(i)) for i in 1:n]

    TaylorSeries.set_variables(U, variable_names, order = 2)

    δx = TaylorN.(1:n, order = 1)
    δy = TaylorN.(1:n, order = 2)

    xaux = copy(δx)
    dx = zero(δx)
    yaux = copy(δy)
    dy = zero(δy)

    ###

    x0 = copy(x_ini)
    lptest0 = 0.0
    htest0  = 0.0
    bptest0 = 0.0

    ###

    x1 = copy(x_ini)
    lptest1 = 0.0
    htest1  = 0.0
    bptest1 = 0.0

    ###

    m = ((n - 1) * (n - 2)) ÷ 2
    C = zeros(U, m, m)
    Jeval = zeros(U, n, n)
    dxeval = zeros(U, n)

    ###
    
    λ = Array{Complex{U}, 1}(undef, n-1)
    ν = Array{U, 1}(undef, n-1)
    Dx = Array{U, 2}(undef, n-1, n-1)

    ###

    y = copy(x_ini)

    ###

    J = zeros(U, n, n)
    F = zeros(U, n)
    Φ = zeros(U, n)
    v = zeros(U, n)
    v[n] = one(U)

    ###

    JLP = zeros(U, 2*n-1, 2*n-1)
    FLP = zeros(U, 2*n-1)
    qLP = zeros(U, 2*n-1)

    ###

    JBP = zeros(U, 2*n, 2*n)
    FBP = zeros(U, 2*n)
    qBP = zeros(U, 2*n)
    ΦBP = zeros(U, n)

    ###

    JH = zeros(U, 3*n-1, 3*n-1)
    FH = zeros(U, 3*n-1)
    qH = zeros(U, 3*n-1)
    w1 = zeros(U, n-1)
    w2 = zeros(U, n-1)

    ###

    x[1, :] .= x0

    for i in 1:n
        xaux[i][0][1] = x1[i]
    end

    f!(dx, xaux, params, 0.0)

    TaylorSeries.jacobian!(J, dx)

    NS = nullspace(J)

    if size(NS, 2) != 2
        @show(NS)
        error("This point is not a branch point. Exiting.")
    end

    a, b = NS \ Φ_ini

    Φ .= - b * NS[:, 1] + a * NS[:, 2]

    EquilibriaJacobian!(J, dx, Φ, n)
    EquilibriaSystem!(F, dx, xzero, x1, x0, Φ, Δs, n)

    for j in 1:n-1
        for i in 1:n-1
            Dx[i,j] = J[i, j]
        end
    end

    λ .= eigvals(Dx)

    BiProduct!(2.0, Dx, LinearAlgebra.I, C, n-1)

    lptest0 = Φ[n]
    bptest0 = det(J)
    htest0  = det(C)

    for j in 1:n-1
        realev[1, j] = real(λ[j])
        imagev[1, j] = imag(λ[j])
    end

    if all(realev[1, :] .<= 0.0)
        stbl[1] = true
    else
        stbl[1] = false
    end

    #

    i = 2

    while i <= maxsteps

        for k in 1:n
            x1[k] = x0[k] + Δs * Φ[k]
            xaux[k][0][1] = x1[k]
        end

        f!(dx, xaux, params, 0.0)
        EquilibriaJacobian!(J, dx, Φ, n)
        EquilibriaSystem!(F, dx, xzero, x1, x0, Φ, Δs, n)

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
            EquilibriaJacobian!(J, dx, Φ, n)
            EquilibriaSystem!(F, dx, xzero, x1, x0, Φ, Δs, n)    

            j += 1

        end

        if norm(F) >= tol
            @warn("La norma del sistema no convergió. Abortando.")
            break
        end

        x[i, :] .= x1

        Φ .= J \ v

        normalize!(Φ)

        for j in 1:n-1
            for k in 1:n-1
                Dx[k,j] = J[k, j]
            end
        end

        λ .= eigvals(Dx)

        BiProduct!(2.0, Dx, LinearAlgebra.I, C, n-1)

        lptest1 = Φ[n]
        bptest1 = det(J)
        htest1  = det(C)

        for j in 1:n-1
            realev[i, j] = real(λ[j])
            imagev[i, j] = imag(λ[j])
        end

        if all(realev[i, :] .<= 0.0)
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

            TaylorSeries.jacobian!(Jeval, dy)

            for j in 1:n-1
                for k in 1:n-1
                    Dx[k,j] = Jeval[k, j]
                end
            end

            Dxeigen = eigen(Dx)

            idx = argmin(abs.(real.(Dxeigen.values)))

            for j in 1:n-1
                ν[j]  = real(Dxeigen.vectors[j, idx])
            end

            LPFinding!(f!, y, params, ν,
                       qLP, JLP, FLP,
                       dxeval, Jeval, xzero,
                       dy, yaux, ite, tol, n)

            for j in 1:n-1
                for k in 1:n-1
                    Dx[k,j] = Jeval[k, j]
                end
            end

            Dxeigen = eigen(Dx)

            push!(lp, copy(y))
            push!(lpeval, Dxeigen.values)
            push!(lpevec, Dxeigen.vectors)

        end

        if sign(htest0) != sign(htest1) && count(!iszero, imagev[i, :]) >= 2

            for j in 1:n
                y[j] = 0.5 * (x1[j] + x0[j])
                yaux[j][0][1] = y[j]
            end

            f!(dy, yaux, params, 0.0)

            TaylorSeries.jacobian!(Jeval, dy)

            for j in 1:n-1
                for k in 1:n-1
                    Dx[k,j] = Jeval[k, j]
                end
            end

            Dxeigen = eigen(Dx)

            indices = findall(x -> imag(x) != 0, Dxeigen.values)

            closest_index = argmin(abs.(real.(Dxeigen.values[indices])))

            idx = indices[closest_index]

            HopfFinding!(f!, y, params, Dxeigen.values[idx], Dxeigen.vectors[:, idx], 
                         w1, w2,
                         qH, JH, FH, 
                         dxeval, Jeval, xzero, 
                         dy, yaux, ite, tol, n)

            for j in 1:n-1
                for k in 1:n-1
                    Dx[k,j] = Jeval[k, j]
                end
            end

            Dxeigen = eigen(Dx)

            push!(h, copy(y))
            push!(heval, Dxeigen.values)
            push!(hevec, Dxeigen.vectors)

        end

        if sign(bptest0) != sign(bptest1)

            for j in 1:n
                y[j] = 0.5 * (x1[j] + x0[j])
                yaux[j][0][1] = y[j]
            end

            f!(dy, yaux, params, 0.0)

            TaylorSeries.jacobian!(Jeval, dy)

            for j in 1:n-1
                for k in 1:n-1
                    Dx[k,j] = Jeval[k, j]
                end
            end

            Dxeigen = eigen(Dx)
            idx = argmin(abs.(real.(Dxeigen.values)))
            for j in 1:n-1
                ν[j]  = real(Dxeigen.vectors[j, idx])
            end

            BPFinding!(f!, y, params, ν, 
                       qBP, JBP, FBP, 
                       dxeval, Jeval, xzero, 
                       dy, yaux, δy, ite, tol, n)

            for j in 1:n-1
                for k in 1:n-1
                    Dx[k,j] = Jeval[k, j]
                end
            end

            Jeval[n, :] .= Φ

            ΦBP .= Jeval \ v

            Dxeigen = eigen(Dx)

            push!(bp, copy(y))
            push!(bpdir, copy(ΦBP))
            push!(bpeval, Dxeigen.values)
            push!(bpevec, Dxeigen.vectors)

        end

        x0 .= x1
        htest0 = htest1
        bptest0 = bptest1
        lptest0 = lptest1

        i += 1

    end

    return EquilibriumBranch(x, realev, imagev, 
                             lp, lpeval, lpevec,
                             bp, bpdir, bpeval, bpevec, 
                             h, heval, hevec, 
                             stbl, n, i-1)

end

function EquilibriaContinuation(f!, x_ini::Array{U, 1}, params, 
                                Δs::U, maxsteps::T, tol::U, ite::T) where {U<:Real, T<:Integer}

    ###

    n = length(x_ini)
    xzero = zero(x_ini)

    ###

    x = Array{U,2}(undef, maxsteps, n)
    realev = Array{U,2}(undef, maxsteps, n-1)
    imagev = Array{U,2}(undef, maxsteps, n-1)
    lp = Array{U, 1}[]
    lpeval= Array{Complex{U}, 1}[]
    lpevec= Array{Complex{U}, 2}[]
    bp = Array{U, 1}[]
    bpdir  = Array{U, 1}[]
    bpeval= Array{Complex{U}, 1}[]
    bpevec= Array{Complex{U}, 2}[]
    h = Array{U, 1}[]
    heval= Array{Complex{U}, 1}[]
    hevec= Array{Complex{U}, 2}[]
    stbl = Array{Bool,1}(undef, maxsteps)

    ###

    variable_names = [string("δx", TaylorSeries.subscriptify(i)) for i in 1:n]

    TaylorSeries.set_variables(U, variable_names, order = 2)

    δx = TaylorN.(1:n, order = 1)
    δy = TaylorN.(1:n, order = 2)

    xaux = copy(δx)
    dx = zero(δx)
    yaux = zero(δy)
    dy = copy(δy)

    ###

    x0 = copy(x_ini)
    lptest0 = 0.0
    htest0  = 0.0
    bptest0 = 0.0

    ###

    x1 = copy(x_ini)
    lptest1 = 0.0
    htest1  = 0.0
    bptest1 = 0.0

    ###

    m = ((n - 1) * (n - 2)) ÷ 2
    C = zeros(U, m, m)
    Jeval = zeros(U, n, n)
    dxeval = zeros(U, n)

    ###
    
    λ = Array{Complex{U}, 1}(undef, n-1)
    ν = Array{U, 1}(undef, n-1)
    Dx = Array{U, 2}(undef, n-1, n-1)

    ###

    y = copy(x_ini)

    ###

    J = zeros(U, n, n)
    F = zeros(U, n)
    Φ = zeros(U, n)
    v = zeros(U, n)
    v[n] = one(U)

    ###

    JLP = zeros(U, 2*n-1, 2*n-1)
    FLP = zeros(U, 2*n-1)
    qLP = zeros(U, 2*n-1)

    ###

    JBP = zeros(U, 2*n, 2*n)
    FBP = zeros(U, 2*n)
    qBP = zeros(U, 2*n)
    ΦBP = zeros(U, n)

    ###

    JH = zeros(U, 3*n-1, 3*n-1)
    FH = zeros(U, 3*n-1)
    qH = zeros(U, 3*n-1)
    w1 = zeros(U, n-1)
    w2 = zeros(U, n-1)

    ###
    x[1, :] .= x0

    for i in 1:n
        xaux[i][0][1] = x1[i]
    end

    f!(dx, xaux, params, 0.0)

    TaylorSeries.jacobian!(J, dx)

    NS = nullspace(J)

    if size(NS, 2) == 0
        error("This point has not a direction branch. Exiting.")
    end

    if size(NS, 2)  > 1
        error("This point is a codimention 2. Exiting.")
    end

    Φ .= NS

    EquilibriaJacobian!(J, dx, Φ, n)
    EquilibriaSystem!(F, dx, xzero, x1, x0, Φ, Δs, n)

    for j in 1:n-1
        for i in 1:n-1
            Dx[i,j] = J[i, j]
        end
    end

    λ .= eigvals(Dx)

    BiProduct!(2.0, Dx, LinearAlgebra.I, C, n-1)

    lptest0 = Φ[n]
    bptest0 = det(J)
    htest0  = det(C)

    for j in 1:n-1
        realev[1, j] = real(λ[j])
        imagev[1, j] = imag(λ[j])
    end

    if all(realev[1, :] .<= 0.0)
        stbl[1] = true
    else
        stbl[1] = false
    end

    #

    i = 2

    while i <= maxsteps

        for k in 1:n
            x1[k] = x0[k] + Δs * Φ[k]
            xaux[k][0][1] = x1[k]
        end

        f!(dx, xaux, params, 0.0)
        EquilibriaJacobian!(J, dx, Φ, n)
        EquilibriaSystem!(F, dx, xzero, x1, x0, Φ, Δs, n)

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
            EquilibriaJacobian!(J, dx, Φ, n)
            EquilibriaSystem!(F, dx, xzero, x1, x0, Φ, Δs, n)    

            j += 1

        end

        if norm(F) >= tol
            @warn("La norma del sistema no convergió. Abortando.")
            break
        end

        x[i, :] .= x1

        Φ .= J \ v

        normalize!(Φ)

        for j in 1:n-1
            for k in 1:n-1
                Dx[k,j] = J[k, j]
            end
        end

        λ .= eigvals(Dx)

        BiProduct!(2.0, Dx, LinearAlgebra.I, C, n-1)

        lptest1 = Φ[n]
        bptest1 = det(J)
        htest1  = det(C)

        for j in 1:n-1
            realev[i, j] = real(λ[j])
            imagev[i, j] = imag(λ[j])
        end

        if all(realev[i, :] .<= 0.0)
            stbl[i] = true
        else
            stbl[i] = false
        end

        if sign(lptest0) != sign(lptest1) && sign(bptest0) == sign(bptest1)

            for j in 1:n
                y[j] = 0.5 * (x1[j] + x0[j])
                # yaux[j] = δy[j]
                yaux[j][0][1] = y[j]
            end

            f!(dy, yaux, params, 0.0)

            TaylorSeries.jacobian!(Jeval, dy)

            for j in 1:n-1
                for k in 1:n-1
                    Dx[k,j] = Jeval[k, j]
                end
            end

            Dxeigen = eigen(Dx)

            idx = argmin(abs.(real.(Dxeigen.values)))

            for j in 1:n-1
                ν[j]  = real(Dxeigen.vectors[j, idx])
            end

            LPFinding!(f!, y, params, ν,
                       qLP, JLP, FLP,
                       dxeval, Jeval, xzero,
                       dy, yaux, ite, tol, n)

            for j in 1:n-1
                for k in 1:n-1
                    Dx[k,j] = Jeval[k, j]
                end
            end

            Dxeigen = eigen(Dx)

            push!(lp, copy(y))
            push!(lpeval, Dxeigen.values)
            push!(lpevec, Dxeigen.vectors)

        end

        if sign(htest0) != sign(htest1) && count(!iszero, imagev[i, :]) >= 2

            for j in 1:n
                y[j] = 0.5 * (x1[j] + x0[j])
                yaux[j][0][1] = y[j]
            end

            f!(dy, yaux, params, 0.0)

            TaylorSeries.jacobian!(Jeval, dy)

            for j in 1:n-1
                for k in 1:n-1
                    Dx[k,j] = Jeval[k, j]
                end
            end

            Dxeigen = eigen(Dx)

            indices = findall(x -> imag(x) != 0, Dxeigen.values)

            closest_index = argmin(abs.(real.(Dxeigen.values[indices])))

            idx = indices[closest_index]

            HopfFinding!(f!, y, params, Dxeigen.values[idx], Dxeigen.vectors[:, idx], 
                         w1, w2,
                         qH, JH, FH, 
                         dxeval, Jeval, xzero, 
                         dy, yaux, ite, tol, n)

            for j in 1:n-1
                for k in 1:n-1
                    Dx[k,j] = Jeval[k, j]
                end
            end

            Dxeigen = eigen(Dx)

            push!(h, copy(y))
            push!(heval, Dxegien.values)
            push!(hevec, Dxeigen.vectors)

        end

        if sign(bptest0) != sign(bptest1)

            for j in 1:n
                y[j] = 0.5 * (x1[j] + x0[j])
                yaux[j][0][1] = y[j]
            end

            f!(dy, yaux, params, 0.0)

            TaylorSeries.jacobian!(Jeval, dy)

            for j in 1:n-1
                for k in 1:n-1
                    Dx[k,j] = Jeval[k, j]
                end
            end

            Dxeigen = eigen(Dx)
            idx = argmin(abs.(real.(Dxeigen.values)))
            for j in 1:n-1
                ν[j]  = real(Dxeigen.vectors[j, idx])
            end

            BPFinding!(f!, y, params, ν, 
                       qBP, JBP, FBP, 
                       dxeval, Jeval, xzero, 
                       dy, yaux, δy, ite, tol, n)

            for j in 1:n-1
                for k in 1:n-1
                    Dx[k,j] = Jeval[k, j]
                end
            end

            Jeval[n, :] .= Φ

            ΦBP .= Jeval \ v

            Dxeigen = eigen(Dx)

            push!(bp, copy(y))
            push!(bpdir, copy(ΦBP))
            push!(bpeval, Dxeigen.values)
            push!(bpevec, Dxeigen.vectors)

        end

        x0 .= x1
        htest0 = htest1
        bptest0 = bptest1
        lptest0 = lptest1

        i += 1

    end

    return EquilibriumBranch(x, realev, imagev, 
                             lp, lpeval, lpevec,
                             bp, bpdir, bpeval, bpevec, 
                             h, heval, hevec, 
                             stbl, n, i-1)

end

