function BPFinding!(f!, x::Array{U, 1}, params, v::Array{U, 1}, 
                    q::Array{U, 1}, J::Array{U, 2}, F::Array{U, 1}, 
                    dxeval::Array{U, 1}, Jeval::Array{U, 2}, xzero::Array{U, 1},
                    dy::Array{TaylorN{U}, 1}, yaux::Array{TaylorN{U},1}, 
                    ite::T, tol::U, n::T) where {U<:Real, T<:Integer}

    @inbounds q[1:n] .= x
    @inbounds q[n+1:2*n-1] .= v
    @inbounds q[2*n] = zero(U)

    for i in 1:n
        yaux[i][0][1] = q[i]
    end

    f!(dy, yaux, params, zero(U))

    evaluate!(dy, xzero, dxeval)
    TaylorSeries.jacobian!(Jeval, dy)

    BPSystem!(F, dxeval, q, Jeval, n)
    BPJacobian!(J, Jeval, q, dy, n)

    k = 1

    while k <= ite && norm(F) > tol

        if det(J) == 0.0
            break
        end

        q .-= J \ F 

        for i in 1:n
            yaux[i][0][1] = q[i]
        end

        f!(dy, yaux, params, zero(U))

        evaluate!(dy, xzero, dxeval)
        TaylorSeries.jacobian!(Jeval, dy)

        BPSystem!(F, dxeval, q, Jeval, n)
        BPJacobian!(J, Jeval, q, dy, n)

        k += 1

    end

    if norm(F) > tol
        @warn("The branch point tolerance was exceded.")
    end 

    for i in 1:n
        x[i] = q[i]
    end
    for i in 1:n-1
            v[i] = q[n+i]
    end

end