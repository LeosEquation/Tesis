function HopfFinding!(f!, x::Array{U, 1}, params,
                      realv::Array{U, 1}, imagv::Array{U, 1}, ω::U,
                      w1::Array{U, 1}, w2::Array{U, 1}, 
                      q::Array{U, 1}, J::Array{U, 2}, F::Array{U, 1}, 
                      dxeval::Array{U, 1}, Jeval::Array{U, 2}, xzero::Array{U, 1},
                      dy::Array{TaylorN{U}, 1}, yaux::Array{TaylorN{U},1}, 
                      ite::T, tol::U, n::T) where {U<:Real, T<:Integer}

    q[1:n] .= x
    q[n+1:2*n-1] .= realv
    q[2*n:3*n-2] .= imagv
    q[3*n-1] = ω

    w1 .= realv
    w2 .= imagv

    for i in 1:n
        yaux[i][0][1] = q[i]
    end

    f!(dy, yaux, params, zero(U))

    evaluate!(dy, xzero, dxeval)
    TaylorSeries.jacobian!(Jeval, dy)

    HopfSystem!(F, dxeval, q, Jeval, w1, w2, n)
    HopfJacobian!(J, Jeval, q, dy, w1, w2, n)
 
    k = 1

    # println(" ite = $k, ||F|| = $(norm(F))")

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

        HopfSystem!(F, dxeval, q, Jeval, w1, w2, n)
        HopfJacobian!(J, Jeval, q, dy, w1, w2, n)

        k += 1

        # println(" ite = $k, ||F|| = $(norm(F))")

    end

    if norm(F) > tol
        @warn("The hopf point tolerance was exceded.")
    end 

    for i in 1:n
        x[i] = q[i]
    end

    for i in 1:n-1
        realv[i] = q[n+i]
        imagv[i] = q[2*n-1+i]
    end

    ω = q[3*n-1]

end