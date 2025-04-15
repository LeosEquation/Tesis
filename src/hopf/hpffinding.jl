function HopfFinding!(f!, x::Array{U, 1}, params, λ::Complex{U}, v::Array{Complex{U}, 1},
                      w1::Array{U, 1}, w2::Array{U, 1}, 
                      q::Array{U, 1}, J::Array{U, 2}, F::Array{U, 1}, 
                      dxeval::Array{U, 1}, Jeval::Array{U, 2}, xzero::Array{U, 1},
                      dy::Array{TaylorN{U}, 1}, yaux::Array{TaylorN{U},1}, 
                      ite::T, tol::U, n::T) where {U<:Real, T<:Integer}

    # @show λ
    # @show v

    q[1:n] .= x

    for i in 1:n-1

        q[n+i] = real(v[i])
        w1[i] = real(v[i])

        q[2*n-1+i] = imag(v[i])
        w2[i] = imag(v[i])

    end

    q[3*n-1] = imag(λ)

    for i in 1:n
        yaux[i][0][1] = q[i]
    end

    f!(dy, yaux, params, zero(U))

    evaluate!(dy, xzero, dxeval)
    TaylorSeries.jacobian!(Jeval, dy)

    HopfSystem!(F, dxeval, q, Jeval, w1, w2, n)
    HopfJacobian!(J, Jeval, q, dy, w1, w2, n)
 
    k = 1

    while k <= ite && norm(F) > tol

        # println(" ite = $k : ||F|| = $(norm(F))")

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

    end

    if norm(F) > tol
        @warn("The hopf point tolerance was exceded.")
    end 

    for i in 1:n
        x[i] = q[i]
    end

end