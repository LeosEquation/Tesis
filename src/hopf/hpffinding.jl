function HopfFinding!(f!, x::Array{U, 1}, params,
                      realv::Array{U, 1}, imagv::Array{U, 1}, ω::U,
                      w1::Array{U, 1}, w2::Array{U, 1}, 
                      q::Array{U, 1}, J::Array{U, 2}, F::Array{U, 1}, 
                      dxeval::Array{U, 1}, Jeval::Array{U, 2},
                      dy::Array{TaylorN{U}, 1}, yaux::Array{TaylorN{U},1}, 
                      ite::T, tol::U, ordtup::Array{NTuple{N, T}, 2}, n::T) where {U<:Real, T<:Integer, N}

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

    for i in 1:n-1
        for j in 1:n
                Jeval[i,j] = dy[i][1][j]
        end
        dxeval[i] = dy[i][0][1]
    end

    HopfSystem!(F, dxeval, q, Jeval, w1, w2, n)
    HopfJacobian!(J, Jeval, q, dy, w1, w2, ordtup, n)
 
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

        for i in 1:n-1
            for j in 1:n
                    Jeval[i,j] = dy[i][1][j]
            end
            dxeval[i] = dy[i][0][1]
        end

        HopfSystem!(F, dxeval, q, Jeval, w1, w2, n)
        HopfJacobian!(J, Jeval, q, dy, w1, w2, ordtup, n)

        k += 1

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