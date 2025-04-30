function LPFinding!(f!, x::Array{U, 1}, params, v::Array{U, 1},
                    q::Array{U, 1}, J::Array{U, 2}, F::Array{U, 1}, 
                    dxeval::Array{U, 1}, Jeval::Array{U, 2},
                    dy::Array{TaylorN{U}, 1}, yaux::Array{TaylorN{U},1}, 
                    ite::T, tol::U, ordtup::Array{NTuple{N, T}, 2}, n::T) where {U<:Real, T<:Integer, N}

        q[1:n] .= x
        q[n+1:2*n-1] .= v

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

        LPSystem!(F, Jeval, dxeval, q, n)
        LPJacobian!(J, Jeval, dy, q, ordtup, n)

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
        
                LPSystem!(F, Jeval, dxeval, q, n)
                LPJacobian!(J, Jeval, dy, q, ordtup, n)

                k += 1

        end

        if norm(F) > tol
                @warn("The limit point tolerance was exceded.")
        end 

        for i in 1:n
                x[i] = q[i]
        end

        for i in 1:n-1
                v[i] = q[n+i]
        end

end