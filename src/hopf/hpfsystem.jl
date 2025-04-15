function HopfJacobian!(J::Array{U, 2}, Jeval::Array{U, 2}, q::Array{U, 1}, 
                       dx::Array{TaylorN{U}, 1}, w1::Array{U, 1}, w2::Array{U, 1}, n::T) where {U<:Real, T<:Integer}
    
    i1 = n-1
    # j1 = n 

    i2 = 2*n-2
    j2 = 2*n-1

    i3 = 3*n-2
    i4 = 3*n-1

    for i in 1:n-1

        for j in 1:n

            J[i, j] = Jeval[i, j]

            J[i1+i, j] = sum(evaluate(differentiate(differentiate(dx[i], k), j)) * q[n+k] for k in 1:n-1)

            J[i2+i, j] = sum(evaluate(differentiate(differentiate(dx[i], k), j)) * q[j2+k] for k in 1:n-1)

            if j < n

                J[i1+i, n+j] = Jeval[i, j]

                J[i2+i, j2+j] = Jeval[i, j]

            end
        end

        J[i1+i, j2+i] = q[i4]
        J[i1+i, i4] = q[j2+i]

        J[i2+i, n+i] = - q[i4]
        J[i2+i, i4] = - q[n+i]

        J[i3, n+i] = w1[i]
        J[i3, j2+i] = w2[i]

        J[i4, n+i] = - w2[i]
        J[i4, j2+i] = w1[i]
    end

end

function HopfSystem!(F::Array{U,1}, dxeval::Array{U, 1}, 
                     q::Array{U, 1}, Jeval::Array{U, 2}, 
                     w1::Array{U, 1}, w2::Array{U, 1}, n::T) where {U<:Real, T<:Integer}

    i1 = n-1
    i2 = 2*n-2
    j2 = 2*n-1

    for i in 1:n-1
        F[i] = dxeval[i]
        F[i1+i] = sum(Jeval[i, k] * q[n+k] for k in 1:n-1) + q[3*n-1] * q[j2+i]
        F[i2+i] = sum(Jeval[i, k] * q[j2+k] for k in 1:n-1) - q[3*n-1] * q[n+i]
    end

    F[3*n-2] = sum(w1[i] * q[n+i] + w2[i] * q[j2+i] for i in 1:n-1) - 1.0
    F[3*n-1] = sum(w1[i] * q[j2+i] - w2[i] * q[n+i] for i in 1:n-1)

end


















































function HopfJacobian!(J::Matrix{Float64}, Fx::Matrix{TaylorN{Float64}}, w::Vector{Float64}, 
                       dx::Vector{TaylorN{Float64}}, x1::Vector{Float64}, Φ::Vector{Float64}, n::Int64)

    J[1:n, 1:n] .= Jeval

    for j in 1:n
    J[n-1:2*n-4, j] .= evaluate(differentiate.(Fx ^ 2, j)) * x1[n+1:2*n-2]
    end

    J[n-1:2*n-4, n+1:2*n-2] .= evaluate(Fx) ^ 2

    for j in 1:n-2
    J[n-2+j, n+j] += x1[2*n-1]
    end

    J[n-1:2*n-4, 2*n-1] .= x1[n+1:2*n-2]

    J[2*n-3, n+1:2*n-2] .= 2.0 * x1[n+1:2*n-2]

    J[2*n-2, n+1:2*n-2] .= w

    J[2*n-1, :] .= Φ

end

function HopfSystem!(F::Vector{Float64}, Fx::Matrix{TaylorN{Float64}}, dx::Vector{TaylorN{Float64}},
x0::Vector{Float64}, x1::Vector{Float64}, w::Vector{Float64},
Φ::Vector{Float64}, Δs::Float64, n::Int64)

F[1:n-2] .= evaluate(dx[1:n-2])
F[(n-1):(2*n-4)] .= evaluate(Fx) ^ 2 * x1[n+1:2*n-2] + x1[2*n-1] * x1[n+1:2*n-2]
F[2*n-3] = dot(x1[n+1:2*n-2], x1[n+1:2*n-2]) - 1.0
F[2*n-2] = dot(w, x1[n+1:2*n-2])
F[2*n-1] = dot(x1 - x0, Φ) - Δs

end