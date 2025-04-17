
function LPJacobian!(J::Matrix{Float64}, Fx::Matrix{TaylorN{Float64}}, dx::Vector{TaylorN{Float64}},
                     x1::Vector{Float64}, Φ::Vector{Float64}, n::Int64)

    @inbounds J[1:n+2, 1:n+2] .= TaylorSeries.jacobian(dx)

    @inbounds for j in 1:n+2
        J[n+1:2*n, j] .= evaluate(differentiate.(Fx, j)) * x1[n+3:2*n+2]
    end

    @inbounds J[n+1:2*n, n+3:2*n+2] .= evaluate(Fx)

    @inbounds J[2*n+1, n+3:2*n+2] .= 2*x1[n+3:2*n+2]

    @inbounds J[2*n+2, :] .= Φ

end

function LPSystem!(F::Vector{Float64}, Fx::Matrix{TaylorN{Float64}}, dx::Vector{TaylorN{Float64}},
                   x0::Vector{Float64}, x1::Vector{Float64},
                   Φ::Vector{Float64}, Δs::Float64, n::Int64)

    @inbounds F[1:n+2] .= evaluate(dx)
    @inbounds F[(n+1):(2*n)] .= evaluate(Fx) * x1[n+3:2*n+2]
    @inbounds F[2*n+1] = dot(x1[n+3:2*n+2], x1[n+3:2*n+2]) - 1.0
    @inbounds F[2*n+2] = dot(x1 - x0, Φ) - Δs

end

function LPJacobian!(J::Array{U, 2}, Jeval::Array{U, 2}, dx::Array{TaylorN{U}, 1},
                     q::Array{U, 1}, n::T) where {U<:Real, T<:Integer}

    i1 = n-1
    i2 = 2*n - 1

    for i in 1:n-1
        for j in 1:n
            J[i, j] = Jeval[i, j]
            J[i1+i, j] = sum(evaluate(differentiate(differentiate(dx[i], k), j)) * q[n+k] for k in 1:n-1)
            if j < n
                J[i1+i, n+j] = Jeval[i,j]
            end
        end
        J[i2, n+i] = 2.0 * q[n+i]
    end

end

function LPSystem!(F::Array{U, 1}, Jeval::Array{U, 2}, dxeval::Array{U, 1},
                   q::Array{U, 1}, n::T) where {U<:Real, T<:Integer}

    i1 = n-1
    for i in 1:n-1
        F[i] = dxeval[i]
        F[i1+i] = sum(Jeval[i, k] * q[n+k] for k in 1:n-1)
    end
    F[2*n-1] = sum(q[n+i]^2 for i in 1:n-1) - 1.0

end