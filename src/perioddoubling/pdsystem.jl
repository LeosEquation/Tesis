
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

function LPJacobian!(J::Matrix{Float64}, Fx::Matrix{TaylorN{Float64}}, dx::Vector{TaylorN{Float64}},
                     x1::Vector{Float64}, n::Int64)

    @inbounds J[1:n+1, 1:n+1] .= TaylorSeries.jacobian(dx)

    @inbounds for j in 1:n+1
        J[n+1:2*n, j] .= evaluate(differentiate.(Fx, j)) * x1[n+2:2*n+1]
    end

    @inbounds J[n+1:2*n, n+2:2*n+1] .= evaluate(Fx)

    @inbounds J[2*n+1, n+2:2*n+1] .= 2*x1[n+2:2*n+1]

end

function LPSystem!(F::Vector{Float64}, Fx::Matrix{TaylorN{Float64}}, dx::Vector{TaylorN{Float64}},
                   x1::Vector{Float64}, n::Int64)

    @inbounds F[1:n+1] .= evaluate(dx)
    @inbounds F[n+1:2*n] .= evaluate(Fx) * x1[n+2:2*n+1]
    @inbounds F[2*n+1] = dot(x1[n+2:2*n+1], x1[n+2:2*n+1]) - 1.0

end