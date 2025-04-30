function EquilibriaJacobian!(J::Array{U, 2}, dx::Array{TaylorN{U}, 1}, 
                             Φ::Array{U, 1}, n::T) where {U<:Real, T<:Integer}

    for j in 1:n
        for i in 1:n-1
            @inbounds J[i,j] = dx[i][1][j]
        end
    end
    @inbounds J[n, :] .= Φ

end

function EquilibriaSystem!(F::Array{U, 1}, dx::Array{TaylorN{U}, 1},
                           x1::Array{U, 1}, x0::Array{U, 1}, Φ::Array{U, 1},
                           Δs::U, n::T) where {U<:Real, T<:Integer}

    for i in 1:n-1
        @inbounds F[i] = dx[i][0][1]
    end
    @inbounds F[n] = sum((x1[i] - x0[i])*Φ[i] for i in 1:n) - Δs

end