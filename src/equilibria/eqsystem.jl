function EquilibriaJacobian!(J::Array{U, 2}, dx::Array{TaylorN{U}, 1}, Φ::Array{U, 1}, n::T) where {U<:Real, T<:Integer}

    TaylorSeries.jacobian!(J, dx)
    @inbounds J[n, :] .= Φ

end

function EquilibriaSystem!(F::Array{U, 1}, dx::Array{TaylorN{U}, 1}, xzero::Array{U, 1},
                           x1::Array{U, 1}, x0::Array{U, 1}, Φ::Array{U, 1},
                           Δs::U, n::T) where {U<:Real, T<:Integer}

    TaylorSeries.evaluate!(dx, xzero, F)
    @inbounds F[n] = dot(x1 - x0, Φ) - Δs

end