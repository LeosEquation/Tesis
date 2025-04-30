
function OrbitSystem!(F::Array{U, 1}, xT::Array{TaylorN{U}, 1}, 
                      x1::Array{U, 1}, x0::Array{U, 1}, 
                      dx0::Array{U, 1}, Φ::Array{U, 1}, Δs::U, n::T) where {U<:Real, T<:Integer}

    for i in 1:n-2
        F[i] = xT[i][0][1] - x1[i]
    end 
    F[n-1] = sum((x1[i] - x0[i]) * dx0[i] for i in 1:n-2)
    F[n]   = sum((x1[i] - x0[i]) * Φ[i]   for i in 1:n) - Δs

end

function OrbitJacobian!(J::Array{U, 2}, xT::Array{TaylorN{U}, 1}, 
                        dx0::Array{U, 1}, Φ::Array{U, 1}, n::T) where {U<:Real, T<:Integer}
    for i in 1:n-2
        for j in 1:n
            J[i, j] = xT[i][1][j]
        end
        J[i, i] -= 1.0
        J[n-1, i] = dx0[i]
    end
    J[n, :] = Φ
end