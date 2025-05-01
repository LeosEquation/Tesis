
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

function OrbitSystem!(F::Array{U, 1}, xT::Array{TaylorN{U}, 2}, 
                      x1::Array{U, 1}, x0::Array{U, 1}, 
                      dx0::Array{U, 1}, Φ::Array{U, 1}, Δs::U, n::T, m::T, N::T) where {U<:Real, T<:Integer}

    comp3 = 0        
    for comp1 in 1:m-1
        for comp2 in 1:n-2
            comp3 += 1
            F[comp3] = xT[comp1, comp2][0][1] - x1[n-2 + comp3]
        end
    end

    for comp2 in 1:n-2
        comp3 += 1
        F[comp3] = xT[m, comp2][0][1] - x1[comp2]
    end

    F[N-1] = sum((x1[i] - x0[i]) * dx0[i] for i in 1:n-2)
    F[N] = sum((x1[i] - x0[i]) * Φ[i]   for i in 1:N) - Δs

    return nothing
end

function OrbitJacobian!(J::Array{U,2}, xT::Array{TaylorN{U},2}, 
                       dx0::Array{U,1}, Φ::Array{U,1}, n::T, m::T, N::T) where {U<:Real, T<:Integer}
    
    n_minus_2 = n - 2
    oneU = one(U)
    
    @inbounds for comp1 in 1:m
        comp1_offset = (comp1-1)*n_minus_2
        for comp2 in 1:n_minus_2
            comp2_idx = comp1_offset + comp2
            
            # Set the main block
            for comp3 in 1:n_minus_2
                J[comp2_idx, comp1_offset + comp3] = xT[comp1, comp2][1][comp3]
            end
            
            # Set the -1 elements
            if comp1 < m
                J[comp2_idx, comp2_idx + n_minus_2] = -oneU
            else 
                J[comp2_idx, comp2] = -oneU
            end
            
            # Set the last two columns
            J[comp2_idx, N-1] = xT[comp1, comp2][1][n-1]
            J[comp2_idx, N] = xT[comp1, comp2][1][n]
        end
    end
    
    # Set the second-to-last row
    @inbounds for i in 1:n_minus_2
        J[N-1, i] = dx0[i]
    end
    
    # Set the last row
    @inbounds J[N, :] = Φ
    
    return nothing
end