function PDSystem!(F::Array{U, 1}, xT::Array{TaylorN{U}, 2}, 
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
        F[comp3] = x1[comp3] - x1[comp2]
    end

    comp4 = comp3

    for comp1 in 1:m-1
        for comp2 in 1:n-2
            comp3 += 1
            F[comp3] = xT[comp1, comp2][0][1] - x1[comp4 + n-2 + comp3]
        end
    end

    for comp2 in 1:n-2
        comp3 += 1
        F[comp3] = x1[comp3] - x1[comp4 + comp2]
    end

    F[N-1] = sum((x1[i] - x0[i]) * dx0[i] for i in 1:n-2)
    F[N] = sum((x1[i] - x0[i]) * Φ[i]   for i in 1:N) - Δs

    return nothing
end

function PDJacobian!(J::Array{U,2}, xT::Array{TaylorN{U},2}, 
                       dx0::Array{U,1}, Φ::Array{U,1}, n::T, m::T, N::T) where {U<:Real, T<:Integer}
    
    n_minus_2 = n - 2
    oneU = one(U)
    
    @inbounds for comp1 in 1:m
        comp1_offset = (comp1-1)*n_minus_2
        for comp2 in 1:n_minus_2
            comp2_idx = comp1_offset + comp2
            
            if comp1 < m
            # Set the main block
                for comp3 in 1:n_minus_2
                    J[comp2_idx, comp1_offset + comp3] = xT[comp1, comp2][1][comp3]
                end

                J[comp2_idx, comp2_idx + n_minus_2] = -oneU

                J[comp2_idx, N-1] = xT[comp1, comp2][1][n-1]
                J[comp2_idx, N] = xT[comp1, comp2][1][n]
            else
                J[comp2_idx, comp2] = -oneU
                J[comp2_idx, comp2_idx] = oneU
                J[comp2_idx, N-1] = 0.0
                J[comp2_idx, N] = 0.0
            end
            
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