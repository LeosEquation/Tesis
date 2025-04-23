function BPJacobian!(J::Array{U, 2}, Jeval::Array{U, 2}, q::Array{U, 1}, 
                     dx::Array{TaylorN{U}, 1}, Φ::Array{U, 1}, n::T) where {U<:Real, T<:Integer}

    for j in 1:n 
        for i in 1:n-1
            J[i, j] = Jeval[i, j]
            J[n-1+i, j] = sum(differentiate(ntuple(l -> count(==(l), (i, j)), n), dx[k]) * q[n+k] for k in 1:n-1)
            J[n+i-1, n+j] = Jeval[j, i]
        end
        J[2*n-1, j] = sum(differentiate(ntuple(l -> count(==(l), (n, j)), n), dx[k]) * q[n+k] for k in 1:n-1)
        J[2*n-1, n+j] = Jeval[j, n]
    end

    for j in 1:n-1
        J[j, n+j] = q[2*n]
        J[j, 2*n] = q[n+j]
        J[2*n, n+j] = 2.0 * q[n+j]
    end

    J[2*n+1, :] .= Φ
end

function BPSystem!(F::Array{U,1}, dxeval::Array{U, 1}, 
                   q::Array{U, 1}, q0::Array{U, 1}, 
                   Φ::Array{U, 1}, Jeval::Array{U, 2}, Δs::U, n::T) where {U<:Real, T<:Integer}

    for i in 1:n-1
        F[i] = dxeval[i] + q[2*n] * q[n+i]
        F[n-1+i] = sum(q[n+k] * Jeval[k, i] for k in 1:n-1)
    end

    F[2*n-1] = sum(q[n+k] * Jeval[k, n] for k in 1:n-1)
    F[2*n] = sum(q[n+k]^2 for k in 1:n-1) - 1.0

    F[2*n+1] = sum((q[k] - q0[k]) * Φ[k] for k in 1:(2*n+1)) - Δs

end


function BPJacobian!(J::Array{U, 2}, Jeval::Array{U, 2}, q::Array{U, 1}, 
                     dx::Array{TaylorN{U}, 1}, n::T) where {U<:Real, T<:Integer}

    for j in 1:n 
        for i in 1:n-1
            J[i, j] = Jeval[i, j]
            J[n-1+i, j] = sum(differentiate(ntuple(l -> count(==(l), (i, j)), n), dx[k]) * q[n+k] for k in 1:n-1)
            J[n+i-1, n+j] = Jeval[j, i]
        end
        J[2*n-1, j] = sum(differentiate(ntuple(l -> count(==(l), (n, j)), n), dx[k]) * q[n+k] for k in 1:n-1)
        J[2*n-1, n+j] = Jeval[j, n]
    end

    for j in 1:n-1
        J[j, n+j] = q[2*n]
        J[j, 2*n] = q[n+j]
        J[2*n, n+j] = 2.0 * q[n+j]
    end

end

function BPSystem!(F::Array{U,1}, dxeval::Array{U, 1}, 
                   q::Array{U, 1}, Jeval::Array{U, 2}, n::T) where {U<:Real, T<:Integer}

    for i in 1:n-1
        F[i] = dxeval[i] + q[2*n] * q[n+i]
        F[n-1+i] = sum(q[n+k] * Jeval[k, i] for k in 1:n-1)
    end

    F[2*n-1] = sum(q[n+k] * Jeval[k, n] for k in 1:n-1)
    F[2*n] = sum(q[n+k]^2 for k in 1:n-1) - 1.0

end