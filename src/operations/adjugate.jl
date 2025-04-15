function adjugate(A::AbstractMatrix)
    ishermitian(A) && return adjugate(Hermitian(A))
    LinearAlgebra.checksquare(A)
    F = svd(A)
    return F.V * (adjugate(Diagonal(F.S)) * det(F.Vt * F.U)) * F.U'
end

function adjugate(A::LinearAlgebra.RealHermSymComplexHerm)
    F = eigen(A)
    return F.vectors * adjugate(Diagonal(F.values)) * F.vectors'
end

function adjugate(D::Diagonal{<:Number})
    d = D.diag
    Base.require_one_based_indexing(d)
    length(d) < 2 && return Diagonal(one.(d))

    # compute diagonal of adj(D): dadj[i] = prod(d) / d[i], but avoid dividing by zero
    dadj = similar(d)
    prod = one(eltype(d))
    for i = 1:length(d)
        dadj[i] = prod
        prod *= d[i]
    end
    prod = one(eltype(d))
    for i = length(d):-1:1
        dadj[i] *= prod
        prod *= d[i]
    end
    
    return Diagonal(dadj)
end