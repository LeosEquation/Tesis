function BiProduct(A::AbstractMatrix{U}, B::AbstractMatrix{U}) where U<:Number
    # Esta función calcula el producto bialternante entre las matrices A y B.
    @assert size(A) == size(B) "A and B must have the same dimensions."
    @assert size(A, 1) == size(A, 2) "A is not a square matrix."
    @assert size(B, 1) == size(B, 2) "B is not a square matrix."

    n = size(A, 1)  # Número de filas (y columnas, ya que A es cuadrada)
    m = (n * (n - 1)) ÷ 2  # Número de elementos en la parte triangular

    # Inicializar la matriz C
    C = Array{U, 2}(undef, m, m)

    ii = 1
    for i in 2:n
        for j in 1:i-1
            kk = 1
            for k in 2:n
                for l in 1:k-1
                    @inbounds C[ii, kk] = 0.5 * (
                        A[i, k] * B[j, l] - B[j, k] * A[i, l] +
                        B[i, k] * A[j, l] - A[j, k] * B[i, l]
                    )
                    kk += 1
                end
            end
            ii += 1
        end
    end

    return C
end


function BiProduct(A::AbstractMatrix{U}, B::AbstractMatrix{T}) where {U<:Number,T<:Number}

    if U<:T
        A_ = T.(A)
        return BiProduct(A_, B)
    else
        B_ = U.(B)
        return BiProduct(A, B_)
    end

end

function BiProduct(A::AbstractMatrix{U}) where U<:Number
    # Esta función calcula el producto bialternante entre las matrices A y B.
    @assert size(A, 1) == size(A, 2) "A is not a square matrix."

    n = size(A, 1)  # Número de filas (y columnas, ya que A es cuadrada)
    m = (n * (n - 1)) ÷ 2  # Número de elementos en la parte triangular

    # Inicializar la matriz C
    C = Array{U, 2}(undef, m, m)

    ii = 1
    for i in 2:n
        for j in 1:i-1
            kk = 1
            for k in 2:n
                for l in 1:k-1
                    @inbounds C[ii, kk] = A[i, k] * A[j, l] - A[j, k] * A[i, l]
                    kk += 1
                end
            end
            ii += 1
        end
    end

    return C
end


function BiProduct(A::AbstractMatrix{U}, I::UniformScaling{Bool}) where U<:Number
    
    @assert size(A, 1) == size(A, 2) "A is not a square matrix."

    n = size(A, 1)  # Número de filas (y columnas, ya que A es cuadrada)
    m = (n * (n - 1)) ÷ 2  # Número de elementos en la parte triangular

    # Inicializar la matriz C
    C = Array{U, 2}(undef, m, m)

    ii = 1
    for i in 2:n
        for j in 1:i-1
            kk = 1
            for k in 2:n
                for l in 1:k-1
                    @inbounds C[ii, kk] = 0.5 * (
                        A[i, k] * I[j, l] - I[j, k] * A[i, l] +
                        I[i, k] * A[j, l] - A[j, k] * I[i, l]
                    )
                    kk += 1
                end
            end
            ii += 1
        end
    end

    return C
end


function BiProduct(a::U, A::AbstractMatrix{U}, I::UniformScaling{Bool}) where U<:Number
    
    @assert size(A, 1) == size(A, 2) "A is not a square matrix."

    n = size(A, 1)  # Número de filas (y columnas, ya que A es cuadrada)
    m = (n * (n - 1)) ÷ 2  # Número de elementos en la parte triangular

    # Inicializar la matriz C
    C = Array{U, 2}(undef, m, m)

    ii = 1
    for i in 2:n
        for j in 1:i-1
            kk = 1
            for k in 2:n
                for l in 1:k-1
                    @inbounds C[ii, kk] = a * 0.5 * (
                        A[i, k] * I[j, l] - I[j, k] * A[i, l] +
                        I[i, k] * A[j, l] - A[j, k] * I[i, l]
                    )
                    kk += 1
                end
            end
            ii += 1
        end
    end

    return C
end

function BiProduct(a::T, A::AbstractMatrix{U}, I::UniformScaling{Bool}) where {U<:Number,T<:Number}

    if U<:T
        A_ = T.(A)
        return BiProduct(a, A_, I)
    else
        a_ = U(a)
        return BiProduct(a_, A, I)
    end

end


function BiProduct!(a::U, A::AbstractMatrix{U}, I::UniformScaling{Bool}, 
                    C::AbstractMatrix{U}, n::T) where {U<:Number, T<:Integer}

    ii = 1
    for i in 2:n
        for j in 1:i-1
            kk = 1
            for k in 2:n
                for l in 1:k-1
                    C[ii, kk] = a * 0.5 * (
                        A[i, k] * I[j, l] - I[j, k] * A[i, l] +
                        I[i, k] * A[j, l] - A[j, k] * I[i, l]
                    )
                    kk += 1
                end
            end
            ii += 1
        end
    end

end