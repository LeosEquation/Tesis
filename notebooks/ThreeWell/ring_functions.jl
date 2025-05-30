################################################################
# Funciones sin parámetro variable

function H(ξ, params, t)

    local γ = params[1]
    local α = params[2]
    local β = params[3]
    local K = params[4]

    J₂ = ξ[1]
    ψ₂ = ξ[2]
    ψ₁ = ξ[3]
    J₁ = ξ[4]
    
    J₃ = ( K - J₁ ) - J₂

    H₀ = γ * K

    V  = α * ( ( ( J₁ ^ 2 ) + ( J₂ ^ 2 ) ) + ( J₃ ^ 2 ) )

    W  = ( ( sqrt(J₁ * J₃) * cos(ψ₁) ) + ( sqrt(J₂ * J₃) * cos(ψ₂) ) ) + ( β * ( sqrt(J₂ * J₁) * cos(ψ₂ - ψ₁) ) )

    return ( H₀ + V ) - W

end

@taylorize function ring!(dξ, ξ, params, t)
    
    local γ = params[1]
    local α = params[2]
    local β = params[3]
    local K = params[4]

    J₂ = ξ[1]
    ψ₂ = ξ[2]
    ψ₁ = ξ[3]
    J₁ = ξ[4]

    ψ₂₁ = ψ₂ - ψ₁
    
    J₃ = ( K - J₁ ) - J₂
    
    sqrtJ₁_J₃ = sqrt(J₁ / J₃)
    sqrtJ₂_J₃ = sqrt(J₂ / J₃)

    sqrtJ₁J₂ = sqrt(J₁ * J₂)

    cosψ₁  = cos(ψ₁)
    cosψ₂  = cos(ψ₂)
    cosψ₂₁ = cos(ψ₂₁)

    sinψ₂₁ = sin(ψ₂₁)

    local α2 = 2.0 * α
    
    local β_2 = 0.5 * β

    dξ[1] = ( sqrt(J₂ * J₃) * sin(ψ₂) ) + ( β * ( sqrtJ₁J₂ * sinψ₂₁ ) )

    dξ[2] = ( α2  * ( J₃ - J₂ ) ) + 
            ( ( 0.5 * ( ( ( sqrt(J₃ / J₂) - sqrtJ₂_J₃ ) *  cosψ₂ ) - ( sqrtJ₁_J₃ * cosψ₁ ) ) ) +
              ( β_2 * ( sqrt(J₁ / J₂) * cosψ₂₁ ) ) )

    dξ[3] = ( α2  * ( J₃ - J₁ ) ) + 
            ( ( 0.5 * ( ( ( sqrt(J₃ / J₁) - sqrtJ₁_J₃ ) *  cosψ₁ ) - ( sqrtJ₂_J₃ * cosψ₂ ) ) )  +
              ( β_2 * ( sqrt(J₂ / J₁) * cosψ₂₁ ) ) )

    dξ[4] = ( sqrt(J₁ * J₃) * sin(ψ₁) ) - ( β * ( sqrtJ₁J₂ * sinψ₂₁ ) )

    return nothing
    
end

################################################################
# Funciones sin parámetro variable

function H_T(ξ, params, t)

    local γ = params[1]
    local α = params[2]
    local β = params[3]
    local K = params[4]

    J₂ = ξ[1]
    ψ₂ = ξ[2]
    ψ₁ = ξ[3]
    J₁ = ξ[4]
    T  = ξ[5]
    
    J₃ = ( K - J₁ ) - J₂

    H₀ = γ * K

    V  = α * ( ( ( J₁ ^ 2 ) + ( J₂ ^ 2 ) ) + ( J₃ ^ 2 ) )

    W  = ( ( sqrt(J₁ * J₃) * cos(ψ₁) ) + ( sqrt(J₂ * J₃) * cos(ψ₂) ) ) + ( β * ( sqrt(J₂ * J₁) * cos(ψ₂ - ψ₁) ) )

    return ( H₀ + V ) - W

end

@taylorize function ring_T!(dξ, ξ, params, t)

    local γ = params[1]
    local α = params[2]
    local β = params[3]
    local K = params[4]

    J₂ = ξ[1]
    ψ₂ = ξ[2]
    ψ₁ = ξ[3]
    J₁ = ξ[4]
    T  = ξ[5]

    ψ₂₁ = ψ₂ - ψ₁
    
    J₃ = ( K - J₁ ) - J₂
    
    sqrtJ₁_J₃ = sqrt(J₁ / J₃)
    sqrtJ₂_J₃ = sqrt(J₂ / J₃)

    sqrtJ₁J₂ = sqrt(J₁ * J₂)

    cosψ₁  = cos(ψ₁)
    cosψ₂  = cos(ψ₂)
    cosψ₂₁ = cos(ψ₂₁)

    sinψ₂₁ = sin(ψ₂₁)

    local α2 = 2.0 * α
    
    local β_2 = 0.5 * β

    dξ[1] = T * ( ( sqrt(J₂ * J₃) * sin(ψ₂) ) + ( β * ( sqrtJ₁J₂ * sinψ₂₁ ) ) )

    dξ[2] = T * ( ( α2  * ( J₃ - J₂ ) ) + 
                  ( ( 0.5 * ( ( ( sqrt(J₃ / J₂) - sqrtJ₂_J₃ ) *  cosψ₂ ) - ( sqrtJ₁_J₃ * cosψ₁ ) ) ) +
                    ( β_2 * ( sqrt(J₁ / J₂) * cosψ₂₁ ) ) ) )

    dξ[3] = T * ( ( α2  * ( J₃ - J₁ ) ) + 
                  ( ( 0.5 * ( ( ( sqrt(J₃ / J₁) - sqrtJ₁_J₃ ) *  cosψ₁ ) - ( sqrtJ₂_J₃ * cosψ₂ ) ) )  +
                    ( β_2 * ( sqrt(J₂ / J₁) * cosψ₂₁ ) ) ) )

    dξ[4] = T * ( ( sqrt(J₁ * J₃) * sin(ψ₁) ) - ( β * ( sqrtJ₁J₂ * sinψ₂₁ ) ) )

    dξ[5] = zero(ξ[1])

    return nothing
    
end

################################################################
# Funciones con parámetro β variable

function H_β(ξ, params, t)

    local γ = params[1]
    local α = params[2]
    local K = params[3]

    J₂ = ξ[1]
    ψ₂ = ξ[2]
    ψ₁ = ξ[3]
    J₁ = ξ[4]
    β  = ξ[5]
    
    J₃ = ( K - J₁ ) - J₂

    H₀ = γ * K

    V  = α * ( ( ( J₁ ^ 2 ) + ( J₂ ^ 2 ) ) + ( J₃ ^ 2 ) )

    W  = ( ( sqrt(J₁ * J₃) * cos(ψ₁) ) + ( sqrt(J₂ * J₃) * cos(ψ₂) ) ) + ( β * ( sqrt(J₂ * J₁) * cos(ψ₂ - ψ₁) ) )

    return ( H₀ + V ) - W

end

@taylorize function ring_β!(dξ, ξ, params, t)
    
    local γ = params[1]
    local α = params[2]
    local K = params[3]

    J₂ = ξ[1]
    ψ₂ = ξ[2]
    ψ₁ = ξ[3]
    J₁ = ξ[4]
    β  = ξ[5]

    ψ₂₁ = ψ₂ - ψ₁
    
    J₃ = ( K - J₁ ) - J₂
    
    sqrtJ₁_J₃ = sqrt(J₁ / J₃)
    sqrtJ₂_J₃ = sqrt(J₂ / J₃)

    sqrtJ₁J₂ = sqrt(J₁ * J₂)

    cosψ₁  = cos(ψ₁)
    cosψ₂  = cos(ψ₂)
    cosψ₂₁ = cos(ψ₂₁)

    sinψ₂₁ = sin(ψ₂₁)

    local α2 = 2.0 * α
    
    β_2 = 0.5 * β

    dξ[1] = ( sqrt(J₂ * J₃) * sin(ψ₂) ) + ( β * ( sqrtJ₁J₂ * sinψ₂₁ ) )

    dξ[2] = ( α2  * ( J₃ - J₂ ) ) + 
            ( 0.5 * ( ( ( sqrt(J₃ / J₂) - sqrtJ₂_J₃ ) *  cosψ₂ ) - ( sqrtJ₁_J₃ * cosψ₁ ) ) +
              ( β_2 * ( sqrt(J₁ / J₂) * cosψ₂₁ ) ) )

    dξ[3] = ( α2  * ( J₃ - J₁ ) ) + 
            ( 0.5 * ( ( ( sqrt(J₃ / J₁) - sqrtJ₁_J₃ ) *  cosψ₁ ) - ( sqrtJ₂_J₃ * cosψ₂ ) )  +
              ( β_2 * ( sqrt(J₂ / J₁) * cosψ₂₁ ) ) )

    dξ[4] = ( sqrt(J₁ * J₃) * sin(ψ₁) ) - ( β * ( sqrtJ₁J₂ * sinψ₂₁ ) )

    dξ[5] = zero(ξ[1])

    return nothing
    
end


################################################################
# Funciones con parámetro β y α variable

function H_βα(ξ, params, t)

    local γ = params[1]
    local K = params[2]

    J₂ = ξ[1]
    ψ₂ = ξ[2]
    ψ₁ = ξ[3]
    J₁ = ξ[4]
    β  = ξ[5]
    α  = ξ[6]
    
    J₃ = ( K - J₁ ) - J₂

    H₀ = γ * K

    V  = α * ( ( ( J₁ ^ 2 ) + ( J₂ ^ 2 ) ) + ( J₃ ^ 2 ) )

    W  = ( ( sqrt(J₁ * J₃) * cos(ψ₁) ) + ( sqrt(J₂ * J₃) * cos(ψ₂) ) ) + ( β * ( sqrt(J₂ * J₁) * cos(ψ₂ - ψ₁) ) )

    return ( H₀ + V ) - W

end

@taylorize function ring_βα!(dξ, ξ, params, t)
    
    local γ = params[1]
    local K = params[2]

    J₂ = ξ[1]
    ψ₂ = ξ[2]
    ψ₁ = ξ[3]
    J₁ = ξ[4]
    β  = ξ[5]
    α  = ξ[6]

    ψ₂₁ = ψ₂ - ψ₁
    
    J₃ = ( K - J₁ ) - J₂
    
    sqrtJ₁_J₃ = sqrt(J₁ / J₃)
    sqrtJ₂_J₃ = sqrt(J₂ / J₃)

    sqrtJ₁J₂ = sqrt(J₁ * J₂)

    cosψ₁  = cos(ψ₁)
    cosψ₂  = cos(ψ₂)
    cosψ₂₁ = cos(ψ₂₁)

    sinψ₂₁ = sin(ψ₂₁)

    α2 = 2.0 * α
    
    β_2 = 0.5 * β

    dξ[1] = ( sqrt(J₂ * J₃) * sin(ψ₂) ) + ( β * ( sqrtJ₁J₂ * sinψ₂₁ ) )

    dξ[2] = ( α2  * ( J₃ - J₂ ) ) + 
            ( 0.5 * ( ( ( sqrt(J₃ / J₂) - sqrtJ₂_J₃ ) *  cosψ₂ ) - ( sqrtJ₁_J₃ * cosψ₁ ) ) +
              ( β_2 * ( sqrt(J₁ / J₂) * cosψ₂₁ ) ) )

    dξ[3] = ( α2  * ( J₃ - J₁ ) ) + 
            ( 0.5 * ( ( ( sqrt(J₃ / J₁) - sqrtJ₁_J₃ ) *  cosψ₁ ) - ( sqrtJ₂_J₃ * cosψ₂ ) )  +
              ( β_2 * ( sqrt(J₂ / J₁) * cosψ₂₁ ) ) )

    dξ[4] = ( sqrt(J₁ * J₃) * sin(ψ₁) ) - ( β * ( sqrtJ₁J₂ * sinψ₂₁ ) )

    dξ[5] = zero(ξ[1])

    dξ[6] = zero(ξ[1])

    return nothing
    
end

################################################################
# Funciones con parámetro λ y periodo T variable

function H_λT(ξ, params, t)

    local γ = params[1]
    local α = params[2]
    local β = params[3]
    local K = params[4]

    J₂ = ξ[1]
    ψ₂ = ξ[2]
    ψ₁ = ξ[3]
    J₁ = ξ[4]
    λ  = ξ[5]
    T  = ξ[6]
    
    J₃ = ( K - J₁ ) - J₂

    H₀ = γ * K

    V  = α * ( ( ( J₁ ^ 2 ) + ( J₂ ^ 2 ) ) + ( J₃ ^ 2 ) )

    W  = ( ( sqrt(J₁ * J₃) * cos(ψ₁) ) + ( sqrt(J₂ * J₃) * cos(ψ₂) ) ) + ( β * ( sqrt(J₂ * J₁) * cos(ψ₂ - ψ₁) ) )

    return ( H₀ + V ) - W

end

@taylorize function ring_λT!(dξ, ξ, params, t)
    
    local γ = params[1]
    local α = params[2]
    local β = params[3]
    local K = params[4]

    J₂ = ξ[1]
    ψ₂ = ξ[2]
    ψ₁ = ξ[3]
    J₁ = ξ[4]
    λ  = ξ[5]
    T  = ξ[6]

    ψ₂₁ = ψ₂ - ψ₁
    
    J₃ = ( K - J₁ ) - J₂
    
    sqrtJ₁_J₃ = sqrt(J₁ / J₃)
    sqrtJ₂_J₃ = sqrt(J₂ / J₃)

    sqrtJ₁J₂ = sqrt(J₁ * J₂)

    cosψ₁  = cos(ψ₁)
    cosψ₂  = cos(ψ₂)
    cosψ₂₁ = cos(ψ₂₁)

    sinψ₂₁ = sin(ψ₂₁)

    local α2 = 2.0 * α
    
    local β_2 = 0.5 * β

    dξ11 = ( sqrt(J₂ * J₃) * sin(ψ₂) ) + ( β * ( sqrtJ₁J₂ * sinψ₂₁ ) )

    dξ21 = ( α2  * ( J₃ - J₂ ) ) + 
           ( ( 0.5 * ( ( ( sqrt(J₃ / J₂) - sqrtJ₂_J₃ ) *  cosψ₂ ) - ( sqrtJ₁_J₃ * cosψ₁ ) ) ) +
             ( β_2 * ( sqrt(J₁ / J₂) * cosψ₂₁ ) ) ) 

    dξ31 = ( α2  * ( J₃ - J₁ ) ) + 
           ( ( 0.5 * ( ( ( sqrt(J₃ / J₁) - sqrtJ₁_J₃ ) *  cosψ₁ ) - ( sqrtJ₂_J₃ * cosψ₂ ) ) ) +
             ( β_2 * ( sqrt(J₂ / J₁) * cosψ₂₁ ) ) )

    dξ41 = ( sqrt(J₁ * J₃) * sin(ψ₁) ) - ( β * ( sqrtJ₁J₂ * sinψ₂₁ ) )

    dξ12 = dξ21
    dξ22 = - dξ11
    dξ32 = - dξ41
    dξ42 = dξ31

    dξ[1] = T * ( dξ11 + ( λ * dξ12 ) )

    dξ[2] = T * ( dξ21 + ( λ * dξ22 ) )

    dξ[3] = T * ( dξ31 + ( λ * dξ32 ) )

    dξ[4] = T * ( dξ41 + ( λ * dξ42 ) ) 

    dξ[5] = zero(ξ[1])

    dξ[6] = zero(ξ[1])

    return nothing
    
end

################################################################
# Funciones con parámetro λ y periodo T variable

function H_λ(ξ, params, t)

    local γ = params[1]
    local α = params[2]
    local β = params[3]
    local K = params[4]

    J₂ = ξ[1]
    ψ₂ = ξ[2]
    ψ₁ = ξ[3]
    J₁ = ξ[4]
    λ  = ξ[5]
    
    J₃ = ( K - J₁ ) - J₂

    H₀ = γ * K

    V  = α * ( ( ( J₁ ^ 2 ) + ( J₂ ^ 2 ) ) + ( J₃ ^ 2 ) )

    W  = ( ( sqrt(J₁ * J₃) * cos(ψ₁) ) + ( sqrt(J₂ * J₃) * cos(ψ₂) ) ) + ( β * ( sqrt(J₂ * J₁) * cos(ψ₂ - ψ₁) ) )

    return ( H₀ + V ) - W

end

@taylorize function ring_λ!(dξ, ξ, params, t)
    
    local γ = params[1]
    local α = params[2]
    local β = params[3]
    local K = params[4]

    J₂ = ξ[1]
    ψ₂ = ξ[2]
    ψ₁ = ξ[3]
    J₁ = ξ[4]
    λ  = ξ[5]

    ψ₂₁ = ψ₂ - ψ₁
    
    J₃ = ( K - J₁ ) - J₂
    
    sqrtJ₁_J₃ = sqrt(J₁ / J₃)
    sqrtJ₂_J₃ = sqrt(J₂ / J₃)

    sqrtJ₁J₂ = sqrt(J₁ * J₂)

    cosψ₁  = cos(ψ₁)
    cosψ₂  = cos(ψ₂)
    cosψ₂₁ = cos(ψ₂₁)

    sinψ₂₁ = sin(ψ₂₁)

    local α2 = 2.0 * α
    
    local β_2 = 0.5 * β

    dξ11 = ( sqrt(J₂ * J₃) * sin(ψ₂) ) + ( β * ( sqrtJ₁J₂ * sinψ₂₁ ) )

    dξ21 = ( α2  * ( J₃ - J₂ ) ) + 
           ( ( 0.5 * ( ( ( sqrt(J₃ / J₂) - sqrtJ₂_J₃ ) *  cosψ₂ ) - ( sqrtJ₁_J₃ * cosψ₁ ) ) ) +
             ( β_2 * ( sqrt(J₁ / J₂) * cosψ₂₁ ) ) ) 

    dξ31 = ( α2  * ( J₃ - J₁ ) ) + 
           ( ( 0.5 * ( ( ( sqrt(J₃ / J₁) - sqrtJ₁_J₃ ) *  cosψ₁ ) - ( sqrtJ₂_J₃ * cosψ₂ ) ) ) +
             ( β_2 * ( sqrt(J₂ / J₁) * cosψ₂₁ ) ) )

    dξ41 = ( sqrt(J₁ * J₃) * sin(ψ₁) ) - ( β * ( sqrtJ₁J₂ * sinψ₂₁ ) )

    dξ12 = dξ21
    dξ22 = - dξ11
    dξ32 = - dξ41
    dξ42 = dξ31

    dξ[1] = dξ11 + ( λ * dξ12 )

    dξ[2] = dξ21 + ( λ * dξ22 )

    dξ[3] = dξ31 + ( λ * dξ32 )

    dξ[4] = dξ41 + ( λ * dξ42 )

    dξ[5] = zero(ξ[1])

    return nothing
    
end

################################################################
# Funcion con perturbación

@taylorize function ring_per!(dξ, ξ, params, t)
    
    local γ = params[1]
    local α = params[2]
    local β = params[3]
    local λ = params[4]
    local K = params[5]

    J₂ = ξ[1]
    ψ₂ = ξ[2]
    ψ₁ = ξ[3]
    J₁ = ξ[4]

    ψ₂₁ = ψ₂ - ψ₁
    
    J₃ = ( K - J₁ ) - J₂
    
    sqrtJ₁_J₃ = sqrt(J₁ / J₃)
    sqrtJ₂_J₃ = sqrt(J₂ / J₃)

    sqrtJ₁J₂ = sqrt(J₁ * J₂)

    cosψ₁  = cos(ψ₁)
    cosψ₂  = cos(ψ₂)
    cosψ₂₁ = cos(ψ₂₁)

    sinψ₂₁ = sin(ψ₂₁)

    local α2 = 2.0 * α
    
    local β_2 = 0.5 * β

    dξ11 = ( sqrt(J₂ * J₃) * sin(ψ₂) ) + ( β * ( sqrtJ₁J₂ * sinψ₂₁ ) )

    dξ21 = ( α2  * ( J₃ - J₂ ) ) + 
           ( ( 0.5 * ( ( ( sqrt(J₃ / J₂) - sqrtJ₂_J₃ ) *  cosψ₂ ) - ( sqrtJ₁_J₃ * cosψ₁ ) ) ) +
             ( β_2 * ( sqrt(J₁ / J₂) * cosψ₂₁ ) ) ) 

    dξ31 = ( α2  * ( J₃ - J₁ ) ) + 
           ( ( 0.5 * ( ( ( sqrt(J₃ / J₁) - sqrtJ₁_J₃ ) *  cosψ₁ ) - ( sqrtJ₂_J₃ * cosψ₂ ) ) ) +
             ( β_2 * ( sqrt(J₂ / J₁) * cosψ₂₁ ) ) )

    dξ41 = ( sqrt(J₁ * J₃) * sin(ψ₁) ) - ( β * ( sqrtJ₁J₂ * sinψ₂₁ ) )

    dξ12 = dξ21
    dξ22 = - dξ11
    dξ32 = - dξ41
    dξ42 = dξ31

    dξ[1] = dξ11 + ( λ * dξ12 )

    dξ[2] = dξ21 + ( λ * dξ22 )

    dξ[3] = dξ31 + ( λ * dξ32 )

    dξ[4] = dξ41 + ( λ * dξ42 )

    return nothing
    
end

################################################################
# Funciones con parámetro β y periodo T variable

function H_βT(ξ, params, t)

    local γ = params[1]
    local α = params[2]
    local K = params[3]

    J₂ = ξ[1]
    ψ₂ = ξ[2]
    ψ₁ = ξ[3]
    J₁ = ξ[4]
    β  = ξ[5]
    T  = ξ[6]
    
    J₃ = ( K - J₁ ) - J₂

    H₀ = γ * K

    V  = α * ( ( ( J₁ ^ 2 ) + ( J₂ ^ 2 ) ) + ( J₃ ^ 2 ) )

    W  = ( ( sqrt(J₁ * J₃) * cos(ψ₁) ) + ( sqrt(J₂ * J₃) * cos(ψ₂) ) ) + ( β * ( sqrt(J₂ * J₁) * cos(ψ₂ - ψ₁) ) )

    return ( H₀ + V ) - W

end

@taylorize function ring_βT!(dξ, ξ, params, t)
    
    local γ = params[1]
    local α = params[2]
    local K = params[3]

    J₂ = ξ[1]
    ψ₂ = ξ[2]
    ψ₁ = ξ[3]
    J₁ = ξ[4]
    β  = ξ[5]
    T  = ξ[6]

    ψ₂₁ = ψ₂ - ψ₁
    
    J₃ = ( K - J₁ ) - J₂
    
    sqrtJ₁_J₃ = sqrt(J₁ / J₃)
    sqrtJ₂_J₃ = sqrt(J₂ / J₃)

    sqrtJ₁J₂ = sqrt(J₁ * J₂)

    cosψ₁  = cos(ψ₁)
    cosψ₂  = cos(ψ₂)
    cosψ₂₁ = cos(ψ₂₁)

    sinψ₂₁ = sin(ψ₂₁)

    local α2 = 2.0 * α
    
    β_2 = 0.5 * β

    dξ[1] = T * ( ( sqrt(J₂ * J₃) * sin(ψ₂) ) + ( β * ( sqrtJ₁J₂ * sinψ₂₁ ) ) )

    dξ[2] = T * ( ( α2  * ( J₃ - J₂ ) ) + 
                 ( ( 0.5 * ( ( ( sqrt(J₃ / J₂) - sqrtJ₂_J₃ ) *  cosψ₂ ) - ( sqrtJ₁_J₃ * cosψ₁ ) ) ) +
                   ( β_2 * ( sqrt(J₁ / J₂) * cosψ₂₁ ) ) ) )

    dξ[3] = T * ( ( α2  * ( J₃ - J₁ ) ) + 
                 ( ( 0.5 * ( ( ( sqrt(J₃ / J₁) - sqrtJ₁_J₃ ) *  cosψ₁ ) - ( sqrtJ₂_J₃ * cosψ₂ ) ) ) +
                   ( β_2 * ( sqrt(J₂ / J₁) * cosψ₂₁ ) ) ) )

    dξ[4] = T * ( ( sqrt(J₁ * J₃) * sin(ψ₁) ) - ( β * ( sqrtJ₁J₂ * sinψ₂₁ ) ) )

    dξ[5] = zero(ξ[1])

    dξ[6] = zero(ξ[1])

    return nothing
    
end

