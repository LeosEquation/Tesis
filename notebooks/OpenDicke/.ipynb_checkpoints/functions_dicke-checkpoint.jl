@taylorize function OpenDicke!(dξ, ξ, p, t)

    local ω = p[1]
    local ω₀ = p[2]
    local κ = p[3]
    local λ₋ = p[4]
    local λ₊ = p[5]

    a₁ = ξ[1]
    a₂ = ξ[2]
    x  = ξ[3]
    y  = ξ[4]

    κa₁ = - (κ * a₁)
    ωa₂ =   ω * a₂
    local λ₋₊ = λ₋ - λ₊
    λ₋₊y = λ₋₊ * y
    x² = x^2
    y² = y^2
    xy = (x² + y²) + 1.0
    λ₋₊y_xy = λ₋₊y / xy
    
    dξ[1] = ( κa₁ + ωa₂ ) + λ₋₊y_xy

    ωa₁ = - (ω * a₁)
    κa₂ = - (κ * a₂)
    local λ⁻⁺ = λ₋ + λ₊
    λ⁻⁺x = λ⁻⁺ * x
    λ⁻⁺x_xy = - (λ⁻⁺x / xy)
    
    dξ[2] = ( ωa₁ + κa₂ ) + λ⁻⁺x_xy

    λ₋₊a₂ = λ₋₊ * a₂
    λ₋₊a₂x² = λ₋₊a₂ * x²
    λ⁻⁺a₁ = λ⁻⁺ * a₁
    λ⁻⁺a₁xy = (λ⁻⁺a₁ * x) * y
    λ⁻⁺a₁xy2 = - (2.0 * λ⁻⁺a₁xy)
    local λ₊₋ = - λ₋₊
    λ₊₋a₂ = λ₊₋ * a₂
    λ₊₋a₂y² = λ₊₋a₂ * y²
    λ₋₊a₂ = λ₋₊ * a₂
    ω₀y = ω₀ * y
    
    dξ[3] = ((λ₋₊a₂x² + λ⁻⁺a₁xy2) + (λ₊₋a₂y² + λ₋₊a₂)) + ω₀y

    λ⁻⁺a₁x² = λ⁻⁺a₁ * x²
    λ₋₊a₂ = λ₋₊ * a₂
    λ₋₊a₂xy = (λ₋₊a₂ * x) * y
    λ₋₊a₂xy2 = 2.0 * λ₋₊a₂xy
    λ⁻⁺a₁y² = - (λ⁻⁺a₁ * y²)
    ω₀x = - (ω₀ * x)
    
    dξ[4] = ((λ⁻⁺a₁x² + λ₋₊a₂xy2) + (λ⁻⁺a₁y² - λ⁻⁺a₁)) + ω₀x 

    return dξ

end

@doc """
    OpenDicke!(dξ, ξ, p, t)

This function evaluate the differential `dξ` in the dynamical variable
`ξ`. None of the parameters of p is evaluate as dynamical variable.
""" OpenDicke!


@taylorize function OpenDicke1!(dξ, ξ, p, t)

    local ω = p[1]
    local ω₀ = p[2]
    local κ = p[3]
    local λ₋ = p[4]

    a₁ = ξ[1]
    a₂ = ξ[2]
    x  = ξ[3]
    y  = ξ[4]
    λ₊ = ξ[5]

    κa₁ = - (κ * a₁)
    ωa₂ =   ω * a₂
    λ₋₊ = λ₋ - λ₊
    λ₋₊y = λ₋₊ * y
    x² = x^2
    y² = y^2
    xy = (x² + y²) + 1.0
    λ₋₊y_xy = λ₋₊y / xy
    
    dξ[1] = ( κa₁ + ωa₂ ) + λ₋₊y_xy

    ωa₁ = - (ω * a₁)
    κa₂ = - (κ * a₂)
    λ⁻⁺ = λ₋ + λ₊
    λ⁻⁺x = λ⁻⁺ * x
    λ⁻⁺x_xy = - (λ⁻⁺x / xy)
    
    dξ[2] = ( ωa₁ + κa₂ ) + λ⁻⁺x_xy

    λ₋₊a₂ = λ₋₊ * a₂
    λ₋₊a₂x² = λ₋₊a₂ * x²
    λ⁻⁺a₁ = λ⁻⁺ * a₁
    λ⁻⁺a₁xy = (λ⁻⁺a₁ * x) * y
    λ⁻⁺a₁xy2 = - (2.0 * λ⁻⁺a₁xy)
    λ₊₋ = - λ₋₊
    λ₊₋a₂ = λ₊₋ * a₂
    λ₊₋a₂y² = λ₊₋a₂ * y²
    λ₋₊a₂ = λ₋₊ * a₂
    ω₀y = ω₀ * y
    
    dξ[3] = ((λ₋₊a₂x² + λ⁻⁺a₁xy2) + (λ₊₋a₂y² + λ₋₊a₂)) + ω₀y

    λ⁻⁺a₁x² = λ⁻⁺a₁ * x²
    λ₋₊a₂ = λ₋₊ * a₂
    λ₋₊a₂xy = (λ₋₊a₂ * x) * y
    λ₋₊a₂xy2 = 2.0 * λ₋₊a₂xy
    λ⁻⁺a₁y² = - (λ⁻⁺a₁ * y²)
    ω₀x = - (ω₀ * x)
    
            
    dξ[4] = ((λ⁻⁺a₁x² + λ₋₊a₂xy2) + (λ⁻⁺a₁y² - λ⁻⁺a₁)) + ω₀x 

    dξ[5] = zero(ξ[1])

    return dξ

end

@doc """
    OpenDicke1!(dξ, ξ, p, t)

This function evaluate the differential `dξ` in the dynamical variable
`ξ`. One parameter of p is evaluate as dynamical variable.
""" OpenDicke1!


@taylorize function OpenDicke2!(dξ, ξ, p, t)

    local ω = p[1]
    local ω₀ = p[2]
    local κ = p[3]

    a₁ = ξ[1]
    a₂ = ξ[2]
    x  = ξ[3]
    y  = ξ[4]
    λ₋ = ξ[5]
    λ₊ = ξ[6]

    κa₁ = - (κ * a₁)
    ωa₂ =   ω * a₂
    λ₋₊ = λ₋ - λ₊
    λ₋₊y = λ₋₊ * y
    x² = x^2
    y² = y^2
    xy = (x² + y²) + 1.0
    λ₋₊y_xy = λ₋₊y / xy
    
    dξ[1] = ( κa₁ + ωa₂ ) + λ₋₊y_xy

    ωa₁ = - (ω * a₁)
    κa₂ = - (κ * a₂)
    λ⁻⁺ = λ₋ + λ₊
    λ⁻⁺x = λ⁻⁺ * x
    λ⁻⁺x_xy = - (λ⁻⁺x / xy)
    
    dξ[2] = ( ωa₁ + κa₂ ) + λ⁻⁺x_xy

    λ₋₊a₂ = λ₋₊ * a₂
    λ₋₊a₂x² = λ₋₊a₂ * x²
    λ⁻⁺a₁ = λ⁻⁺ * a₁
    λ⁻⁺a₁xy = (λ⁻⁺a₁ * x) * y
    λ⁻⁺a₁xy2 = - (2.0 * λ⁻⁺a₁xy)
    λ₊₋ = - λ₋₊
    λ₊₋a₂ = λ₊₋ * a₂
    λ₊₋a₂y² = λ₊₋a₂ * y²
    λ₋₊a₂ = λ₋₊ * a₂
    ω₀y = ω₀ * y
    
    dξ[3] = ((λ₋₊a₂x² + λ⁻⁺a₁xy2) + (λ₊₋a₂y² + λ₋₊a₂)) + ω₀y

    λ⁻⁺a₁x² = λ⁻⁺a₁ * x²
    λ₋₊a₂ = λ₋₊ * a₂
    λ₋₊a₂xy = (λ₋₊a₂ * x) * y
    λ₋₊a₂xy2 = 2.0 * λ₋₊a₂xy
    λ⁻⁺a₁y² = - (λ⁻⁺a₁ * y²)
    ω₀x = - (ω₀ * x)
    
            
    dξ[4] = ((λ⁻⁺a₁x² + λ₋₊a₂xy2) + (λ⁻⁺a₁y² - λ⁻⁺a₁)) + ω₀x 

    dξ[5] = zero(ξ[1])
    dξ[6] = zero(ξ[1])

    return dξ

end

@doc """
    OpenDicke2!(dξ, ξ, p, t)

This function evaluate the differential `dξ` in the dynamical variable
`ξ`. Two parameters of p are evaluate as dynamical variable.
""" OpenDicke2!

@taylorize function OpenDickeP!(dξ, ξ, p, t)

    local ω = p[1]
    local ω₀ = p[2]
    local κ = p[3]
    local λ₋ = p[4]
    local λ₊ = p[5]
    local T  = p[6]

    a₁ = ξ[1]
    a₂ = ξ[2]
    x  = ξ[3]
    y  = ξ[4]
   
    κa₁ = - (κ * a₁)
    ωa₂ =   ω * a₂
    λ₋₊ = λ₋ - λ₊
    λ₋₊y = λ₋₊ * y
    x² = x^2
    y² = y^2
    xy = (x² + y²) + 1.0
    λ₋₊y_xy = λ₋₊y / xy

    dξ1 = ( κa₁ + ωa₂ ) + λ₋₊y_xy
    
    dξ[1] = T * dξ1

    ωa₁ = - (ω * a₁)
    κa₂ = - (κ * a₂)
    λ⁻⁺ = λ₋ + λ₊
    λ⁻⁺x = λ⁻⁺ * x
    λ⁻⁺x_xy = - (λ⁻⁺x / xy)
    
    dξ2 = ( ωa₁ + κa₂ ) + λ⁻⁺x_xy

    dξ[2] = T * dξ2

    λ₋₊a₂ = λ₋₊ * a₂
    λ₋₊a₂x² = λ₋₊a₂ * x²
    λ⁻⁺a₁ = λ⁻⁺ * a₁
    λ⁻⁺a₁xy = (λ⁻⁺a₁ * x) * y
    λ⁻⁺a₁xy2 = - (2.0 * λ⁻⁺a₁xy)
    λ₊₋ = - λ₋₊
    λ₊₋a₂ = λ₊₋ * a₂
    λ₊₋a₂y² = λ₊₋a₂ * y²
    λ₋₊a₂ = λ₋₊ * a₂
    ω₀y = ω₀ * y

    dξ3 = ((λ₋₊a₂x² + λ⁻⁺a₁xy2) + (λ₊₋a₂y² + λ₋₊a₂)) + ω₀y
    
    dξ[3] = T * dξ3

    λ⁻⁺a₁x² = λ⁻⁺a₁ * x²
    λ₋₊a₂ = λ₋₊ * a₂
    λ₋₊a₂xy = (λ₋₊a₂ * x) * y
    λ₋₊a₂xy2 = 2.0 * λ₋₊a₂xy
    λ⁻⁺a₁y² = - (λ⁻⁺a₁ * y²)
    ω₀x = - (ω₀ * x)
            
    dξ4 = ((λ⁻⁺a₁x² + λ₋₊a₂xy2) + (λ⁻⁺a₁y² - λ⁻⁺a₁)) + ω₀x 

    dξ[4] = T *  dξ4
    
    return dξ

end

@doc """
    OpenDickeP!(dξ, ξ, p, t)

This function evaluate the differential `dξ` in the dynamical variable
`ξ`. None parameter of p is evaluate as dynamical variable and the
maximum time of integration is agree in params p.
""" OpenDickeP!


@taylorize function OpenDickeP1!(dξ, ξ, p, t)
    
    local ω = p[1]
    local ω₀ = p[2]
    local κ = p[3]
    local λ₋ = p[4]

    a₁ = ξ[1]
    a₂ = ξ[2]
    x  = ξ[3]
    y  = ξ[4]
    λ₊ = ξ[5]
    T  = ξ[6]

    κa₁ = - (κ * a₁)
    ωa₂ =   ω * a₂
    λ₋₊ = λ₋ - λ₊
    λ₋₊y = λ₋₊ * y
    x² = x^2
    y² = y^2
    xy = (x² + y²) + 1.0
    λ₋₊y_xy = λ₋₊y / xy

    dξ1 = ( κa₁ + ωa₂ ) + λ₋₊y_xy
    
    dξ[1] = T * dξ1

    ωa₁ = - (ω * a₁)
    κa₂ = - (κ * a₂)
    λ⁻⁺ = λ₋ + λ₊
    λ⁻⁺x = λ⁻⁺ * x
    λ⁻⁺x_xy = - (λ⁻⁺x / xy)
    
    dξ2 = ( ωa₁ + κa₂ ) + λ⁻⁺x_xy

    dξ[2] = T * dξ2

    λ₋₊a₂ = λ₋₊ * a₂
    λ₋₊a₂x² = λ₋₊a₂ * x²
    λ⁻⁺a₁ = λ⁻⁺ * a₁
    λ⁻⁺a₁xy = (λ⁻⁺a₁ * x) * y
    λ⁻⁺a₁xy2 = - (2.0 * λ⁻⁺a₁xy)
    λ₊₋ = - λ₋₊
    λ₊₋a₂ = λ₊₋ * a₂
    λ₊₋a₂y² = λ₊₋a₂ * y²
    λ₋₊a₂ = λ₋₊ * a₂
    ω₀y = ω₀ * y

    dξ3 = ((λ₋₊a₂x² + λ⁻⁺a₁xy2) + (λ₊₋a₂y² + λ₋₊a₂)) + ω₀y
    
    dξ[3] = T * dξ3

    λ⁻⁺a₁x² = λ⁻⁺a₁ * x²
    λ₋₊a₂ = λ₋₊ * a₂
    λ₋₊a₂xy = (λ₋₊a₂ * x) * y
    λ₋₊a₂xy2 = 2.0 * λ₋₊a₂xy
    λ⁻⁺a₁y² = - (λ⁻⁺a₁ * y²)
    ω₀x = - (ω₀ * x)
            
    dξ4 = ((λ⁻⁺a₁x² + λ₋₊a₂xy2) + (λ⁻⁺a₁y² - λ⁻⁺a₁)) + ω₀x 

    dξ[4] = T *  dξ4

    dξ[5] = zero(ξ[1])

    dξ[6] = zero(ξ[1])

    return dξ

end

@doc """
    OpenDickeP1!(dξ, ξ, p, t)

This function evaluate the differential `dξ` in the dynamical variable
`ξ`. Two parameter of p is evaluate as dynamical variable, where one of this
is the maximum time of integration.
""" OpenDickeP1!


@taylorize function OpenDickeP2!(dξ, ξ, p, t)
    
    local ω = p[1]
    local ω₀ = p[2]
    local κ = p[3]
    
    a₁ = ξ[1]
    a₂ = ξ[2]
    x  = ξ[3]
    y  = ξ[4]
    λ₋ = ξ[5]
    λ₊ = ξ[6]
    T  = ξ[7]

    κa₁ = - (κ * a₁)
    ωa₂ =   ω * a₂
    λ₋₊ = λ₋ - λ₊
    λ₋₊y = λ₋₊ * y
    x² = x^2
    y² = y^2
    xy = (x² + y²) + 1.0
    λ₋₊y_xy = λ₋₊y / xy

    dξ1 = ( κa₁ + ωa₂ ) + λ₋₊y_xy
    
    dξ[1] = T * dξ1

    ωa₁ = - (ω * a₁)
    κa₂ = - (κ * a₂)
    λ⁻⁺ = λ₋ + λ₊
    λ⁻⁺x = λ⁻⁺ * x
    λ⁻⁺x_xy = - (λ⁻⁺x / xy)
    
    dξ2 = ( ωa₁ + κa₂ ) + λ⁻⁺x_xy

    dξ[2] = T * dξ2

    λ₋₊a₂ = λ₋₊ * a₂
    λ₋₊a₂x² = λ₋₊a₂ * x²
    λ⁻⁺a₁ = λ⁻⁺ * a₁
    λ⁻⁺a₁xy = (λ⁻⁺a₁ * x) * y
    λ⁻⁺a₁xy2 = - (2.0 * λ⁻⁺a₁xy)
    λ₊₋ = - λ₋₊
    λ₊₋a₂ = λ₊₋ * a₂
    λ₊₋a₂y² = λ₊₋a₂ * y²
    λ₋₊a₂ = λ₋₊ * a₂
    ω₀y = ω₀ * y

    dξ3 = ((λ₋₊a₂x² + λ⁻⁺a₁xy2) + (λ₊₋a₂y² + λ₋₊a₂)) + ω₀y
    
    dξ[3] = T * dξ3

    λ⁻⁺a₁x² = λ⁻⁺a₁ * x²
    λ₋₊a₂ = λ₋₊ * a₂
    λ₋₊a₂xy = (λ₋₊a₂ * x) * y
    λ₋₊a₂xy2 = 2.0 * λ₋₊a₂xy
    λ⁻⁺a₁y² = - (λ⁻⁺a₁ * y²)
    ω₀x = - (ω₀ * x)
            
    dξ4 = ((λ⁻⁺a₁x² + λ₋₊a₂xy2) + (λ⁻⁺a₁y² - λ⁻⁺a₁)) + ω₀x 

    dξ[4] = T *  dξ4

    dξ[5] = zero(ξ[1])

    dξ[6] = zero(ξ[1])

    return dξ

end

@doc """
    OpenDickeP2!(dξ, ξ, p, t)

This function evaluate the differential `dξ` in the dynamical variable
`ξ`. Three parameter of p is evaluate as dynamical variable, where one of this
is the maximum time of integration.
""" OpenDickeP2!