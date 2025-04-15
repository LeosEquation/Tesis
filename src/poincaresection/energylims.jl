


function _CriticalEnergy!(f!, h!, x::Vector{Float64}, params, tol::Float64, ite::Int64, 
                          S::Matrix{Taylor1{Float64}}, J::Matrix{Float64}, dx::Vector{Taylor1{Float64}})

    @inbounds for j in 1:4

        f!(dx, x + S[j, :], params, 0.0)

        for i in 1:4

            J[i, j] = evaluate(differentiate(dx[i]))

        end

    end

    k = 1

    while k <= ite && norm(evaluate(dx)) > tol

        if det(J) == 0.0
            @warn("El jacobiano se volvió singular.")
            return norm(evaluate(dx))
        end

        x[1:4] .-= J \ evaluate(dx[1:4])

        if h!(x, params)
            @warn("El método Newton Raphson salió de la zona permitida.")
            return norm(evaluate(dx))
        end

        @inbounds for j in 1:4

            f!(dx, x + S[j, :], params, 0.0)
    
            for i in 1:4
    
                J[i, j] = evaluate(differentiate(dx[i]))
    
            end
    
        end

        k += 1

    end 
    
    return norm(evaluate(dx))

end


"""
    CriticalEnergy(H, f!, test_function, Q1min::Float64, Q1max::Float64, Pmin::Union{Float64, Irrational}, 
                   Pmax::Union{Float64, Irrational}, Q2lims, params, tol::Float64, ite::Int64, nrand::Int64)

Esta función calcula los valores máximos y mínimos de energía de un hamiltoniano `H` mediante el uso de un método 
iterativo de Newton-Raphson para obtener los puntos críticos de energía. Los puntos críticos corresponden a los 
valores de las coordenadas y momento que corresponden a los máximos y mínimos de energía dentro de un dominio permitido.

### Parámetros:
- `H`: Función que calcula el valor del hamiltoniano en un punto dado.
- `f!`: Función que define las ecuaciones dinámicas para el sistema.
- `test_function`: Función que verifica si un punto está dentro de los límites permitidos definidos por las coordenadas 
   y momentos.
- `Q1min`, `Q1max`: Límites para la coordenada generalizada `Q1`.
- `Pmin`, `Pmax`: Límites para la coordenada generalizada `P`, que pueden ser números irracionales.
- `Q2lims`: Función que determina los límites de la coordenada `Q2` en función de las coordenadas `Q1` y `P`.
- `params`: Parámetros adicionales necesarios para las funciones `H` y `f!`.
- `tol`: Tolerancia para la convergencia del método de Newton-Raphson.
- `ite`: Número máximo de iteraciones permitidas para el método de Newton-Raphson.
- `nrand`: Número de puntos aleatorios generados para la búsqueda de máximos y mínimos.

### Salida:
La función devuelve una tupla con tres elementos:
1. Una tupla con las coordenadas de los puntos donde se encuentran los máximos y mínimos de energía.
2. Una tupla con los valores de la energía en esos puntos.
3. Una tupla con los valores de la norma de la diferencia entre las soluciones obtenidas y los valores exactos de la energía.

### Descripción:
La función genera puntos aleatorios dentro de los límites definidos por `Q1min`, `Q1max`, `Pmin`, `Pmax`, y `Q2lims`. Luego, 
calcula el valor del hamiltoniano en esos puntos y encuentra el punto con la energía mínima y máxima. Para cada uno de estos puntos, 
se aplica un método de Newton-Raphson para encontrar la solución precisa dentro del dominio permitido. La función `test_function` asegura 
que los puntos resultantes estén dentro de la zona permitida.
"""
function CriticalEnergy(H, f!, h!, Q1min::Float64, Q1max::Float64, 
                        Pmin::Union{Float64, Irrational}, Pmax::Union{Float64, Irrational}, Q2lims, 
                        params, λ0::Float64, tol::Float64, ite::Int64, nrand::Int64)

    RandomPoints = Matrix{Float64}(undef, nrand, 5)

    s = Taylor1(Float64, 1)

    S = zeros(Taylor1{Float64}, 4, 5)

    @inbounds for i in 1:4
            S[i, i] = s
    end

    J = zeros(4, 4)

    dx = zeros(Taylor1{Float64}, 5)

    ΔQ1 = Q1max - Q1min
    ΔP = Pmax - Pmin

    ΔQ2 = 0.0
    Q2min = 0.0
    Q2max = 0.0

    for i in 1:nrand

        RandomPoints[i,1] = Q1min + rand()*ΔQ1

        RandomPoints[i,2] = Pmin + rand()*ΔP

        Q2min, Q2max = Q2lims(RandomPoints[i, :], params)
        ΔQ2 = Q2max - Q2min

        RandomPoints[i,3] = Q2min + rand()*ΔQ2

        RandomPoints[i,4] = Pmin + rand()*ΔP

        RandomPoints[i,5] = λ0

    end

    H_values = map(x -> H(x, params), eachrow(RandomPoints))

    imin = findmin(H_values)[2]

    xmin = RandomPoints[imin, :]

    # @show(xmin)

    imax = findmax(H_values)[2]

    xmax = RandomPoints[imax, :]

    # @show(xmax)

    maxnorm = _CriticalEnergy!(f!, h!, xmax, params, tol, ite, S, J, dx)
    
    minnorm = _CriticalEnergy!(f!, h!, xmin, params, tol, ite, S, J, dx)

    @show(maxnorm)
    @show(minnorm)
    
    return xmin, xmax

end
















function MinimalEnergy(H, f!, h!, Q1min::Float64, Q1max::Float64, 
                        Pmin::Union{Float64, Irrational}, Pmax::Union{Float64, Irrational}, Q2lims, 
                        params, λ0::Float64, tol::Float64, ite::Int64, nrand::Int64)

    RandomPoints = Matrix{Float64}(undef, nrand, 5)

    s = Taylor1(Float64, 1)

    S = zeros(Taylor1{Float64}, 4, 5)

    @inbounds for i in 1:4
        S[i, i] = s
    end

    J = zeros(4, 4)

    dx = zeros(Taylor1{Float64}, 5)

    ΔQ1 = Q1max - Q1min
    ΔP = Pmax - Pmin

    ΔQ2 = 0.0
    Q2min = 0.0
    Q2max = 0.0

    for i in 1:nrand

        RandomPoints[i,1] = Q1min + rand()*ΔQ1

        RandomPoints[i,2] = Pmin + rand()*ΔP

        Q2min, Q2max = Q2lims(RandomPoints[i, :], params)
        ΔQ2 = Q2max - Q2min

        RandomPoints[i,3] = Q2min + rand()*ΔQ2

        RandomPoints[i,4] = Pmin + rand()*ΔP

        RandomPoints[i,5] = λ0

    end

    H_values = map(x -> H(x, params), eachrow(RandomPoints))

    i = findmin(H_values)[2]

    x = RandomPoints[i, :]

    norm = _CriticalEnergy!(f!, h!, x, params, tol, ite, S, J, dx)

    if norm > tol
        @show(evaluate(dx))
        @warn("La energía no convergió en el randgo deseado")
    end

    return x

end