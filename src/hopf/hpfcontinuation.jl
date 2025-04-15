

function HopfJacobian!(J::Matrix{Float64}, Fx::Matrix{TaylorN{Float64}}, w::Vector{Float64}, 
                       dx::Vector{TaylorN{Float64}}, x1::Vector{Float64}, Φ::Vector{Float64}, n::Int64)

    J[1:n, 1:n] .= TaylorSeries.jacobian(dx)

    for j in 1:n
        J[n-1:2*n-4, j] .= evaluate(differentiate.(Fx ^ 2, j)) * x1[n+1:2*n-2]
    end

    J[n-1:2*n-4, n+1:2*n-2] .= evaluate(Fx) ^ 2

    for j in 1:n-2
        J[n-2+j, n+j] += x1[2*n-1]
    end

    J[n-1:2*n-4, 2*n-1] .= x1[n+1:2*n-2]

    J[2*n-3, n+1:2*n-2] .= 2.0 * x1[n+1:2*n-2]

    J[2*n-2, n+1:2*n-2] .= w

    J[2*n-1, :] .= Φ

end

function HopfSystem!(F::Vector{Float64}, Fx::Matrix{TaylorN{Float64}}, dx::Vector{TaylorN{Float64}},
                   x0::Vector{Float64}, x1::Vector{Float64}, w::Vector{Float64},
                   Φ::Vector{Float64}, Δs::Float64, n::Int64)

    F[1:n-2] .= evaluate(dx[1:n-2])
    F[(n-1):(2*n-4)] .= evaluate(Fx) ^ 2 * x1[n+1:2*n-2] + x1[2*n-1] * x1[n+1:2*n-2]
    F[2*n-3] = dot(x1[n+1:2*n-2], x1[n+1:2*n-2]) - 1.0
    F[2*n-2] = dot(w, x1[n+1:2*n-2])
    F[2*n-1] = dot(x1 - x0, Φ) - Δs

end


function HopfContinuation(f!, x_ini::Vector{Float64}, v_ini::Vector{Float64}, ω_ini::Float64, params, Δs::Float64, maxsteps::Int64, tol::Float64, ite::Int64)

#-

n = length(x_ini) # Dimensión del sistema

x = Array{Float64,2}(undef, maxsteps, 2*n - 1) # Rama de equilibrio

###

# Inicializando TaylorN

δx = set_variables("δx", numvars = n, order = 2)

###

# Inicializando valores previos

x0 = zeros(Float64, 2*n - 1)
x0[1:n] .= copy(x_ini)
x0[n+1:2*n-2] .= copy(v_ini) / norm(v_ini)
x0[2*n-1] = ω_ini ^ 2

###

# Inicializando valores posteriores

x1 = copy(x0)

###

# Inicializando variables para bifurcaciones

dx = Array{TaylorN{Float64},1}(undef, n)
Fx = Array{TaylorN{Float64},2}(undef, n-2, n-2)
w = zeros(Float64, n-2)

###

# Inicializando variables para la continuación

J = zeros(Float64, 2*n - 1, 2*n - 1)
F = zeros(Float64, 2*n - 1)
Φ = zeros(Float64, 2*n - 1)
v = zeros(Float64, 2*n - 1)
v[2*n - 1] = 1.0

###

# Iniciamos el primer paso

x[1, :] .= x1

f!(dx, x1[1:n] + δx, params, 0.0)

for j in 1:n-2
    Fx[:, j] .= differentiate.(dx[1:n-2], j)
end

w .= nullspace([(evaluate(Fx)^2 + x1[2*n-1] * LinearAlgebra.I(n-2))' ; v_ini'], atol = 1.e-10)

HopfJacobian!(J, Fx, w, dx, x1, Φ, n)
HopfSystem!(F, Fx, dx, x0, x1, w, Φ, Δs, n)

NS = nullspace(J)

if size(NS, 2) == 0
@show(J)
error("No hay dirección a seguir")
end

if size(NS, 2) > 1
display(NS)
error("Hay más de una dirección inicial a seguir")
end

x[1, :] .= x1

Φ .= NS

#

i = 2

while i <= maxsteps
# println(" i = $i :")

x1 .= x0 .+ (Φ * Δs)

f!(dx, x1[1:n] + δx, params, 0.0)

for j in 1:n-2
    Fx[:, j] .= differentiate.(dx[1:n-2], j)
end

HopfJacobian!(J, Fx, w, dx, x1, Φ, n)
HopfSystem!(F, Fx, dx, x0, x1, w, Φ, Δs, n)

# return F, J

j = 1

while j <= ite && norm(F) > tol

# println("j = $j : norm = $(norm(F))")

if abs(det(J)) == 0.0
break
end

x1 .-= J \ F

f!(dx, x1[1:n] + δx, params, 0.0)

for j in 1:n-2
    Fx[:, j] .= differentiate.(dx[1:n-2], j)
end

HopfJacobian!(J, Fx, w, dx, x1, Φ, n)
HopfSystem!(F, Fx, dx, x0, x1, w, Φ, Δs, n)

j += 1

end

if norm(F) >= tol
@warn("La norma del sistema no convergió. Abortando.")
break
end

if x1[2*n-1] <= tol
    @warn("La parte imaginaria se anuló. Abortando")
    break
    end

x[i, :] .= x1

w .= nullspace([(evaluate(Fx)^2 + x1[2*n-1] * LinearAlgebra.I(n-2))' ; x1[n+1:2*n-2]'], atol = 1.e-10)

Φ .= J \ v

normalize!(Φ)

print("\r i = $i : norm = $(norm(F)) \t")

x0 .= x1

i += 1

end

return x[1:i-1, 1:n]

end