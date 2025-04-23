

function BPContinuation(f!, x_ini::Vector{Float64}, params, Δs::Float64, maxsteps::Int64, tol::Float64, ite::Int64)

#-

    n = length(x_ini) # Dimensión del sistema

    x = Array{Float64,2}(undef, maxsteps, n) # Rama de equilibrio

###

# Inicializando TaylorN

    variable_names = [string("δx", TaylorSeries.subscriptify(i)) for i in 1:n]

    TaylorSeries.set_variables(U, variable_names, order = 2)

    δx = TaylorN.(1:n, order = 2)

    xaux = copy(δx)
    dx = zero(δx)

    q0 = zeros(U, 2*n + 1)
    q0[1:n] .= x_ini
    q0[n+1:2*n-2] .= v_ini

    q1 = zeros(U, 2*n + 1)
    q1[1:n] .= x_ini
    q1[n+1:2*n-2] .= v_ini

    ###

    dx = Array{TaylorN{U},1}(undef, n)

    ###

    # Inicializando variables para la continuación

    J = zeros(U, 2*n + 1, 2*n - 2)
    F = zeros(U, 2*n - 2)
    Φ = zeros(U, 2*n - 2)
    v = zeros(U, 2*n - 2)
    v[2*n - 2] = one(U)

    Jeval = zeros(U, n, n)
    dxeval = zeros(U, n)
    xzero = zeros(U, n) 

###

# Iniciamos el primer paso

x[1, :] .= x1#[1:n]

f!(dx, x1[1:n] + δx, params, 0.0)

for j in 1:n-1
    Fx[:, j] .= differentiate.(dx[1:n-2], j)
end

@assert rank(evaluate(Fx)) == n-3 "La bifurcación no es un Branch Point (BP)"

x1[n+1:2*n-2] .= nullspace(evaluate(Fx)')

BPJacobian!(J, Fx, dx, x1, Φ, n)
BPSystem!(F, Fx, dx, x0, x1, Φ, Δs, n)

# @show F

# return J

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

for j in 1:n-1
    Fx[:, j] .= differentiate.(dx[1:n-2], j)
end

BPJacobian!(J, Fx, dx, x1, Φ, n)
BPSystem!(F, Fx, dx, x0, x1, Φ, Δs, n)

# return F, J

j = 1

while j <= ite && norm(F) > tol

# println("j = $j : norm = $(norm(F))")

if abs(det(J)) == 0.0
break
end

x1 .-= J \ F

f!(dx, x1[1:n] + δx, params, 0.0)

for j in 1:n-1
    Fx[:, j] .= differentiate.(dx[1:n-2], j)
end

BPJacobian!(J, Fx, dx, x1, Φ, n)
BPSystem!(F, Fx, dx, x0, x1, Φ, Δs, n)


j += 1

end

if norm(F) >= tol
@warn("La norma del sistema no convergió. Abortando.")
break
end

x[i, :] .= x1#[1:n]

Φ .= J \ v

normalize!(Φ)

print("\r i = $i : norm = $(norm(F)) \t")

x0 .= x1

i += 1

end

return x[1:i-1, 1:n]

end