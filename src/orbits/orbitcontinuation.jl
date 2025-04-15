
function OrbitSystem!(F::Vector{Float64}, evalxT::Vector{Float64}, x1::Vector{Float64},
                    x0::Vector{Float64}, dx0::Vector{Float64},
                    Φ::Vector{Float64}, Δs::Float64, n::Int64)

    @inbounds F[1:n-2] .= evalxT .- x1[1:n-2]
    @inbounds F[n-1] = dot(dx0, x1 .- x0)
    @inbounds F[n] = dot(Φ, x1 .- x0) - Δs

end
function OrbitJacobian!(J::Matrix{Float64},
                        M::Matrix{Float64}, ϕ_p::Vector{Float64}, dx1::Vector{Float64}, 
                        dx0::Vector{Float64}, Φ::Vector{Float64},
                        n::Int64)

    @inbounds J[1:n-2, 1:n-2] .= M
    @inbounds for i in 1:n-2
        J[i, i] -= 1.0
    end
    @inbounds J[1:n-2, n-1] .= ϕ_p
    @inbounds J[1:n-2, n] .= dx1
    @inbounds J[n-1, 1:n-2] .= dx0[1:n-2]
    @inbounds J[n, :] .= Φ

end


function PDSystem!(F::Vector{Float64}, evalxT::Vector{Float64}, x1::Vector{Float64},
                   x0::Vector{Float64}, dx0::Vector{Float64},
                   PDtest::Float64, n::Int64)

@inbounds F[1:n-2] .= evalxT .- x1[1:n-2]
@inbounds F[n-1] = dot(dx0, x1 .- x0)
@inbounds F[n] = PDtest

end
function PDJacobian!(J::Matrix{Float64},
        M::Matrix{TaylorN{Float64}}, ϕ_p::Vector{Float64}, dx1::Vector{Float64}, 
        dx0::Vector{Float64}, I_n,
        n::Int64)

@inbounds J[1:n-2, 1:n-2] .= M
@inbounds for i in 1:n-2
J[i, i] -= 1.0
end
@inbounds J[1:n-2, n-1] .= ϕ_p
@inbounds J[1:n-2, n] .= dx1
@inbounds J[n-1, 1:n-2] .= dx0[1:n-2]
@inbounds for j in 1:n
    J[n, j] = tr( adjugate(evaluate(M) .- I_n) * evaluate(differentiate.(M, j)) )
end

end


function OrbitContinuation(f!, x_ini::Vector{Float64}, λ_ini::ComplexF64, v_ini::Vector{ComplexF64}, 
                           params, Δs::Float64, maxsteps::Int64, tol::Float64, ite::Int64,
                           integorder::Int64, integtol::Float64, integmaxsteps::Int64)

    ### Inicialización -----------------------------------------------------------------------------

    ### Dimensión

    n = length(x_ini)

    ### TaylorSeries

    δx = set_variables("δx", numvars = n, order = 1)
    δx2 =  set_variables("δx", numvars = n, order = 2)

    ### TaylorInteg

    xT = Array{TaylorN{Float64},1}(undef, n)
    xT2 = Array{TaylorN{Float64},1}(undef, n)
    t_integ = 0.0 + Taylor1( Float64 , integorder )
    x_integ = Array{Taylor1{TaylorN{Float64}}}(undef, n)
    dx_integ = Array{Taylor1{TaylorN{Float64}}}(undef, n)
    xaux_integ = Array{Taylor1{TaylorN{Float64}}}(undef, n)

    @inbounds for i in eachindex(x_ini)
        xT[i] = x_ini[i] + δx[i]
        xT2[i] = x_ini[i] + δx2[i]
        x_integ[i] = Taylor1( xT[i], integorder )
        dx_integ[i] = Taylor1( zero(xT[i]), integorder )
    end

    ### Almacenamiento

    x = Array{Float64,2}(undef, maxsteps, n)
    freal = Array{Float64,2}(undef, maxsteps, n-2)
    fimag = Array{Float64,2}(undef, maxsteps, n-2)
    PDtest = Array{Float64,1}(undef, maxsteps)
    LPCtest = Array{Float64,1}(undef, maxsteps)
    M = Array{Float64, 2}(undef, n-2, n-2)
    Mb = Array{TaylorN{Float64}, 2}(undef, n-2, n-2)
    evalxT = Array{Float64, 1}(undef, n-2)
    ϕ_p = Array{Float64, 1}(undef, n-2)
    dx1 = Array{Float64, 1}(undef, n-2)
    I_n = LinearAlgebra.I(n-2)

    μ = Array{ComplexF64, 1}(undef, n-2)
    pd = Vector{Float64}[] # Bifurcaciones de rama
    pddir  = Vector{Float64}[] # Dirección tangente de la rama en la bifurcación de ram
    pdeval  = Vector{ComplexF64}[]
    pdevec  = Matrix{ComplexF64}[]
    
    ### Continuación
    
    J = zeros(Float64, n, n)
    F = zeros(Float64, n)
    Φ = zeros(Float64, n)
    v = zeros(Float64, n)
    v[n] = 1.0

    ### Pasos iniciales

    x0 = copy(x_ini)
    dx0 = zero(x0)
    x1 = copy(x_ini)
    PDtest0 = 0.0
    PDtest1 = 0.0
    PDtest2 = 0.0
    LPCtest0 = 0.0
    LPCtest1 = 0.0
    x2 = copy(x0)

    ### Primer almacenamiento

    x[1, :] .= x0
    dx0[1:n-2] .= abs(imag(λ_ini)) * imag(v_ini)
    Φ[1:n-2] .= real(v_ini)

    taylorinteg_optim!(f!, xT, 0.0, 1.0, integorder, integtol,
                       t_integ, x_integ, dx_integ, xaux_integ, paramsP; maxsteps = integmaxsteps)

    @inbounds for i in 1:n-2
        dx1[i] = evaluate(differentiate(xT[i], n))
        evalxT[i] = evaluate(xT[i])
    end

    @inbounds for j in 1:n-2
        ϕ_p[j] = evaluate(differentiate(xT[j], n-1))
        for i in 1:n-2
            M[i, j] = evaluate(differentiate(xT[i], j))
        end
    end

    μ .= eigvals(M)

    freal[1, :] .= real.(μ)
    fimag[1, :] .= imag.(μ)
    PDtest0 = det(M + I_n)
    LPCtest0 = 0.0

    PDtest[1]  = PDtest0

    LPCtest[1] = LPCtest0

    ### Continuación ------------------------------------------------------------------------------------

    i = 2

    while i <= maxsteps

        @. x1 = x0 + (Δs * Φ)

        xT .= x1 .+ δx

        taylorinteg_optim!(f!, xT, 0.0, 1.0, integorder, integtol,
                           t_integ, x_integ, dx_integ, xaux_integ, paramsP; maxsteps = integmaxsteps)

        @inbounds for i in 1:n-2
            dx1[i] = evaluate(differentiate(xT[i], n))
            evalxT[i] = evaluate(xT[i])
        end
    
        @inbounds for j in 1:n-2
            ϕ_p[j] = evaluate(differentiate(xT[j], n-1))
            for i in 1:n-2
                M[i, j] = evaluate(differentiate(xT[i], j))
            end
        end

        OrbitJacobian!(J, M, ϕ_p, dx1, dx0, Φ, n)
        OrbitSystem!(F, evalxT, x1, x0, dx0, Φ, Δs, n)

        j = 1

        while j <= ite && norm(F) > tol

            if abs(det(J)) == 0.0
                break
            end

            x1 .-= J \ F

            xT .= x1 .+ δx

            taylorinteg_optim!(f!, xT, 0.0, 1.0, integorder, integtol,
                            t_integ, x_integ, dx_integ, xaux_integ, paramsP; maxsteps = integmaxsteps)

            @inbounds for i in 1:n-2
                dx1[i] = evaluate(differentiate(xT[i], n))
                evalxT[i] = evaluate(xT[i])
            end
        
            @inbounds for j in 1:n-2
                ϕ_p[j] = evaluate(differentiate(xT[j], n-1))
                for i in 1:n-2
                    M[i, j] = evaluate(differentiate(xT[i], j))
                end
            end
    
            OrbitJacobian!(J, M, ϕ_p, dx1, dx0, Φ, n)
            OrbitSystem!(F, evalxT, x1, x0, dx0, Φ, Δs, n)

            j += 1

        end

        if norm(F) >= tol
            @warn("La norma del sistema no convergió. Abortando.")
            break
        end

        x[i, :] .= x1

        Φ .= J \ v

        normalize!(Φ)

        μ .= eigvals(M)

        freal[i, :] .= real.(μ)
        fimag[i, :] .= imag.(μ)
        PDtest0 = det(M + I_n)
        LPCtest0 = 0.0

        PDtest[i]  = PDtest0

        LPCtest[i] = LPCtest0

        x0 .= x1

        f!(dx0, x0, params, 0.0)

        print("\r Progreso: $(round(i * 100 / maxsteps, digits = 2)) % \t")

        if false # sign(PDtest0) != sign(PDtest1) && sign(LPCtest0) == sign(LPCtest1)

            x2 .= 0.5 * (x1 + x0)

            xT2 .= x2 .+ δx2

            taylorinteg_optim!(f!, xT2, 0.0, 1.0, integorder, integtol,
                               t_integ, x_integ, dx_integ, xaux_integ, paramsP; maxsteps = integmaxsteps)
    
            @inbounds for i in 1:n-2
                dx1[i] = evaluate(differentiate(xT2[i], n))
                evalxT[i] = evaluate(xT2[i])
            end
        
            @inbounds for j in 1:n-2
                ϕ_p[j] = evaluate(differentiate(xT2[j], n-1))
                for i in 1:n-2
                    Mb[i, j] = differentiate(xT[i], j)
                end
            end

            PDtest2 = det(evaluate(Mb) + I_n)
    
            PDJacobian!(J, Mb, ϕ_p, dx1, dx0, I_n, n)
            PDSystem!(F, evalxT, x1, x0, dx0, PDtest2, Δs, n)

            j = 1

            while j <= ite && norm(F) > tol

                if abs(det(J)) == 0.0
                    break
                end

                x2 .-= J \ F

                xT2 .= x2 .+ δx2

                taylorinteg_optim!(f!, xT2, 0.0, 1.0, integorder, integtol,
                               t_integ, x_integ, dx_integ, xaux_integ, paramsP; maxsteps = integmaxsteps)
    
                @inbounds for i in 1:n-2
                    dx1[i] = evaluate(differentiate(xT2[i], n))
                    evalxT[i] = evaluate(xT2[i])
                end
            
                @inbounds for j in 1:n-2
                    ϕ_p[j] = evaluate(differentiate(xT2[j], n-1))
                    for i in 1:n-2
                        Mb[i, j] = differentiate(xT[i], j)
                    end
                end

                PDtest2 = det(evaluate(Mb) + I_n)
        
                PDJacobian!(J, Mb, ϕ_p, dx1, dx0, I_n, n)
                PDSystem!(F, evalxT, x1, x0, dx0, PDtest2, Δs, n)

                j += 1

            end

            J[n, :] .= Φ

            Φ2 .= J \ v
            normalize!(Φ2)

            evec .= eigvecs(evaluate(Fx2))
            eval .= eigvals(evaluate(Fx2))
            push!(h, copy(x2))
            push!(hdir, copy(Φ2))
            push!(hevec, copy(evec))
            push!(heval, copy(eval))


            if norm(F) > tol
                @warn("det(2Fₓ(x,λ)⊛Iₙ) = $(F[end]). La bifurcación de Hopf $(length(h)) no convergió.")
            end

        end

        i += 1

    end

    return x[1:i-1, :], PDtest[1:i-1], LPCtest[1:i-1], freal[1:i-1, :], fimag[1:i-1, :]

end

#=

function OrbitContinuation(f!, x_ini::Vector{Float64}, 
                            params, Δs::Float64, maxsteps::Int64, tol::Float64, ite::Int64,
                            integorder::Int64, integtol::Float64, integmaxsteps::Int64)

### Inicialización -----------------------------------------------------------------------------

### Dimensión

n = length(x_ini)

### TaylorSeries

δx = set_variables("δx", numvars = n, order = 1)
δx_zero = zero(x_ini)

### TaylorInteg

xT = Array{TaylorN{Float64},1}(undef, n)
t_integ = 0.0 + Taylor1( Float64 , integorder )
x_integ = Array{Taylor1{TaylorN{Float64}}}(undef, n)
dx_integ = Array{Taylor1{TaylorN{Float64}}}(undef, n)
xaux_integ = Array{Taylor1{TaylorN{Float64}}}(undef, n)

@inbounds for i in eachindex(x_ini)
xT[i] = x_ini[i] + δx[i]
x_integ[i] = Taylor1( xT[i], integorder )
dx_integ[i] = Taylor1( zero(xT[i]), integorder )
end

### Almacenamiento

x = Array{Float64,2}(undef, maxsteps, n)
PDtest = Array{Float64,1}(undef, maxsteps)
LPCtest = Array{Float64,1}(undef, maxsteps)
JxT = Array{Float64, 2}(undef, n, n)
evalxT = Array{Float64, 1}(undef, n)
# dx = deepcopy(δx)

### Continuación

J = zeros(Float64, n, n)
F = zeros(Float64, n)
Φ = zeros(Float64, n)
v = zeros(Float64, n)
v[n] = 1.0

### Matrix Identidad

I_n = LinearAlgebra.I(n)
I_n2 = LinearAlgebra.I(n-2)
M = zeros(Float64, n-2, n-2)

### Pasos iniciales

x0 = copy(x_ini)
dx0 = zero(x0)
x1 = copy(x_ini)

x[1, :] .= x0

f!(dx0, x0, params, 0.0)
# dx0[1:n-2] .= abs(imag(λ_ini)) * imag(v_ini)
# Φ[1:n-2] .= real(v_ini)

xT .= x1 .+ δx

taylorinteg_optim!(f!, xT, 0.0, 1.0, integorder, integtol,
                   t_integ, x_integ, dx_integ, xaux_integ, paramsP; maxsteps = integmaxsteps)
TaylorSeries.evaluate!(xT, δx_zero, evalxT)
TaylorSeries.jacobian!(JxT, xT)
JxT .-= I_n
OrbitJacobian!(J, JxT, dx0, Φ, n)
OrbitSystem!(F, evalxT, x1, x0, dx0, Φ, Δs, n)

NS = nullspace(J)

if size(NS, 2) == 0
    @show(J)
    error("No hay dirección a seguir")
end

if size(NS, 2) > 1
    @show(NS)
    error("Hay más de una dirección inicial a seguir")
end

Φ .= NS

### Continuación ------------------------------------------------------------------------------------

i = 2

while i <= maxsteps

@. x1 = x0 + (Δs * Φ)

xT .= x1 .+ δx

taylorinteg_optim!(f!, xT, 0.0, 1.0, integorder, integtol,
    t_integ, x_integ, dx_integ, xaux_integ, paramsP; maxsteps = integmaxsteps)
TaylorSeries.evaluate!(xT, δx_zero, evalxT)
TaylorSeries.jacobian!(JxT, xT)
JxT .-= I_n
OrbitJacobian!(J, JxT, dx0, Φ, n)
OrbitSystem!(F, evalxT, x1, x0, dx0, Φ, Δs, n)

j = 1

while j <= ite && norm(F) > tol

if abs(det(J)) == 0.0
break
end

x1 .-= J \ F

xT .= x1 .+ δx

taylorinteg_optim!(f!, xT, 0.0, 1.0, integorder, integtol,
     t_integ, x_integ, dx_integ, xaux_integ, paramsP; maxsteps = integmaxsteps)
TaylorSeries.evaluate!(xT, δx_zero, evalxT)
TaylorSeries.jacobian!(JxT, xT)
JxT .-= I_n
OrbitJacobian!(J, JxT, dx0, Φ, n)
OrbitSystem!(F, evalxT, x1, x0, dx0, Φ, Δs, n)

j += 1

end

if norm(F) >= tol
@warn("La norma del sistema no convergió. Abortando.")
break
end

M .= JxT[1:n-2, 1:n-2]

PDtest[i] = det(M + I_n2)
LPCtest[i] = det(M - I_n2)

x[i, :] .= x1

Φ .= J \ v

normalize!(Φ)

x0 .= x1

f!(dx0, x0, params, 0.0)

print("\r i = $i : norm = $(norm(F)) \t")

i += 1

end

return x[1:i-1, :], PDtest[1:i-1], LPCtest[1:i-1]

end

=#