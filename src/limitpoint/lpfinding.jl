

# x1 .= x0 .+ (Φ * Δs)

#         f!(dx, x1[1:n+2] + δx, params, 0.0)
        
#         for j in 1:n
#             Fx[:, j] .= differentiate.(dx[1:n], j)
#         end

#         LPJacobian!(J, Fx, dx, x0, x1, Φ, n)
#         LPSystem!(F, Fx, dx, x0, x1, Φ, Δs, n)

#         j = 1

#         while j <= ite && norm(F) > tol

#             # println("$j : norm = $(norm(F))")

#             if abs(det(J)) == 0.0
#                 break
#             end

#             x1 .-= J \ F
    
#             f!(dx, x1[1:n+2] + δx, params, 0.0)
#             for j in 1:n
#                 Fx[:, j] .= differentiate.(dx[1:n], j)
#             end

#             LPJacobian!(J, Fx, dx, x0, x1, Φ, n)
#             LPSystem!(F, Fx, dx, x0, x1, Φ, Δs, n)

#             j += 1

#         end

        # if norm(F) >= tol
        #     @warn("La norma del sistema no convergió. Abortando.")
        #     break
        # end