g_bt(evalA, adjA, B, C) = - det(evalA) / dot(C, (adjA * B))

function ∇g_bt!(∇, A, adjA, B, C, n)

    for j in 1:n 

        ∇[j] = - dot(C, (((adjA * B) * evaluate(differentiate.(A, j))) * adjA)) / dot(C, (adjA * B))^2

    end
    
end