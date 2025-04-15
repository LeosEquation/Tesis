
function closest_eigenpair_by_real_part(F::Eigen)
    idx = argmin(abs.(real.(F.values)))
    return F.values[idx], F.vectors[:, idx]
end