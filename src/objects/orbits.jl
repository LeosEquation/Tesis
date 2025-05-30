struct PeriodicBranch{U<:Real}
    t::LinRange{U}
    x::Array{U,3}    
    μ::Array{Complex{U},2}
    stbl::Array{Bool, 1}
end