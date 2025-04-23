struct LimitPoint{U<:Real}
    x::Array{U, 1}
    vector::Array{U, 1}
    value::U
end

struct HopfPoint{U<:Real}
    x::Array{U, 1}
    vector::Array{Complex{U}, 1}
    value::Complex{U}
end

struct BranchPoint{U<:Real}
    x::Array{U, 1}
    vector::Array{U, 1}
    value::U
    dir::Array{U, 1}
end

struct EquilibriumBranch{U<:Real}
    x::Array{U,2}    
    realeigen::Array{U,2}
    imageigen::Array{U,2}
    lp::Array{LimitPoint{U}, 1}
    hp::Array{HopfPoint{U}, 1}
    bp::Array{BranchPoint{U}, 1}
    stbl::Array{Bool, 1}
end

function EquilibriumBranch(x::Array{U,2}, realev::Array{U,2}, imagev::Array{U,2},
                           lp::Array{LimitPoint{U}, 1}, hp::Array{HopfPoint{U}, 1}, bp::Array{BranchPoint{U}, 1},
                           stbl::Array{Bool,1}, i::T) where {U<:Real, T<:Integer}

    return EquilibriumBranch(x[1:i-1, :], realev[1:i-1, :], imagev[1:i-1, :],
                             lp, hp, bp, stbl[1:i-1])
                             
end

function Base.show(io::IO, obj::EquilibriumBranch)
    println(io, " Equilibrium Branch")
    println(io, " # steps : $(size(obj.x, 1))")
    println(io, " # Limit Point Bifurcations : $(length(obj.lp))")
    println(io, " # Branch Point Bifurcations : $(length(obj.bp))")
    println(io, " # Hopf Point Bifurcations : $(length(obj.hp))")
end