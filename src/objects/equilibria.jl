struct LimitPoint{U<:Real}
    x::Array{U, 1}
    vector::Array{Complex{U}, 1}
    value::Complex{U}
end

function LimitPoint(x::Vector{T}, v::Vector{T}, p::T) where T<:Real
    LimitPoint(x, complex.(v), complex(p))
end

struct HopfPoint{U<:Real}
    x::Array{U, 1}
    vector::Array{Complex{U}, 1}
    value::Complex{U}
end

function HopfPoint(x::Vector{T}, v::Vector{T}, p::T) where T<:Real
    HopfPoint(x, complex.(v), complex(p))
end

struct BranchPoint{U<:Real}
    x::Array{U, 1}
    vector::Array{Complex{U}, 1}
    value::Complex{U}
    dir::Array{U, 1}
end

function BranchPoint(x::Vector{T}, v::Vector{T}, p::T, dir::Vector{T}) where T<:Real
    BranchPoint(x, complex.(v), complex(p), dir)
end

struct EquilibriumBranch{U<:Real}
    x::Array{U,2}    
    λ::Array{Complex{U},2}
    lp::Array{LimitPoint{U}, 1}
    hp::Array{HopfPoint{U}, 1}
    bp::Array{BranchPoint{U}, 1}
    stbl::Array{Bool, 1}
end

function Base.show(io::IO, obj::EquilibriumBranch)
    println(io, " Equilibrium Branch")
    println(io, " # steps : $(size(obj.x, 1))")
    println(io, " # Limit Point Bifurcations : $(length(obj.lp))")
    println(io, " # Branch Point Bifurcations : $(length(obj.bp))")
    println(io, " # Hopf Point Bifurcations : $(length(obj.hp))")
end