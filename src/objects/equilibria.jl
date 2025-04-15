struct EquilibriumBranch{U<:Real}
    x::Array{U,2}    
    realeval::Array{U,2}
    imageval::Array{U,2}
    lp::Array{U,2}
    lpeval::Array{Complex{U},2}
    lpevec::Array{Complex{U},3}
    bp::Array{U,2}
    bpdir::Array{U,2}
    bpeval::Array{Complex{U},2}
    bpevec::Array{Complex{U},3}
    h::Array{U,2}
    heval::Array{Complex{U},2}
    hevec::Array{Complex{U},3}
    stbl::Array{Bool, 1}
end

function EquilibriumBranch(x::Array{U,2}, realev::Array{U,2}, imagev::Array{U,2},
                           lp::Array{Array{U, 1}, 1},
                           lpeval::Array{Array{Complex{U}, 1}, 1}, lpevec::Array{Array{Complex{U}, 2}, 1},
                           bp::Array{Array{U, 1}, 1}, bpdir::Array{Array{U, 1}, 1}, 
                           bpeval::Array{Array{Complex{U}, 1}, 1}, bpevec::Array{Array{Complex{U}, 2}, 1},
                           h::Array{Array{U, 1}, 1},
                           heval::Array{Array{Complex{U}, 1}, 1}, hevec::Array{Array{Complex{U}, 2}, 1},
                           stbl::Array{Bool,1}, n::T, i::T) where {U<:Real, T<:Integer}

    lp_ = [i[j] for i in lp, j in 1:n]
    lpeval_ = [i[j] for i in lpeval, j in 1:n-1]
    lpevec_ = [i[j, k] for i in lpevec, j in 1:n-1, k in 1:n-1]

    bp_ = [i[j] for i in bp, j in 1:n]
    bpdir_ = [i[j] for i in bpdir, j in 1:n]
    bpeval_ = [i[j] for i in bpeval, j in 1:n-1]
    bpevec_ = [i[j, k] for i in bpevec, j in 1:n-1, k in 1:n-1]

    h_ = [i[j] for i in h, j in 1:n]
    heval_ = [i[j] for i in heval, j in 1:n-1]
    hevec_ = [i[j, k] for i in hevec, j in 1:n-1, k in 1:n-1]

    return EquilibriumBranch(x[1:i-1, :], realev[1:i-1, :], imagev[1:i-1, :],
                             lp_,
                             lpeval_, lpevec_,
                             bp_, bpdir_, bpeval_, bpevec_,
                             h_,
                             heval_, hevec_, stbl[1:i-1])
                             
end

function Base.show(io::IO, obj::EquilibriumBranch)
    println(io, " Equilibrium Branch")
    println(io, "  → # Limit Point Bifurcations : $(size(obj.lp, 1))")
    println(io, "  → # Branch Point Bifurcations : $(size(obj.bp, 1))")
    println(io, "  → # Hopf Point Bifurcations : $(size(obj.h, 1))")
end