using TaylorIntegration, TaylorSeries, LinearAlgebra

# Inicializando objetos
include("./objects/equilibria.jl")
include("./inicialization/taylornini.jl")

# Inicializando operaciones extra
include("./operations/adjugate.jl")
include("./operations/biproduct.jl")

# Equilibrio
include("./equilibria/eqfinding.jl")
include("./equilibria/eqsystem.jl")
include("./equilibria/eqcontinuation.jl")
include("./equilibria/branchswitching.jl")

# Limit Points
include("./limitpoint/lpfinding.jl")
include("./limitpoint/lpsystem.jl")
include("./limitpoint/lpcontinuation.jl")

# Hopf
include("./hopf/hpffinding.jl")
include("./hopf/hpfsystem.jl")
include("./hopf/hpfcontinuation.jl")

# Branch Points
include("./branchpoint/bpfinding.jl")
include("./branchpoint/bpsystem.jl")
include("./branchpoint/bpcontinuation.jl")

# Orbits 
include("./orbits/orbitfinding.jl")
include("./orbits/orbitsystem.jl")
include("./orbits/orbitcontinuation.jl")
include("./orbits/multishooting.jl")

# Limit Point Cycle
include("./limitpointcycle/lpcfinding.jl")
include("./limitpointcycle/lpcsystem.jl")
include("./limitpointcycle/lpccontinuation.jl")

# PeriodDoubling
include("./perioddoubling/pdfinding.jl")
include("./perioddoubling/pdsystem.jl")
include("./perioddoubling/pdcontinuation.jl")

# Poincare Section 
include("./poincaresection/poincaresection.jl")
include("./poincaresection/energylims.jl")

# Hamiltonian Orbits 
include("./hamiltonianorbits/orbitfinding.jl")
include("./hamiltonianorbits/orbitsystem.jl")
include("./hamiltonianorbits/orbitcontinuation.jl")