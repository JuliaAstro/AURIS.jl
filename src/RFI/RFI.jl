module RFI

include("Bezier/BBasis.jl")
include("Bezier/BProbability.jl")
include("Bezier/BFitting.jl")
include("Bezier/BAPI.jl")

include("Splines/SBasis.jl")
include("Splines/SFitting.jl")
include("Splines/SOptimisation.jl")
include("Splines/SAPI.jl")

using .BBasis
using .BProbability
using .BFitting
using .BAPI

using .SBasis
using .SFitting
using .SOptimisation
using .SAPI

export BAPI
export SAPI

end