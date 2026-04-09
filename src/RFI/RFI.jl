module RFI

include("Bezier/BBasis.jl")
include("Bezier/BProbability.jl")
include("Bezier/BFitting.jl")
include("Bezier/BAPI.jl")

include("Splines/SBasis.jl")
include("Splines/SProbability.jl")
include("Splines/SFitting.jl")
include("Splines/SOptimisation.jl")
include("Splines/SProjection.jl")
include("Splines/SAPI.jl")

using .BBasis
using .BProbability
using .BFitting
using .BAPI

using .SBasis
using .SProbability
using .SFitting
using .SOptimisation
using .SProjection
using .SAPI

export BAPI
export SAPI

end