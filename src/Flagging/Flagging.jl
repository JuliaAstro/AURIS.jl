# To be replaced by RFI module
module Flagging

using Statistics

using ..Data
using ..Data: VisibilityDataset

include("sumthreshold.jl")

export flag_rfi!, flag_edges!, sumthreshold!

end
