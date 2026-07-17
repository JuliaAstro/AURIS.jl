"""
    Data
  * `uvw.jl`     – TAI→UTC time handling and ITRF→J2000 UVW geometry
                   (SatelliteToolboxTransformations for Earth orientation).
  * `dataset.jl` – `VisibilityDataset`, intent/field/scan selection helpers,
                   and `load_visibility_dataset`.
  * `statwt.jl`  – data-driven inverse-variance weights (`statwt!`).
  * `tsysweights.jl` – switched-power Tsys weights (`tsys_weights!`).
"""
module Data

using LinearAlgebra
using Statistics
using StaticArrays
using SatelliteToolboxTransformations

using ..Components
using ..Components: Main, Flag
using ..Loader: SDM

include("uvw.jl")
include("dataset.jl")
include("statwt.jl")
include("tsysweights.jl")

export VisibilityDataset, load_visibility_dataset, statwt!, tsys_weights!
export scans_with_intent, mains_for_scans, mains_for_field, field_by_id
export baseline_pairs, uvw!, uvw_basis, itrf_to_j2000, tai_ns_to_jd_utc
export n_time, n_baseline, n_spw, n_chan, n_pol

end 
