"""
    Calibration

  * `stefcal.jl`   – the core alternating least-squares solver.
  * `gaintable.jl` – `GainTable` solutions and `applycal!`.
  * `solve.jl`     – `solve_bandpass` (`:B`, per channel) and `solve_gains`
                     (`:G`, per solution interval), point-source model
                     (scalar or frequency-dependent).
  * `fluxscale.jl` – Perley–Butler 2017 calibrator spectra (`setjy_flux`,
                     `setjy_model`) and secondary-calibrator flux
                     bootstrapping (`fluxscale`).
"""
module Calibration

using Statistics

using ..Data
using ..Data: VisibilityDataset

include("stefcal.jl")
include("gaintable.jl")
include("solve.jl")
include("fluxscale.jl")

export GainTable, stefcal, applycal!, merge_spw
export solve_bandpass, solve_gains, solve_selfcal
export setjy_flux, setjy_model, fluxscale

end 
