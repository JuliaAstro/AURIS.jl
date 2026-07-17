"""
    Pipeline
  * `refant.jl`   – `rank_refants` automatic reference-antenna ranking
                    (unflagged fraction × array-centre proximity).
  * `pipeline.jl` – `run_pipeline` stage driver → `PipelineResult`
                    (cal tables, per-target `TargetImage`s, `StageLog`s).
  * `report.jl`   – `pipeline_report` markdown run summary and
                    `export_gains` CSV export for CASA cross-validation.
"""
module Pipeline

using Dates
using LinearAlgebra: norm
using Statistics

using ..Components
import ..Components: Main as CMain
using ..Loader: SDM
using ..Data
using ..Data: VisibilityDataset, n_baseline
using ..Flagging
using ..Calibration
using ..Calibration: GainTable
using ..Imaging

include("refant.jl")
include("pipeline.jl")
include("report.jl")

export rank_refants, run_pipeline, pipeline_report, export_gains, export_summary
export PipelineResult, TargetImage, StageLog

end
