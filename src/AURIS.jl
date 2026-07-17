module AURIS

include("Components.jl")
using .Components

include("Loader/Loader.jl")
using .Loader

include("Data/Data.jl")
using .Data

include("Flagging/Flagging.jl")
using .Flagging

include("Calibration/Calibration.jl")
using .Calibration

include("Imaging/Imaging.jl")
using .Imaging

include("Pipeline/Pipeline.jl")
using .Pipeline

export Components, Loader, SDM, Data, Flagging, Calibration, Imaging, Pipeline
export VisibilityDataset, load_visibility_dataset, statwt!, tsys_weights!
export flag_rfi!, flag_edges!, sumthreshold!
export GainTable, stefcal, applycal!, merge_spw
export solve_bandpass, solve_gains, solve_selfcal
export setjy_flux, setjy_model, fluxscale
export uv_samples, uvw_samples, dirty_image, dirty_image_dft, nyquist_cell, briggs_weights
export dirty_beam, hogbom_clean, fit_beam, restore
export predict_vis, predict_dataset_vis, cs_clean
export dirty_image_w, dirty_image_wdft, predict_vis_w, cs_clean_w
export write_fits
export rank_refants, run_pipeline, pipeline_report, export_gains, export_summary
export PipelineResult, TargetImage

end