"""
    Imaging
  * `gridding.jl` – `uv_samples`, `dirty_image` (NFFT), `dirty_image_dft`
                    (direct-sum reference used to validate the fast path).
  * `clean.jl`    – Högbom CLEAN, PSF main-lobe beam fit, restoration.
  * `csclean.jl`  – `predict_vis` (forward NFFT degridding) and `cs_clean`
                    (Cotton–Schwab major/minor cycles).
  * `weighting.jl`– `briggs_weights` (Briggs 1995 robust density weighting).
  * `wstack.jl`   – w-stacking: `uvw_samples`-driven wide-field imaging
                    (`dirty_image_w`, `predict_vis_w`, `cs_clean_w`).
  * `fits.jl`     – `write_fits` standard radio FITS output (FITSFiles.jl).
"""
module Imaging

using NFFT
using PrecompileTools: @compile_workload
using FITSFiles
using FITSFiles: Card, HDU

using ..Data
using ..Data: VisibilityDataset

include("gridding.jl")
include("weighting.jl")
include("clean.jl")
include("csclean.jl")
include("wstack.jl")
include("fits.jl")

export uv_samples, uvw_samples, dirty_image, dirty_image_dft, nyquist_cell
export briggs_weights
export dirty_beam, hogbom_clean, fit_beam, restore
export predict_vis, predict_dataset_vis, cs_clean
export dirty_image_w, dirty_image_wdft, predict_vis_w, cs_clean_w
export write_fits

# Loading this package's dependency stack invalidates NFFT's precompiled
# convolution kernels, so they recompile at first use. This helps.
@compile_workload begin
    let n = 48
        θ  = [2π * k / n for k in 1:n]
        r  = 40 .* θ ./ 2π
        uv = permutedims([r .* cos.(3θ) r .* sin.(3θ)])
        old = NFFT._use_threads[]
        NFFT._use_threads[] = false
        try
            cs_clean(uv, ones(ComplexF64, n), ones(n);
                     npix = 32, cell = 2e-3, niter = 10, nmajor = 2)
        finally
            NFFT._use_threads[] = old
        end
    end
end

end 
