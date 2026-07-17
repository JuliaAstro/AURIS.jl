"""
    predict_vis(model, uv, cell) → Vector{ComplexF64}

Model visibilities of a CLEAN component image `model` (npix×npix, same grid
conventions as `dirty_image`) at the uv points (wavelengths).
"""
function predict_vis(model::AbstractMatrix{<:Real}, uv::AbstractMatrix, cell::Real)
    x = Matrix{Float64}(uv .* cell)
    p = plan_nfft(x, size(model))
    p * Matrix{ComplexF64}(model)
end

"""
    predict_dataset_vis(v::VisibilityDataset, model, cell)
        → Array{ComplexF32,4}

Model visibilities of the component image `model` (npix×npix, grid conventions
of `dirty_image`) for every sample of `v`, shaped `(chan, spw, baseline,
time)` — the input `solve_selfcal` needs.  Each channel scales the metre
baselines by ν/c, exactly as `uv_samples` does.  Prediction runs in time
chunks to bound NFFT memory.
"""
function predict_dataset_vis(v::VisibilityDataset, model::AbstractMatrix{<:Real},
                             cell::Real; chunk_samples::Int=4_000_000)
    nchan, nspw = Data.n_chan(v), Data.n_spw(v)
    nbl, nt = Data.n_baseline(v), Data.n_time(v)
    per_t = nchan * nspw * nbl
    tchunk = max(1, chunk_samples ÷ per_t)

    out = Array{ComplexF32,4}(undef, nchan, nspw, nbl, nt)
    x = Matrix{Float64}(undef, 2, per_t * tchunk)
    mimg = Matrix{ComplexF64}(model)
    for t0 in 1:tchunk:nt
        ts = t0:min(t0 + tchunk - 1, nt)
        k = 0
        for t in ts, bl in 1:nbl, s in 1:nspw, c in 1:nchan
            k += 1
            scale = v.freqs[c, s] / C_LIGHT * cell
            x[1, k] = v.uvw[1, bl, t] * scale
            x[2, k] = v.uvw[2, bl, t] * scale
        end
        xk = view(x, :, 1:k)
        mx = maximum(abs, xk)
        mx < 0.5 || throw(ArgumentError(
            "uv samples exceed the model grid Nyquist limit; use cell ≤ $(cell * 0.499 / mx)"))
        p = plan_nfft(Matrix(xk), size(model))
        mv = p * mimg
        copyto!(view(out, :, :, :, ts), reshape(mv, nchan, nspw, nbl, length(ts)))
    end
    out
end

"""
    cs_clean(uv, vis, wt; npix, cell, gain=0.1, niter=2000, threshold=0.0,
             nmajor=10, major_frac=0.15)
        → (model, residual, restored, beam, iters)

Cotton–Schwab CLEAN: Högbom minor cycles run on the residual image down to
`major_frac × current peak` (or `threshold`), then the accumulated component
model is degridded with `predict_vis` and subtracted from the visibilities
before the next major cycle regrids the true residual.  Stops after `nmajor`
major cycles, `niter` total components, or when the residual peak reaches
`threshold`.
"""
function cs_clean(uv::AbstractMatrix, vis::AbstractVector, wt::AbstractVector;
                  npix::Int, cell::Real, gain::Real=0.1, niter::Int=2000,
                  threshold::Real=0.0, nmajor::Int=10, major_frac::Real=0.15)
    iseven(npix) || throw(ArgumentError("npix must be even"))
    x = Matrix{Float64}(uv .* cell)
    mx = maximum(abs, x)
    mx < 0.5 || throw(ArgumentError(
        "uv samples exceed the grid Nyquist limit (max node $mx ≥ 0.5); " *
        "use cell ≤ $(cell * 0.499 / mx)"))

    p     = plan_nfft(x, (npix, npix))
    p2    = plan_nfft(x, (2 * npix, 2 * npix))
    sumwt = sum(wt)
    cwt   = Vector{ComplexF64}(complex.(wt))
    psf2  = real.(adjoint(p2) * cwt) ./ sumwt
    psf   = real.(adjoint(p) * cwt) ./ sumwt
    beam  = fit_beam(psf, cell)

    model = zeros(npix, npix)
    res   = real.(adjoint(p) * Vector{ComplexF64}(wt .* vis)) ./ sumwt
    total = 0
    for _ in 1:nmajor
        peak = maximum(abs, res)
        peak <= threshold && break
        mthr = max(threshold, major_frac * peak)
        dm, res, it = hogbom_clean(res, psf2; gain, niter = niter - total,
                                   threshold = mthr)
        it == 0 && break
        total += it
        model .+= dm
        resid_vis = Vector{ComplexF64}(vis) .- p * Matrix{ComplexF64}(model)
        res = real.(adjoint(p) * (wt .* resid_vis)) ./ sumwt
        total >= niter && break
    end
    restored = restore(model, res, beam, cell)
    (model = model, residual = res, restored = restored, beam = beam,
     iters = total)
end