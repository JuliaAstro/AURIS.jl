function n_minus_1(npix::Int, cell::Real)
    c0 = npix ÷ 2 + 1
    [sqrt(max(0.0, 1 - ((il - c0) * cell)^2 - ((im - c0) * cell)^2)) - 1
     for il in 1:npix, im in 1:npix]
end

function nw_auto(w::AbstractVector, npix::Int, cell::Real; tol::Real=0.05)
    isempty(w) && return 1
    wspan = maximum(w) - minimum(w)
    r2 = 2 * (npix ÷ 2 * cell)^2                 
    n1max = 1 - sqrt(max(0.0, 1 - r2))
    clamp(ceil(Int, π * wspan * n1max / tol), 1, 256)
end

function w_bins(w::AbstractVector, nw::Int)
    wmin, wmax = extrema(w)
    Δ = (wmax - wmin) / nw
    idx = Δ > 0 ? [clamp(floor(Int, (x - wmin) / Δ) + 1, 1, nw) for x in w] :
                  ones(Int, length(w))
    centres = [wmin + (b - 0.5) * Δ for b in 1:nw]
    idx, centres
end

"""
    dirty_image_w(uvw, vis, wt; npix=256, cell=nyquist_cell(uvw)/3, nw=0)
        → (image, psf, cell, sumwt, nw)

Naturally-weighted Stokes-I dirty image and PSF with the w-term, by
w-stacking over `nw` planes.
"""
function dirty_image_w(uvw::AbstractMatrix, vis::AbstractVector,
                       wt::AbstractVector;
                       npix::Int=256, cell::Real=nyquist_cell(uvw[1:2, :]) / 3,
                       nw::Int=0)
    size(uvw, 1) == 3 || throw(DimensionMismatch("uvw must be 3×M"))
    iseven(npix) || throw(ArgumentError("npix must be even"))
    x = Matrix{Float64}(uvw[1:2, :] .* cell)
    mx = maximum(abs, x)
    mx < 0.5 || throw(ArgumentError(
        "uv samples exceed the grid Nyquist limit (max node $mx ≥ 0.5); " *
        "use cell ≤ $(cell * 0.499 / mx)"))
    nw == 0 && (nw = nw_auto(vec(uvw[3, :]), npix, cell))

    xh   = hcat(x, -x)
    wh   = vcat(Vector{Float64}(uvw[3, :]), -Vector{Float64}(uvw[3, :]))
    valh = vcat(wt .* vis, wt .* conj.(vis))
    wth  = vcat(Vector{ComplexF64}(complex.(wt)), Vector{ComplexF64}(complex.(wt)))

    idx, centres = w_bins(wh, nw)
    n1 = n_minus_1(npix, cell)
    acc  = zeros(ComplexF64, npix, npix)
    pacc = zeros(ComplexF64, npix, npix)
    for b in 1:nw
        sel = findall(==(b), idx)
        isempty(sel) && continue
        p = plan_nfft(xh[:, sel], (npix, npix))
        ph = cispi.(2 .* centres[b] .* n1)
        acc  .+= (adjoint(p) * valh[sel]) .* ph
        pacc .+= (adjoint(p) * wth[sel]) .* ph
    end
    sumwt = 2 * sum(wt)
    (image = real.(acc) ./ sumwt, psf = real.(pacc) ./ sumwt,
     cell = Float64(cell), sumwt = sumwt, nw = nw)
end

"""
    dirty_image_wdft(uvw, vis, wt; npix, cell) → (image, psf, cell, sumwt)

Direct-summation reference for `dirty_image_w`, the exact w-term transform
with no plane quantisation.
"""
function dirty_image_wdft(uvw::AbstractMatrix, vis::AbstractVector,
                          wt::AbstractVector; npix::Int=64, cell::Real)
    img = zeros(npix, npix)
    psf = zeros(npix, npix)
    c0  = npix ÷ 2 + 1
    for il in 1:npix, im in 1:npix
        l = (il - c0) * cell
        m = (im - c0) * cell
        n1 = sqrt(max(0.0, 1 - l^2 - m^2)) - 1
        acc  = 0.0
        pacc = 0.0
        for j in eachindex(vis)
            ph = cispi(2 * (uvw[1, j] * l + uvw[2, j] * m + uvw[3, j] * n1))
            acc  += wt[j] * real(vis[j] * ph)
            pacc += wt[j] * real(ph)
        end
        img[il, im] = acc
        psf[il, im] = pacc
    end
    sumwt = sum(wt)
    (image=img ./ sumwt, psf=psf ./ sumwt, cell=Float64(cell), sumwt=sumwt)
end

"""
    predict_vis_w(model, uvw, cell; nw=0) → Vector{ComplexF64}

Model visibilities of a CLEAN component image at the `uvw` points
(wavelengths), including the w-term.
"""
function predict_vis_w(model::AbstractMatrix{<:Real}, uvw::AbstractMatrix,
                       cell::Real; nw::Int=0)
    size(uvw, 1) == 3 || throw(DimensionMismatch("uvw must be 3×M"))
    npix = size(model, 1)
    w = Vector{Float64}(uvw[3, :])
    nw == 0 && (nw = nw_auto(w, npix, cell))
    x = Matrix{Float64}(uvw[1:2, :] .* cell)

    idx, centres = w_bins(w, nw)
    n1 = n_minus_1(npix, cell)
    out = Vector{ComplexF64}(undef, size(uvw, 2))
    for b in 1:nw
        sel = findall(==(b), idx)
        isempty(sel) && continue
        p = plan_nfft(x[:, sel], size(model))
        out[sel] = p * (Matrix{ComplexF64}(model) .* cispi.(-2 .* centres[b] .* n1))
    end
    out
end

"""
    cs_clean_w(uvw, vis, wt; npix, cell, gain=0.1, niter=2000, threshold=0.0,
               nmajor=10, major_frac=0.15, nw=0)
        → (model, residual, restored, beam, iters, nw)

Cotton–Schwab CLEAN with the w-term: as `cs_clean`, but gridding with
`dirty_image_w` and degridding with `predict_vis_w` over `nw` w-planes.
"""
function cs_clean_w(uvw::AbstractMatrix, vis::AbstractVector, wt::AbstractVector;
                    npix::Int, cell::Real, gain::Real=0.1, niter::Int=2000,
                    threshold::Real=0.0, nmajor::Int=10, major_frac::Real=0.15,
                    nw::Int=0)
    nw == 0 && (nw = nw_auto(vec(uvw[3, :]), npix, cell))
    psf2  = dirty_image_w(uvw, ones(ComplexF64, length(wt)), wt;
                          npix = 2 * npix, cell, nw).psf
    d0    = dirty_image_w(uvw, vis, wt; npix, cell, nw)
    beam  = fit_beam(d0.psf, cell)

    model = zeros(npix, npix)
    res   = d0.image
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
        resid_vis = Vector{ComplexF64}(vis) .- predict_vis_w(model, uvw, cell; nw)
        res = dirty_image_w(uvw, resid_vis, wt; npix, cell, nw).image
        total >= niter && break
    end
    restored = restore(model, res, beam, cell)
    (model = model, residual = res, restored = restored, beam = beam,
     iters = total, nw = nw)
end
