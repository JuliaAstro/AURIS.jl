"""
    dirty_beam(uv, wt; npix, cell) → Matrix{Float64}

PSF only (adjoint NFFT of the weights), normalised to 1 at the centre.  Grid it
at `npix = 2 × image npix` for use with `hogbom_clean`.
"""
function dirty_beam(uv::AbstractMatrix, wt::AbstractVector;
                    npix::Int, cell::Real)
    d = dirty_image(uv, ones(ComplexF64, length(wt)), wt; npix, cell)
    d.psf
end

"""
    hogbom_clean(dirty, psf; gain=0.1, niter=1000, threshold=0.0)
        → (model, residual, iters)

Högbom CLEAN: repeatedly locate the absolute peak of the residual, add
`gain × peak` to the model at that pixel, and subtract the centred-PSF copy.
Stops after `niter` iterations or when |peak| < `threshold`.  `psf` must be at
least `2×size(dirty) − 1` in each axis (see `dirty_beam`).
"""
function hogbom_clean(dirty::AbstractMatrix, psf::AbstractMatrix;
                      gain::Real=0.1, niter::Int=1000, threshold::Real=0.0)
    all(size(psf) .>= 2 .* size(dirty) .- 1) ||
        throw(ArgumentError("psf must be ≥ 2*size(dirty)-1 (grid it at 2×npix)"))
    res   = Matrix{Float64}(dirty)
    model = zeros(size(res))
    pc    = CartesianIndex(size(psf) .÷ 2 .+ 1)   # psf centre
    iters = 0
    pk    = _abspeak(res)
    for _ in 1:niter
        val = res[pk]
        abs(val) <= threshold && break
        iters += 1
        model[pk] += gain * val
        pk = _subtract_peak!(res, psf, pc - pk, gain * val)
    end
    model, res, iters
end

function _abspeak(res::AbstractMatrix)
    bv, bi = -1.0, CartesianIndex(1, 1)
    @inbounds for I in CartesianIndices(res)
        a = abs(res[I])
        a > bv && ((bv, bi) = (a, I))
    end
    bi
end

function _subtract_peak!(res::Matrix{Float64}, psf::AbstractMatrix,
                         off::CartesianIndex{2}, a::Float64)
    n2 = size(res, 2)
    o1, o2 = Tuple(off)
    nchunk = clamp(length(res) ÷ 1_000_000, 1, Threads.nthreads())
    nchunk == 1 && return _subtract_peak_cols!(res, psf, o1, o2, a, 1:n2)[2]
    tasks = map(Iterators.partition(1:n2, cld(n2, nchunk))) do js
        Threads.@spawn _subtract_peak_cols!(res, psf, o1, o2, a, js)
    end
    best = fetch.(tasks)
    best[argmax(first.(best))][2]
end

function _subtract_peak_cols!(res::Matrix{Float64}, psf::AbstractMatrix,
                              o1::Int, o2::Int, a::Float64, js)
    bv, bi = -1.0, CartesianIndex(1, first(js))
    @inbounds for j in js, i in 1:size(res, 1)
        r = res[i, j] - a * psf[i + o1, j + o2]
        res[i, j] = r
        ar = abs(r)
        ar > bv && ((bv, bi) = (ar, CartesianIndex(i, j)))
    end
    (bv, bi)
end

"""
    fit_beam(psf, cell) → (bmaj, bmin, bpa)

Elliptical-Gaussian fit to the PSF main lobe (least squares on log(psf) over
the pixels above 0.35 of the peak around the centre).  Returns FWHM major/minor
axes in radians and the position angle in radians.
"""
function fit_beam(psf::AbstractMatrix, cell::Real)
    c0 = size(psf) .÷ 2 .+ 1
    box = max(3, round(Int, size(psf, 1) / 16))
    xs = Float64[]; ys = Float64[]; zs = Float64[]
    for i in max(1, c0[1]-box):min(size(psf,1), c0[1]+box),
        j in max(1, c0[2]-box):min(size(psf,2), c0[2]+box)
        p = psf[i, j]
        p > 0.35 || continue
        push!(xs, (i - c0[1]) * cell)   
        push!(ys, (j - c0[2]) * cell)   
        push!(zs, -log(p))              
    end
    length(zs) >= 6 || error("Too few main-lobe pixels to fit a beam; oversample the PSF")
    A = [xs .^ 2 xs .* ys ys .^ 2]
    a, b, c = A \ zs
    tr, dt = a + c, a * c - b^2 / 4
    dt > 0 || error("Beam fit is not elliptical (unphysical PSF core)")
    λ1 = tr / 2 - sqrt(tr^2 / 4 - dt)   
    λ2 = tr / 2 + sqrt(tr^2 / 4 - dt)
    fwhm(λ) = 2 * sqrt(log(2) / λ)
    vl, vm = b / 2, λ1 - a
    (vl == 0 && vm == 0) && ((vl, vm) = (λ1 - c, b / 2))
    bpa = atan(vl, vm)
    fwhm(λ1), fwhm(λ2), bpa
end

"""
    restore(model, residual, beam, cell) → Matrix{Float64}

Convolve the CLEAN component `model` with the elliptical-Gaussian restoring
`beam` (from `fit_beam`) and add the `residual`.
"""
function restore(model::AbstractMatrix, residual::AbstractMatrix,
                 beam::Tuple{<:Real,<:Real,<:Real}, cell::Real)
    bmaj, bmin, bpa = beam
    out = Matrix{Float64}(residual)
    sx = bmaj / (2 * sqrt(2 * log(2)))   
    sy = bmin / (2 * sqrt(2 * log(2)))
    rad = ceil(Int, 4 * max(sx, sy) / cell)
    s, c = sincos(bpa)
    for I in CartesianIndices(model)
        f = model[I]
        f == 0 && continue
        for di in -rad:rad, dj in -rad:rad
            i, j = Tuple(I) .+ (di, dj)
            (1 <= i <= size(out, 1) && 1 <= j <= size(out, 2)) || continue
            l, m = di * cell, dj * cell
            xr =  l * s + m * c    
            yr =  l * c - m * s
            out[i, j] += f * exp(-0.5 * ((xr / sx)^2 + (yr / sy)^2))
        end
    end
    out
end