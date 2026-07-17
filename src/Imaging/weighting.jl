"""
    briggs_weights(uv, wt; npix, cell, robust=0.0) → Vector{Float64}

Briggs (1995) robust reweighting of the imaging weights `wt` for samples at
`uv` (2×M, wavelengths), using the weight density on the same npix×npix grid
of `cell` radians that the image will be made on.  `robust` ∈ [−2, 2]:
−2 ≈ uniform, +2 ≈ natural.  Returns the new weight vector (the input is not
modified); pass it to `dirty_image`/`cs_clean` in place of `wt`.
"""
function briggs_weights(uv::AbstractMatrix, wt::AbstractVector;
                        npix::Int=256, cell::Real=nyquist_cell(uv) / 3,
                        robust::Real=0.0)
    size(uv, 1) == 2 && size(uv, 2) == length(wt) ||
        throw(DimensionMismatch("uv must be 2×M with M = length(wt)"))
    iseven(npix) || throw(ArgumentError("npix must be even"))
    mx = maximum(abs, uv) * cell
    mx < 0.5 || throw(ArgumentError("uv·cell must lie in [-0.5, 0.5); " *
        "use cell ≤ $(cell * 0.499 / mx)"))

    dens = zeros(npix, npix)
    cellindex(x) = clamp(floor(Int, (x + 0.5) * npix) + 1, 1, npix)
    M = length(wt)
    iu = Vector{Int}(undef, M)
    iv = Vector{Int}(undef, M)
    for i in 1:M
        iu[i] = cellindex(uv[1, i] * cell)
        iv[i] = cellindex(uv[2, i] * cell)
        dens[iu[i], iv[i]] += wt[i]
        dens[cellindex(-uv[1, i] * cell), cellindex(-uv[2, i] * cell)] += wt[i]
    end

    sumw  = 2 * sum(wt)
    sumd2 = sum(abs2, dens)
    f2 = (5.0 * 10.0^(-robust))^2 / (sumd2 / sumw)
    [wt[i] / (1 + dens[iu[i], iv[i]] * f2) for i in 1:M]
end
