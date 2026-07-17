const C_LIGHT = 299_792_458.0   # m/s

"""
    uvw_samples(v::VisibilityDataset; time_stride=1, chan_stride=1)
        → (uvw, vis, wt)

Flatten a `VisibilityDataset` into Stokes-I continuum samples for imaging:
`uvw` (3×M, wavelengths — each channel scales the metre baseline by ν/c),
`vis` (ComplexF64, I = (RR+LL)/2) and `wt` (Float64).  A sample is included
only if both parallel-hand products are unflagged.  `time_stride`/`chan_stride`
subsample for quick looks.
"""
function uvw_samples(v::VisibilityDataset; time_stride::Int=1, chan_stride::Int=1)
    npol = Data.n_pol(v)
    rr, ll = 1, npol            
    nt, nbl, nspw, nchan = Data.n_time(v), Data.n_baseline(v),
                           Data.n_spw(v), Data.n_chan(v)

    us  = Float64[];  vs = Float64[];  ws = Float64[]
    vis = ComplexF64[]; wt = Float64[]
    sizehint = length(1:time_stride:nt) * nbl * nspw * length(1:chan_stride:nchan)
    foreach(a -> Base.sizehint!(a, sizehint), (us, vs, ws, vis, wt))

    for t in 1:time_stride:nt, bl in 1:nbl, s in 1:nspw
        u_m, v_m, w_m = v.uvw[1, bl, t], v.uvw[2, bl, t], v.uvw[3, bl, t]
        w_s = Float64(v.weights[s, bl, t])
        for c in 1:chan_stride:nchan
            (v.flags[rr, c, s, bl, t] || v.flags[ll, c, s, bl, t]) && continue
            z = 0.5 * (ComplexF64(v.vis[rr, c, s, bl, t]) +
                       ComplexF64(v.vis[ll, c, s, bl, t]))
            isfinite(real(z)) && isfinite(imag(z)) || continue   
            scale = v.freqs[c, s] / C_LIGHT
            push!(us, u_m * scale)
            push!(vs, v_m * scale)
            push!(ws, w_m * scale)
            push!(vis, z)
            push!(wt, 2 * w_s)   
        end
    end
    permutedims([us vs ws]), vis, wt
end

function uv_samples(v::VisibilityDataset; kw...)
    uvw, vis, wt = uvw_samples(v; kw...)
    uvw[1:2, :], vis, wt
end

nyquist_cell(uv::AbstractMatrix) = 1 / (2 * maximum(abs, uv))

"""
    dirty_image(uv, vis, wt; npix=256, cell=nyquist_cell(uv)/3)
        → (image, psf, cell, sumwt)

Naturally-weighted dirty image and PSF on an npix×npix grid of `cell` radians
via adjoint NFFT.  The PSF is normalised to 1 at the phase centre. The image
carries the same normalisation (Σw), so a point source at the phase centre
appears with its correlated flux.
"""
function dirty_image(uv::AbstractMatrix, vis::AbstractVector, wt::AbstractVector;
                     npix::Int=256, cell::Real=nyquist_cell(uv) / 3)
    iseven(npix) || throw(ArgumentError("npix must be even"))
    x = Matrix{Float64}(uv .* cell)             
    mx = maximum(abs, x)
    mx < 0.5 || throw(ArgumentError(
        "uv samples exceed the grid Nyquist limit (max node $mx ≥ 0.5); " *
        "use cell ≤ $(cell * 0.499 / mx)"))

    p     = plan_nfft(x, (npix, npix))
    sumwt = sum(wt)
    img   = real.(adjoint(p) * Vector{ComplexF64}(wt .* vis)) ./ sumwt
    psf   = real.(adjoint(p) * Vector{ComplexF64}(complex.(wt))) ./ sumwt
    (image=img, psf=psf, cell=Float64(cell), sumwt=sumwt)
end

function dirty_image_dft(uv::AbstractMatrix, vis::AbstractVector, wt::AbstractVector;
                         npix::Int=64, cell::Real=nyquist_cell(uv) / 3)
    img = zeros(npix, npix)
    psf = zeros(npix, npix)
    c0  = npix ÷ 2 + 1
    for il in 1:npix, im in 1:npix
        l = (il - c0) * cell
        m = (im - c0) * cell
        acc  = 0.0 + 0.0im
        pacc = 0.0 + 0.0im
        for j in eachindex(vis)
            ph = cispi(2 * (uv[1, j] * l + uv[2, j] * m))
            acc  += wt[j] * vis[j] * ph
            pacc += wt[j] * ph
        end
        img[il, im] = real(acc)
        psf[il, im] = real(pacc)
    end
    sumwt = sum(wt)
    (image=img ./ sumwt, psf=psf ./ sumwt, cell=Float64(cell), sumwt=sumwt)
end
