module BAPI

using ..BFitting: fit_bezier_surface_bayesian_weighted, fit_bezier_surface_bayesian_weighted_window

export flag_tile, flag_window

function flag_tile(vis;
    surf_iter::Int=8,
    deg_time::Int=3,
    deg_freq::Int=3,
    tiles_time::Int=6,
    tiles_freq::Int=8,
    overlap_frac::Float64=0.5,
    p_thresh::Float64=0.5
)

    if ndims(vis) == 2
        vis3 = reshape(vis, 1, size(vis, 1), size(vis, 2))
    elseif ndims(vis) == 3
        vis3 = vis
    else
        throw(ArgumentError("vis must be a 2D (ntxnchan) or 3D (nblxntxnchan) array"))
    end

    nbl, nt, nchan = size(vis3)
    mask = falses(nbl, nt, nchan)
    p_good_all = ones(Float64, nbl, nt, nchan)
    yfit_surfaces = zeros(Float64, nbl, nt, nchan)

    for b in 1:nbl
        y = abs.(vis3[b, :, :])

        yfit, p_good = fit_bezier_surface_bayesian_weighted(
            y;
            deg_time=deg_time,
            deg_freq=deg_freq,
            tiles_time=tiles_time,
            tiles_freq=tiles_freq,
            overlap_frac=overlap_frac,
            maxiter=surf_iter
        )

        yfit_surfaces[b, :, :] .= yfit
        p_good_all[b, :, :] .= p_good
        mask[b, :, :] .= p_good .< p_thresh
    end

    return mask, p_good_all, yfit_surfaces
end

function flag_window(vis;
    surf_iter::Int=8,
    deg_time::Int=3,
    deg_freq::Int=3,
    window_time::Int=64,
    window_freq::Int=32,
    p_thresh::Float64=0.1
)

    vis3 = ndims(vis) == 2 ? reshape(vis, 1, size(vis, 1), size(vis, 2)) : vis
    nbl, nt, nf = size(vis3)

    mask = falses(nbl, nt, nf)
    pgood_all = ones(Float64, nbl, nt, nf)
    yfits = zeros(Float64, nbl, nt, nf)

    for b in 1:nbl
        y = abs.(vis3[b, :, :])
        yfit, pg = fit_bezier_surface_bayesian_weighted_window(y;
            deg_time=deg_time, deg_freq=deg_freq,
            window_time=window_time, window_freq=window_freq,
            maxiter=surf_iter)

        mask[b, :, :] .= pg .< p_thresh
        pgood_all[b, :, :] .= pg
        yfits[b, :, :] .= yfit
    end

    return mask, pgood_all, yfits
end

end