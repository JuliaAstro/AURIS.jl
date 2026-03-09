module BFitting

using LinearAlgebra, Statistics
using ..BBasis: reflect_pad_2d, bezier_basis, safe_median, get_tile_indices, get_sliding_window_indices, tile_window
using ..BProbability: bayesian_prob_good

export fit_tile_local_pad, fit_window_local_pad, fit_bezier_surface_bayesian_weighted, fit_bezier_surface_bayesian_weighted_window

### Tiles
function fit_tile_local_pad(tile::AbstractMatrix{T}, w_tile; pad_t::Int=6, pad_f::Int=6, deg_time::Int=3, deg_freq::Int=3) where T
    tilef = Array{Float64}(tile)
    tile_nt, tile_nf = size(tilef)
    pad_t_eff = min(pad_t, tile_nt - 1)
    pad_f_eff = min(pad_f, tile_nf - 1)

    padded = reflect_pad_2d(tilef, pad_t_eff, pad_f_eff)
    if w_tile === nothing
        w_pad = nothing
    else
        wpad_tmp = Array{Float64}(w_tile)
        w_pad = reflect_pad_2d(wpad_tmp, pad_t_eff, pad_f_eff)
        w_pad .= clamp.(w_pad, 1e-12, 1.0)
    end

    ntp, ncp = size(padded)
    u = range(0.0, 1.0, length=ncp)
    v = range(0.0, 1.0, length=ntp)
    Bu = [bezier_basis(vj, deg_freq) for vj in u]
    Bv = [bezier_basis(vi, deg_time) for vi in v]

    nparam = (deg_time + 1) * (deg_freq + 1)
    A = zeros(Float64, ntp * ncp, nparam)
    row = 1
    for i in 1:ntp, j in 1:ncp
        idx = 1
        for ii in 1:(deg_time+1), jj in 1:(deg_freq+1)
            A[row, idx] = Bv[i][ii] * Bu[j][jj]
            idx += 1
        end
        row += 1
    end

    yvec = vec(padded)
    if !all(isfinite, yvec)
        fillval = safe_median(yvec[isfinite.(yvec)])
        yvec .= ifelse.(isfinite.(yvec), yvec, fillval)
    end

    if w_pad === nothing
        λ = 1e-10 * maximum(diag(A' * A))
        M = A' * A + λ * I(size(A, 2))
        Pvec = M \ (A' * yvec)
    else
        W = Diagonal(vec(w_pad))
        λ = 1e-8 * maximum(diag(A' * W * A))
        M = A' * W * A + λ * I(size(A, 2))
        RHS = A' * W * yvec
        if any(.!isfinite.(M)) || any(.!isfinite.(RHS))
            fillval = safe_median(yvec)
            yvec .= ifelse.(isfinite.(yvec), yvec, fillval)
            M = A' * W * A + λ * I(size(A, 2))
            RHS = A' * W * yvec
        end
        Pvec = M \ RHS
    end

    padded_fit = reshape(A * Pvec, ntp, ncp)
    crop = padded_fit[pad_t_eff+1:pad_t_eff+tile_nt, pad_f_eff+1:pad_f_eff+tile_nf]
    return crop
end

function fit_bezier_surface_bayesian_weighted(y::AbstractMatrix{T};
    deg_time::Int=3,
    deg_freq::Int=3,
    tiles_time::Int=6,
    tiles_freq::Int=8,
    overlap_frac::Float64=0.5,
    maxiter::Int=8,
    tol::Float64=1e-5,
    β_floor::Float64=1e-4,
    global_pad_t::Int=6,
    global_pad_f::Int=6,
    local_pad_t::Int=2,
    local_pad_f::Int=2
) where T

    ymat = Array{Float64}(copy(y))
    nt, nchan = size(ymat)

    if nt < 3 || nchan < 3
        return copy(ymat), ones(size(ymat))
    end

    if !all(isfinite, ymat)
        fillval = safe_median(vec(ymat)[isfinite.(vec(ymat))])
        ymat .= ifelse.(isfinite.(ymat), ymat, fillval)
    end

    y_padded = reflect_pad_2d(ymat, global_pad_t, global_pad_f)
    nt_pad, nchan_pad = size(y_padded)

    yfit_surface = zeros(Float64, nt_pad, nchan_pad)
    weights_accum = zeros(Float64, nt_pad, nchan_pad)

    t_tiles = get_tile_indices(nt_pad, tiles_time, overlap_frac)
    f_tiles = get_tile_indices(nchan_pad, tiles_freq, overlap_frac)

    for (ti, ti1) in t_tiles, (fj, fj1) in f_tiles
        tile = y_padded[ti:ti1, fj:fj1]
        yfit_tile = fit_tile_local_pad(tile, nothing; pad_t=local_pad_t, pad_f=local_pad_f, deg_time=deg_time, deg_freq=deg_freq)

        win = tile_window(size(yfit_tile, 1), size(yfit_tile, 2))
        yfit_surface[ti:ti1, fj:fj1] .+= yfit_tile .* win
        weights_accum[ti:ti1, fj:fj1] .+= win
    end

    weights_accum .= ifelse.(weights_accum .== 0.0, 1.0, weights_accum)
    yfit_surface ./= weights_accum
    prev_surface = copy(yfit_surface)

    for iter in 1:maxiter
        resid = y_padded .- yfit_surface
        rvec = vec(resid)

        rvec_finite = rvec[isfinite.(rvec)]
        if isempty(rvec_finite)
            σ_global = 1.0
        else
            σ_global = 1.4826 * median(abs.(rvec_finite .- median(rvec_finite))) + eps()
        end

        p_good = zeros(Float64, nt_pad, nchan_pad)
        wtemp = zeros(Float64, nt_pad, nchan_pad)

        for (ti, ti1) in t_tiles, (fj, fj1) in f_tiles
            resid_tile = resid[ti:ti1, fj:fj1]
            rvec_tile = vec(resid_tile)

            if !all(isfinite, rvec_tile)
                fillval = safe_median(rvec_tile[isfinite.(rvec_tile)])
                rvec_tile .= ifelse.(isfinite.(rvec_tile), rvec_tile, fillval)
            end

            med_r_tile = median(rvec_tile)
            absr = abs.(rvec_tile .- med_r_tile)
            q80 = quantile(absr, 0.8)
            tr = absr[absr.<q80]
            tr = isempty(tr) ? absr : tr
            σ_local = 1.4826 * median(tr) + eps()
            σ_local = isfinite(σ_local) && σ_local > 0 ? σ_local : σ_global

            ac = 0.0
            if length(rvec_tile) > 2
                x = rvec_tile .- mean(rvec_tile)
                denom = sum(x .^ 2)
                if denom > 0
                    ac = sum(x[1:end-1] .* x[2:end]) / denom
                end
                ac = clamp(ac, -0.99, 0.99)
            end
            σ_eff = σ_local * sqrt(1 + max(ac, 0.0))

            tail_frac = mean(absr .> 3.0 * σ_eff)
            β_tile = clamp(tail_frac, 0.001, 0.2)

            p_tile = bayesian_prob_good(rvec_tile .- med_r_tile, σ_eff; β=β_tile, γ=5.0)
            p_tile = reshape(p_tile, size(resid_tile))

            p_good[ti:ti1, fj:fj1] .+= p_tile
            wtemp[ti:ti1, fj:fj1] .+= 1.0
        end

        p_good ./= ifelse.(wtemp .== 0.0, 1.0, wtemp)
        p_good .= clamp.(p_good .^ 2.5, 1e-8, 1.0)

        denom = max(σ_global^2, 1e-12)
        wmat = clamp.(p_good ./ denom, 1e-8, 1.0)

        yfit_surface .= 0.0
        weights_accum .= 0.0

        for (ti, ti1) in t_tiles, (fj, fj1) in f_tiles
            tile = y_padded[ti:ti1, fj:fj1]
            w_tile = wmat[ti:ti1, fj:fj1]
            yfit_tile = fit_tile_local_pad(tile, w_tile; pad_t=local_pad_t, pad_f=local_pad_f, deg_time=deg_time, deg_freq=deg_freq)

            win = tile_window(size(yfit_tile, 1), size(yfit_tile, 2))
            yfit_surface[ti:ti1, fj:fj1] .+= yfit_tile .* win
            weights_accum[ti:ti1, fj:fj1] .+= win
        end

        weights_accum .= ifelse.(weights_accum .== 0.0, 1.0, weights_accum)
        yfit_surface ./= weights_accum

        relchg = norm(yfit_surface - prev_surface) / (norm(prev_surface) + 1e-12)
        prev_surface .= yfit_surface
        if relchg < tol
            break
        end
    end

    resid = y_padded .- yfit_surface
    rvec = vec(resid)
    rvec_finite = rvec[isfinite.(rvec)]
    if isempty(rvec_finite)
        med_r = 0.0
        σ_global = 1.0
    else
        med_r = median(rvec_finite)
        σ_global = 1.4826 * median(abs.(rvec_finite .- med_r)) + eps()
    end

    tail_frac = mean(abs.(rvec_finite .- med_r) .> 3.0 * σ_global)
    β = clamp(tail_frac, β_floor, 0.2)
    γ = clamp(quantile(abs.(rvec_finite .- med_r), 0.9) / σ_global, 3.0, 50.0)

    pvals = bayesian_prob_good(rvec .- med_r, σ_global; β=β, γ=γ)
    p_good = reshape(clamp.(pvals, 1e-12, 1.0), nt_pad, nchan_pad)
    p_good[resid.≤0] .= 1.0

    yfit_surface = yfit_surface[global_pad_t+1:end-global_pad_t, global_pad_f+1:end-global_pad_f]
    p_good = p_good[global_pad_t+1:end-global_pad_t, global_pad_f+1:end-global_pad_f]

    return yfit_surface, p_good
end

### Windows
function fit_window_local_pad(tile::AbstractMatrix{T}, w_tile;
    deg_time::Int=3, deg_freq::Int=3
) where T

    y = Array{Float64}(tile)
    nt, nf = size(y)

    u = range(0.0, 1.0, length=nf)
    v = range(0.0, 1.0, length=nt)
    Bu = [bezier_basis(uj, deg_freq) for uj in u]
    Bv = [bezier_basis(vi, deg_time) for vi in v]

    nparam = (deg_time + 1) * (deg_freq + 1)
    A = zeros(Float64, nt * nf, nparam)

    row = 1
    for i in 1:nt, j in 1:nf
        idx = 1
        for ii in 1:(deg_time+1), jj in 1:(deg_freq+1)
            A[row, idx] = Bv[i][ii] * Bu[j][jj]
            idx += 1
        end
        row += 1
    end

    yvec = vec(y)
    if !all(isfinite, yvec)
        fillval = safe_median(yvec[isfinite.(yvec)])
        yvec .= ifelse.(isfinite.(yvec), yvec, fillval)
    end

    if w_tile === nothing
        λ = 1e-10 * maximum(diag(A' * A))
        M = A' * A + λ * I(size(A, 2))
        P = M \ (A' * yvec)
    else
        W = Diagonal(vec(w_tile))
        λ = 1e-8 * maximum(diag(A' * W * A))
        M = A' * W * A + λ * I(size(A, 2))
        RHS = A' * W * yvec
        P = M \ RHS
    end

    reshape(A * P, nt, nf)
end

function fit_bezier_surface_bayesian_weighted_window(
    y::AbstractMatrix{T};
    deg_time::Int=3,
    deg_freq::Int=3,
    window_time::Int=64,
    window_freq::Int=32,
    maxiter::Int=8,
    tol::Float64=1e-5,
    β_floor::Float64=1e-4
) where T

    ymat = Array{Float64}(y)
    nt, nf = size(ymat)
    fillval = safe_median(vec(ymat)[isfinite.(vec(ymat))])
    ymat .= ifelse.(isfinite.(ymat), ymat, fillval)

    t_windows = get_sliding_window_indices(nt, window_time)
    f_windows = get_sliding_window_indices(nf, window_freq)

    yfit = zeros(nt, nf)
    wacc = zeros(nt, nf)

    for (ti, ti1) in t_windows, (fj, fj1) in f_windows
        tile = ymat[ti:ti1, fj:fj1]
        f = fit_window_local_pad(tile, nothing; deg_time=deg_time, deg_freq=deg_freq)
        yfit[ti:ti1, fj:fj1] .+= f
        wacc[ti:ti1, fj:fj1] .+= 1
    end
    yfit ./= max.(wacc, 1)
    prev = copy(yfit)

    for iter in 1:maxiter
        resid = ymat .- yfit
        rvec = vec(resid)
        med_r = median(rvec)
        σg = 1.4826 * median(abs.(rvec .- med_r)) + eps()

        p_good = zeros(nt, nf)
        wcnt = zeros(nt, nf)

        for (ti, ti1) in t_windows, (fj, fj1) in f_windows
            rt = resid[ti:ti1, fj:fj1]
            rv = vec(rt)
            medt = median(rv)
            σloc = 1.4826 * median(abs.(rv .- medt)) + eps()

            tail_frac = mean(abs.(rv .- medt) .> 3σloc)
            β = clamp(tail_frac, 0.001, 0.2)

            pt = bayesian_prob_good(rv .- medt, σloc; β=β, γ=5.0)
            p_good[ti:ti1, fj:fj1] .+= reshape(pt, size(rt))
            wcnt[ti:ti1, fj:fj1] .+= 1
        end

        p_good ./= max.(wcnt, 1)
        p_good .= clamp.(p_good .^ 1.0, 1e-8, 1.0)

        wmat = clamp.(p_good ./ (σg^2), 1e-8, 1.0)

        yfit .= 0
        wacc .= 0

        for (ti, ti1) in t_windows, (fj, fj1) in f_windows
            tile = ymat[ti:ti1, fj:fj1]
            wt = wmat[ti:ti1, fj:fj1]
            f = fit_window_local_pad(tile, wt; deg_time=deg_time, deg_freq=deg_freq)
            yfit[ti:ti1, fj:fj1] .+= f
            wacc[ti:ti1, fj:fj1] .+= 1
        end

        yfit ./= max.(wacc, 1)
        rel = norm(yfit - prev) / (norm(prev) + 1e-12)
        prev .= yfit
        if rel < tol
            break
        end
    end

    resid = ymat .- yfit
    rvec = vec(resid)
    med_r = median(rvec)
    σ = 1.4826 * median(abs.(rvec .- med_r)) + eps()

    tail_frac = mean(abs.(rvec .- med_r) .> 3σ)
    β = clamp(tail_frac, β_floor, 0.2)
    γ = clamp(quantile(abs.(rvec .- med_r), 0.9) / σ, 3, 50)

    p = bayesian_prob_good(rvec .- med_r, σ; β=β, γ=γ)
    p = reshape(p, nt, nf)
    p[resid.≤0] .= 1.0

    return yfit, p
end

end