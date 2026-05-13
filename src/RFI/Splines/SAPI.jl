module SAPI

using LinearAlgebra, Statistics

using ..SBasis: bases, knot_vector
using ..SFitting: fit_ls, fit_em, fit_masked, evaluate
using ..SOptimisation: select_knots
using ..SProbability: prob_good, fit_beta_gamma
using ..SProjection: flag_broadband_times, flag_narrowband_freqs
using ..SFourier: flag_broadband_times_fourier

export fit_surface, apply_flags, SplineFlagBasis,
       flag_broadband_times, flag_narrowband_freqs,
       flag_broadband_times_fourier, flag_broadband_times_surface

struct SplineFlagBasis
    p             :: Int
    kx            :: Vector{Float64}
    ky            :: Vector{Float64}
    C             :: Matrix{Float64}
    σ             :: Float64
    β             :: Float64
    γ             :: Float64
    prob_threshold:: Float64
    nx            :: Int
    ny            :: Int
end

function apply_flags(basis::SplineFlagBasis, Z::AbstractMatrix)
    @assert size(Z) == (basis.nx, basis.ny) "Z size $(size(Z)) does not match basis ($(basis.nx) × $(basis.ny))"
    ε  = eps(Float64)
    ξ  = [clamp(xi, ε, 1 - ε) for xi in range(0.0, 1.0, length=basis.nx)]
    η  = [clamp(yj, ε, 1 - ε) for yj in range(0.0, 1.0, length=basis.ny)]
    Bx, By  = bases(basis.kx, basis.ky, basis.p, ξ, η)
    Z_recon = evaluate(Bx, By, basis.C)
    resid   = Z .- Z_recon
    col_med = median(resid, dims=2)
    Z_fit_adj = max.(Z_recon .+ col_med, 1e-10)
    p = prob_good(Z, Z_fit_adj, basis.σ; β=basis.β, γ=basis.γ)
    p[resid .≤ 0] .= 1.0
    return p .> basis.prob_threshold
end

function fit_surface(
    Z               :: AbstractMatrix;
    p               :: Int     = 3,
    nk_time         :: Int     = 0,
    nk_freq         :: Int     = 0,
    bic_range               = nothing,
    bic_range_freq          = nothing,
    λ                        = 1e-6,
    em_maxiters     :: Int     = 30,
    em_tol                   = 1e-6,
    β               :: Float64 = 0.01,
    γ               :: Float64 = 10.0,
    prob_threshold  :: Float64 = 0.35,
    β_range                    = exp.(range(log(0.001), log(0.5),  length=40)),
    γ_range                    = exp.(range(log(2.0),   log(5000.0), length=40)),
    bad_times                  = nothing,
    bad_freqs                  = nothing
)
    nx, ny = size(Z)
    ε      = eps(Float64)
    ξ = [clamp(xi, ε, 1 - ε) for xi in range(0.0, 1.0, length=nx)]
    η = [clamp(yj, ε, 1 - ε) for yj in range(0.0, 1.0, length=ny)]

    proj_mask = trues(nx, ny)
    if !isnothing(bad_times)
        @assert length(bad_times) == nx "bad_times length $(length(bad_times)) ≠ nx=$nx"
        proj_mask[bad_times, :] .= false
    end
    if !isnothing(bad_freqs)
        @assert length(bad_freqs) == ny "bad_freqs length $(length(bad_freqs)) ≠ ny=$ny"
        proj_mask[:, bad_freqs] .= false
    end
    has_proj = any(.!proj_mask)

    if nk_time > 0 && nk_freq > 0
        nkx, nky = nk_time, nk_freq
    else
        kx0  = knot_vector(0.0, 1.0, p+1, p)
        ky0  = knot_vector(0.0, 1.0, p+1, p)
        Bx0, By0 = bases(kx0, ky0, p, ξ, η)
        C0   = fit_ls(Bx0, By0, Z; λ=λ)
        Z0   = evaluate(Bx0, By0, C0)
        r0   = Z .- Z0
        σ0   = max(1.4826 * median(abs.(r0 .- median(r0))), 1e-10)
        pg0  = prob_good(Z, max.(Z0 .+ median(r0), 1e-10), σ0; β=β, γ=γ)
        pg0[r0 .≤ 0] .= 1.0
        rough_mask = pg0 .> prob_threshold
        has_proj && (rough_mask .&= proj_mask)

        bic_nkx, bic_nky = select_knots(Z, ξ, η, p;
            mask=rough_mask, k_range=bic_range, k_range_freq=bic_range_freq, λ=λ)
        nkx = nk_time > 0 ? nk_time : bic_nkx
        nky = nk_freq > 0 ? nk_freq : bic_nky
    end

    kx = knot_vector(0.0, 1.0, nkx, p)
    ky = knot_vector(0.0, 1.0, nky, p)
    Bx, By = bases(kx, ky, p, ξ, η)

    C     = fit_em(Bx, By, Z; λ=λ, β=β, γ=γ, maxiters=em_maxiters, tol=em_tol,
                   forced_zero = has_proj ? .!proj_mask : nothing)
    Z_fit = evaluate(Bx, By, C)

    resid = Z .- Z_fit

    col_med = zeros(nx, 1)
    for i in 1:nx
        (has_proj && !isnothing(bad_times) && bad_times[i]) && continue
        clean_cols = (!isnothing(bad_freqs) && any(bad_freqs)) ? findall(.!bad_freqs) : 1:ny
        col_med[i] = isempty(clean_cols) ? 0.0 : median(resid[i, clean_cols])
    end

    resid_c   = resid .- col_med
    Z_fit_adj = max.(Z_fit .+ col_med, 1e-10)

    σ_cells = has_proj ? resid_c[proj_mask] : vec(resid_c)
    σ       = max(1.4826 * median(abs.(σ_cells)), 1e-10)

    active = resid .> 0
    β_fit, γ_fit = fit_beta_gamma(Z, Z_fit_adj, σ, active; β_range=β_range, γ_range=γ_range)

    pg = prob_good(Z, Z_fit_adj, σ; β=β_fit, γ=γ_fit)
    pg[resid .≤ 0] .= 1.0
    mask = pg .> prob_threshold

    clean_mask = has_proj ? mask .& proj_mask : mask
    C_clean = fit_masked(Bx, By, Z, clean_mask; λ=λ)

    basis = SplineFlagBasis(p, kx, ky, C_clean, σ, β_fit, γ_fit, prob_threshold, nx, ny)

    return Z_fit, mask, basis
end

function flag_broadband_times_surface(
    Z        :: AbstractMatrix{<:Real};
    p        :: Int     = 3,
    λ        :: Float64 = 1e-6,
    σ_thresh :: Float64 = 3.0
) :: BitVector
    n_t, n_f = size(Z)
    n_t < 2*(p+1) && return falses(n_t)

    ε  = eps(Float64)
    ξ  = [clamp(xi, ε, 1-ε) for xi in range(0.0, 1.0, length=n_t)]
    η  = [clamp(yi, ε, 1-ε) for yi in range(0.0, 1.0, length=n_f)]
    kx = knot_vector(0.0, 1.0, p+1, p)
    ky = knot_vector(0.0, 1.0, p+1, p)
    Bx, By   = bases(kx, ky, p, ξ, η)
    C        = fit_ls(Bx, By, Z; λ=λ)
    Z_coarse = evaluate(Bx, By, C)

    resid   = Z .- Z_coarse
    row_med = [median(view(resid, t, :)) for t in 1:n_t]

    neg_r = row_med[row_med .< 0]
    σ     = isempty(neg_r) ? max(1.4826 * median(abs.(row_med)), 1e-10) :
                             max(1.4826 * median(abs.(neg_r)), 1e-10)
    row_med .> σ_thresh * σ
end

end
