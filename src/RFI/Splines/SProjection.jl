module SProjection

using LinearAlgebra, Statistics
using ..SBasis: knot_vector, basis

export flag_broadband_times, flag_narrowband_freqs

function fit_1d_bic(
    z       :: AbstractVector,
    ξ       :: AbstractVector,
    p       :: Int;
    k_range = nothing,
    λ       = 1e-6
) :: Tuple{Vector{Float64}, Float64}
    n = length(z)
    n < 2 * (p + 1) && return (zeros(n), 1.0)

    win  = max(2*(p+1), n ÷ 5)
    half = win ÷ 2
    bg   = similar(z)
    for i in 1:n
        lo = max(1, i - half); hi = min(n, i + half)
        bg[i] = median(z[lo:hi])
    end

    r0  = z .- bg
    med = median(r0)
    σ0  = max(1.4826 * median(abs.(r0 .- med)), 1e-10)
    w   = r0 .- med .< 3.0 * σ0
    sum(w) < 2 * (p + 1) && (w = trues(n))

    k_max    = max(p + 1, min(n ÷ 4, 20))
    rng      = isnothing(k_range) ? (p+1 : 2 : k_max) : k_range
    best_bic = Inf
    best_nk  = first(rng)

    for nk in rng
        kv  = knot_vector(0.0, 1.0, nk, p)
        B   = [basis(kv, a, p, ti) for ti in ξ, a in 1:nk]
        Bw  = B[w, :]
        c   = (Bw'Bw + λ * I(nk)) \ (Bw'z[w])
        r   = z[w] .- Bw * c
        med = median(r)
        σ   = max(1.4826 * median(abs.(r .- med)), 1e-10)
        nin = sum(abs.(r .- med) .< 3.0 * σ)
        nin < 5 && continue
        rss = max(sum(r[abs.(r .- med) .< 3.0 * σ].^2), 1e-30)
        bic = nin * log(rss / nin) + nk * log(nin)
        if bic < best_bic
            best_bic = bic
            best_nk  = nk
        end
    end

    kv    = knot_vector(0.0, 1.0, best_nk, p)
    B     = [basis(kv, a, p, ti) for ti in ξ, a in 1:best_nk]
    Bw    = B[w, :]
    c     = (Bw'Bw + λ * I(best_nk)) \ (Bw'z[w])
    z_fit = B * c
    r_w   = z[w] .- Bw * c
    med_w = median(r_w)
    σ_mad = max(1.4826 * median(abs.(r_w .- med_w)), 1e-10)
    return z_fit, σ_mad
end

function per_baseline_residuals(
    M       :: AbstractMatrix{<:Real},
    ξ       :: AbstractVector,
    p       :: Int;
    k_range = nothing,
    λ       = 1e-6
) :: Matrix{Float64}
    n_bl, n = size(M)
    δ = zeros(n_bl, n)
    for b in 1:n_bl
        z_fit, σ_b = fit_1d_bic(M[b, :], ξ, p; k_range=k_range, λ=λ)
        δ[b, :] = (M[b, :] .- z_fit) ./ σ_b
    end
    δ
end

function coherence_flag(
    δ_norm   :: AbstractMatrix{<:Real},
    σ_thresh :: Float64,
    min_rank :: Int
) :: BitVector
    n_bl, n = size(δ_norm)
    flags = falses(n)
    for i in 1:n
        flags[i] = sort(δ_norm[:, i])[min_rank] > σ_thresh
    end
    flags
end

function flag_broadband_times(
    Z_all    :: AbstractArray{<:Real, 3};
    p        :: Int     = 3,
    k_range             = nothing,
    λ                   = 1e-6,
    σ_thresh :: Float64 = 2.0,
    min_rank :: Int     = 2
) :: BitVector
    n_bl, n_t, n_f = size(Z_all)
    n_t < 2 * (p + 1) && return falses(n_t)
    M_bt = dropdims(median(Z_all, dims=3), dims=3)
    ξ    = collect(range(0.0, 1.0, length=n_t))
    δ    = per_baseline_residuals(M_bt, ξ, p; k_range=k_range, λ=λ)
    coherence_flag(δ, σ_thresh, min_rank)
end

function flag_narrowband_freqs(
    Z_all    :: AbstractArray{<:Real, 3};
    p        :: Int     = 3,
    k_range             = nothing,
    λ                   = 1e-6,
    σ_thresh :: Float64 = 2.0,
    min_rank :: Int     = 2,
    q_time   :: Float64 = 0.95
) :: BitVector
    n_bl, n_t, n_f = size(Z_all)
    n_f < 2 * (p + 1) && return falses(n_f)
    η  = collect(range(0.0, 1.0, length=n_f))
    δ  = zeros(n_bl, n_f)
    for b in 1:n_bl
        z_med      = [median(view(Z_all, b, :, f)) for f in 1:n_f]
        z_fit, _   = fit_1d_bic(z_med, η, p; k_range=k_range, λ=λ)
        r_tf       = zeros(n_t, n_f)
        for t in 1:n_t
            for f in 1:n_f
                r_tf[t, f] = Z_all[b, t, f] - z_fit[f]
            end
            med_t = median(view(r_tf, t, :))
            for f in 1:n_f
                r_tf[t, f] -= med_t
            end
            σ_t = max(1.4826 * median(abs.(view(r_tf, t, :))), 1e-10)
            for f in 1:n_f
                r_tf[t, f] /= σ_t
            end
        end
        for f in 1:n_f
            δ[b, f] = quantile(view(r_tf, :, f), q_time)
        end
    end
    coherence_flag(δ, σ_thresh, min_rank)
end

end
