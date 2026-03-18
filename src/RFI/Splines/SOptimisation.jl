module SOptimisation

using ..SBasis: knot_vector, bases, basis
using LinearAlgebra, Statistics

export select_knots

function select_knots(
    Z::AbstractMatrix,
    ξ::AbstractVector,
    η::AbstractVector,
    p::Int;
    k_range=nothing,
    k_range_freq=nothing,
    λ=1e-6,
    bic_penalty::Float64=1.0
)
    z_t = vec(median(Z, dims=2))
    z_f = vec(median(Z, dims=1))

    nkx = bic_1d(z_t, ξ, p; k_range=k_range, λ=λ, freq_axis=false, bic_penalty=bic_penalty)
    nky = bic_1d(z_f, η, p; k_range=k_range_freq, λ=λ, freq_axis=true, bic_penalty=bic_penalty)

    return nkx, nky
end

function bic_1d(z, t, p; k_range=nothing, λ=1e-6, freq_axis=false, bic_penalty::Float64=1.0)
    n = length(z)
    k_max = if freq_axis
        max(p + 1, min(n ÷ 4, 20))
    else
        max(p + 1, n ÷ 4)
    end
    rng = isnothing(k_range) ? (p+1:2:k_max) : k_range

    best_bic = Inf
    best_nk = first(rng)

    for nk in rng
        kv = knot_vector(0.0, 1.0, nk, p)
        B = [basis(kv, a, p, ti) for ti in t, a in 1:nk]
        c = (B'B + λ * I(nk)) \ (B'z)

        resid = z .- B * c
        med_r = median(resid)
        σ_mad = max(1.4826 * median(abs.(resid .- med_r)), 1e-10)
        inliers = abs.(resid .- med_r) .< 3.0 * σ_mad
        n_in = sum(inliers)
        n_in < 5 && continue

        rss_in = max(sum(resid[inliers] .^ 2), 1e-30)
        bic = n_in * log(rss_in / n_in) + bic_penalty * nk * log(n_in)

        if bic < best_bic
            best_bic = bic
            best_nk = nk
        end
    end

    return best_nk
end

end