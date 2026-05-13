module SFourier

using LinearAlgebra, Statistics

using ..SBasis: knot_vector, basis

export rdft, rdft_inv, flag_broadband_times_fourier

function rdft(z::AbstractVector{<:Real})
    n   = length(z)
    n_k = n ÷ 2 + 1
    F_re = zeros(n_k)
    F_im = zeros(n_k)
    @inbounds for k in 0:n_k-1
        ωk = 2π * k / n
        for t in 1:n
            φ = ωk * (t - 1)
            F_re[k+1] += z[t] * cos(φ)
            F_im[k+1] -= z[t] * sin(φ)
        end
    end
    F_re, F_im
end

function rdft_inv(F_re::AbstractVector, F_im::AbstractVector, n::Int)
    n_k = length(F_re)
    z   = zeros(n)
    @inbounds for t in 1:n
        z[t] = F_re[1]
        for k in 1:n_k-1
            ωk  = 2π * k / n
            φ   = ωk * (t - 1)
            fac = (n % 2 == 0 && k == n_k - 1) ? 1.0 : 2.0
            z[t] += fac * (F_re[k+1] * cos(φ) - F_im[k+1] * sin(φ))
        end
        z[t] /= n
    end
    z
end

function fit_1d_robust(
    z       :: AbstractVector{<:Real},
    ξ       :: AbstractVector{<:Real},
    p       :: Int;
    k_range         = nothing,
    λ       :: Float64 = 1e-6
) :: Tuple{Vector{Float64}, Float64}
    n = length(z)
    n < 2*(p+1) && return (fill(median(z), n), max(1.4826*median(abs.(z.-median(z))),1e-10))

    win  = max(2*(p+1), n ÷ 5)
    half = win ÷ 2
    bg   = similar(z, Float64)
    for i in 1:n
        lo, hi = max(1, i - half), min(n, i + half)
        bg[i] = median(view(z, lo:hi))
    end

    r0  = z .- bg
    med = median(r0)
    σ0  = max(1.4826 * median(abs.(r0 .- med)), 1e-10)
    w   = r0 .- med .< 2.5 * σ0
    sum(w) < 2*(p+1) && (w = trues(n))

    k_max    = max(p+1, min(n ÷ 4, 20))
    rng      = isnothing(k_range) ? (p+1 : 2 : k_max) : k_range
    best_bic = Inf
    best_nk  = first(rng)

    for nk in rng
        kv  = knot_vector(0.0, 1.0, nk, p)
        B   = [basis(kv, a, p, ηi) for ηi in ξ, a in 1:nk]
        Bw  = B[w, :]
        c   = (Bw'Bw + λ * I(nk)) \ (Bw' * z[w])
        r   = z[w] .- Bw * c
        med_r = median(r)
        σ_r   = max(1.4826 * median(abs.(r .- med_r)), 1e-10)
        nin   = sum(abs.(r .- med_r) .< 3.0 * σ_r)
        nin < 5 && continue
        rss   = max(sum(r[abs.(r .- med_r) .< 3.0 * σ_r] .^ 2), 1e-30)
        bic   = nin * log(rss / nin) + nk * log(nin)
        if bic < best_bic
            best_bic = bic
            best_nk  = nk
        end
    end

    kv    = knot_vector(0.0, 1.0, best_nk, p)
    B     = [basis(kv, a, p, ηi) for ηi in ξ, a in 1:best_nk]
    Bw    = B[w, :]
    c     = (Bw'Bw + λ * I(best_nk)) \ (Bw' * z[w])
    z_fit = B * c
    r_w   = z[w] .- Bw * c
    σ     = max(1.4826 * median(abs.(r_w .- median(r_w))), 1e-10)
    z_fit, σ
end

function flag_broadband_times_fourier(
    Z_all    :: AbstractArray{<:Real, 3};
    p        :: Int     = 3,
    λ        :: Float64 = 1e-6,
    k_range             = nothing,
    σ_thresh :: Float64 = 3.0
) :: BitVector
    n_bl, n_t, n_f = size(Z_all)
    n_t < 2*(p+1) && return falses(n_t)

    z = [median(vec(view(Z_all, :, t, :))) for t in 1:n_t]

    ξ = collect(range(0.0, 1.0, length=n_t))
    z_fit, _ = fit_1d_robust(z, ξ, p; k_range=k_range, λ=λ)
    r = z .- z_fit

    neg_r = r[r .< 0]
    σ     = isempty(neg_r) ? max(1.4826 * median(abs.(r)), 1e-10) :
                             max(1.4826 * median(abs.(neg_r)), 1e-10)

    r .> σ_thresh * σ
end

end
