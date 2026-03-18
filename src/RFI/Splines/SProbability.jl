module SProbability

using Statistics

export prob_good

function log_i0(x::Float64)
    if x < 3.75
        t = (x / 3.75)^2
        I0 = 1.0 + t * (3.5156229 + t * (3.0899424 + t * (1.2067492 +
                                                          t * (0.2659732 + t * (0.0360768 + t * 0.0045813)))))
        return log(I0)
    else
        return x - 0.5 * log(2π * x) + 1.0 / (8.0 * x)
    end
end

function prob_good(
    Z::AbstractMatrix,
    Z_fit::AbstractMatrix,
    σ::Float64;
    β::Float64=0.01,
    γ::Float64=10.0
)
    σ2 = σ^2
    γ2σ2 = (γ * σ)^2
    log_β_ratio = log(β / (1.0 - β))
    log_2γ = 2.0 * log(γ)
    p = similar(Z, Float64)

    for idx in eachindex(Z)
        z = Z[idx]
        ν = Z_fit[idx]

        if z <= 0.0 || ν <= 0.0
            p[idx] = 1.0
            continue
        end

        x_g = z * ν / σ2
        x_b = z * ν / γ2σ2
        log_fg = log_i0(x_g) - (z^2 + ν^2) / (2.0 * σ2)
        log_fb = log_i0(x_b) - (z^2 + ν^2) / (2.0 * γ2σ2) - log_2γ

        log_odds = log_β_ratio + log_fb - log_fg
        pij = 1.0 / (1.0 + exp(log_odds))
        p[idx] = isfinite(pij) ? clamp(pij, 1e-12, 1.0) :
                 (log_odds > 0.0 ? 0.0 : 1.0)
    end

    return p
end

end