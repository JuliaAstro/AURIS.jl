module SProbability

using Statistics

export prob_good, fit_beta_gamma

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

function fit_beta_gamma(
    Z::AbstractMatrix,
    Z_fit::AbstractMatrix,
    σ::Float64,
    active::AbstractMatrix{Bool};
    β_range=exp.(range(log(0.001), log(0.5), length=40)),
    γ_range=exp.(range(log(2.0), log(5000.0), length=40))
)
    σ2 = σ^2

    log_fg = Float64[]
    for idx in eachindex(Z)
        active[idx] || continue
        z = Float64(Z[idx])
        ν = Float64(Z_fit[idx])
        push!(log_fg, log_i0(z * ν / σ2) - (z^2 + ν^2) / (2.0 * σ2) - 2.0 * log(σ))
    end
    isempty(log_fg) && return (β_range[1], γ_range[1])

    best_ll = -Inf
    β_opt = Float64(first(β_range))
    γ_opt = Float64(first(γ_range))

    for γ in γ_range
        γ = Float64(γ)
        γ2σ2 = (γ * σ)^2
        log_γσ2 = 2.0 * log(γ * σ)

        log_fb = Float64[]
        k = 0
        for idx in eachindex(Z)
            active[idx] || continue
            k += 1
            z = Float64(Z[idx])
            ν = Float64(Z_fit[idx])
            push!(log_fb, log_i0(z * ν / γ2σ2) - (z^2 + ν^2) / (2.0 * γ2σ2) - log_γσ2)
        end

        for β in β_range
            β = Float64(β)
            log1mβ = log(1.0 - β)
            logβ = log(β)
            ll = 0.0
            for k in eachindex(log_fg)
                a = log1mβ + log_fg[k]
                b = logβ + log_fb[k]

                if a >= b
                    ll += a + log1p(exp(b - a))
                else
                    ll += b + log1p(exp(a - b))
                end
            end
            if ll > best_ll
                best_ll = ll
                β_opt = β
                γ_opt = γ
            end
        end
    end

    return β_opt, γ_opt
end

end
