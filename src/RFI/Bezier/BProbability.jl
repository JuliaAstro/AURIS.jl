module BProbability

export bayesian_prob_good

function bayesian_prob_good(resid, σ; β=0.01, γ=10.0)
    σ = max(σ, eps())
    R = resid ./ σ

    finite_mask = isfinite.(R)
    if !all(finite_mask)
        R = copy(R)
        R .= ifelse.(finite_mask, R, 0.0)
    end

    num = (1 - β) * exp.(-0.5 .* R .^ 2)
    den = num .+ (β / γ) * exp.(-0.5 .* (R ./ γ) .^ 2)
    p = num ./ den
    p .= ifelse.(isfinite.(p), p, 0.0)
    clamp.(p, 1e-12, 1.0)
end

end