"""
    stefcal(V, W; model=1.0+0im, refant=1, maxiter=100, tol=1e-10)
        → (gains, converged, snr)

Solve `V[p,q] ≈ g[p] * conj(g[q]) * model` for the antenna gains.

- `V` — n×n visibility matrix, `V[p,q] = ⟨E_p E_q*⟩`
- `W` — n×n non-negative weights
- `model` — point-source model visibility 
- `refant` — phases are referenced so `angle(g[refant]) = 0`
"""
function stefcal(V::AbstractMatrix, W::AbstractMatrix;
                 model::Number=1.0 + 0im, refant::Int=1,
                 maxiter::Int=100, tol::Real=1e-10)
    n = size(V, 1)
    size(V) == (n, n) && size(W) == (n, n) ||
        throw(DimensionMismatch("V and W must be square and equal-sized"))

    present = [any(q -> q != p && W[p, q] > 0, 1:n) for p in 1:n]
    g = ComplexF64[present[p] ? 1.0 : NaN for p in 1:n]
    any(present) || return g, false, zeros(n)
    gprev = copy(g)
    converged = false

    for it in 1:maxiter
        for p in 1:n
            present[p] || continue
            num = 0.0 + 0.0im
            den = 0.0
            for q in 1:n
                (q == p && continue)
                present[q] || continue
                w = W[p, q]
                w > 0 || continue
                c = conj(gprev[q]) * model      
                num += w * conj(c) * V[p, q]
                den += w * abs2(c)
            end
            g[p] = den > 0 ? num / den : gprev[p]
        end
        if iseven(it)
            δ = sqrt(sum(abs2(g[p] - gprev[p]) for p in 1:n if present[p]) /
                     max(sum(abs2(g[p]) for p in 1:n if present[p]), eps()))
            @. g = 0.5 * (g + gprev)
            if δ < tol
                converged = true
                break
            end
        end
        copyto!(gprev, g)
    end

    r = present[refant] && isfinite(abs(g[refant])) && abs(g[refant]) > 0 ?
        refant : findfirst(present)
    if r !== nothing && abs(g[r]) > 0
        ph = conj(g[r]) / abs(g[r])
        for p in 1:n
            present[p] && (g[p] *= ph)
        end
    end

    snr = zeros(n)
    for p in 1:n
        present[p] || continue
        rss, info = 0.0, 0.0
        nobs = 0
        for q in 1:n
            (q == p || !present[q]) && continue
            w = W[p, q]
            w > 0 || continue
            c = conj(g[q]) * model
            rss  += w * abs2(V[p, q] - g[p] * c)
            info += w * abs2(c)
            nobs += 1
        end
        info > 0 || continue
        ν = max(nobs - 1, 1)
        snr[p] = rss > 0 ? abs(g[p]) / sqrt(rss / (ν * info)) : Inf
    end
    g, converged, snr
end