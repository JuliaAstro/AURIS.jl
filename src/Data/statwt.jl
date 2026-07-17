"""
    statwt!(v; minsamples=30) → v

Replace `v.weights` with the inverse variance 1/σ² measured from the data
itself, per (spw, baseline) — cf. CASA `statwt`.
"""
function statwt!(v::VisibilityDataset; minsamples::Int=30)
    npol  = n_pol(v)
    prods = npol >= 2 ? (1, npol) : (1,)
    d2 = Float64[]
    for s in 1:n_spw(v), bl in 1:n_baseline(v)
        empty!(d2)
        for t in 2:n_time(v), p in prods, c in 1:n_chan(v)
            (v.flags[p, c, s, bl, t] || v.flags[p, c, s, bl, t-1]) && continue
            δ = ComplexF64(v.vis[p, c, s, bl, t]) -
                ComplexF64(v.vis[p, c, s, bl, t-1])
            push!(d2, abs2(δ))
        end
        w = if length(d2) < minsamples
            0.0f0
        else
            σ2 = median(d2) / (2 * log(2))
            σ2 > 0 ? Float32(1 / σ2) : 0.0f0
        end
        v.weights[s, bl, :] .= w
    end
    v
end
