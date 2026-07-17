function rank_refants(flagfrac::AbstractVector{<:Real},
                      positions::AbstractMatrix{<:Real})
    n = length(flagfrac)
    size(positions, 2) == n ||
        throw(DimensionMismatch("positions must be 3 × $(n)"))
    centroid = vec(mean(positions; dims=2))
    d = [norm(positions[:, a] .- centroid) for a in 1:n]
    d0 = max(median(d), eps())
    score = [(1 - flagfrac[a]) / (1 + d[a] / d0) for a in 1:n]
    sortperm(score; rev=true)
end

function rank_refants(v::VisibilityDataset, ds::SDM.SDMDataset)
    ants = sort(unique(vcat(v.antenna1, v.antenna2)))
    frac = map(ants) do a
        bad, tot = 0, 0
        for bl in 1:n_baseline(v)
            (v.antenna1[bl] == a || v.antenna2[bl] == a) || continue
            f = view(v.flags, :, :, :, bl, :)
            bad += count(f)
            tot += length(f)
        end
        tot == 0 ? 1.0 : bad / tot
    end
    pos = reduce(hcat, ds.antennas[a].position for a in ants)
    ants[rank_refants(frac, pos)]
end
