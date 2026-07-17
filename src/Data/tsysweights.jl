# Nearest-in-time value of a sorted series; 0.0 if the series is empty.
function nearest_value(ts::Vector{Int64}, vs::Vector{Float64}, t::Int64)
    isempty(ts) && return 0.0
    k = searchsortedfirst(ts, t)
    k <= 1 && return vs[1]
    k > length(ts) && return vs[end]
    (t - ts[k-1]) <= (ts[k] - t) ? vs[k-1] : vs[k]
end

"""
    tsys_weights!(v, series) → v

Sets `v.weights[s, bl, t] = Δν(s)·τ(t) / (Tsys_i(t)·Tsys_j(t))` 
using the nearest Tsys sample in time. Baselines without a series for 
both antennas get weight 0.
"""
function tsys_weights!(v::VisibilityDataset,
                       series::Dict{Tuple{Int,Int},Tuple{Vector{Int64},Vector{Float64}}})
    for si in 1:n_spw(v)
        s = v.spw_index[si]
        for bl in 1:n_baseline(v)
            ki = get(series, (v.antenna1[bl], s), nothing)
            kj = get(series, (v.antenna2[bl], s), nothing)
            for t in 1:n_time(v)
                w = 0.0
                if ki !== nothing && kj !== nothing
                    Ti = nearest_value(ki..., v.time_ns[t])
                    Tj = nearest_value(kj..., v.time_ns[t])
                    Ti > 0 && Tj > 0 &&
                        (w = v.chan_width[si] * v.interval[t] / (Ti * Tj))
                end
                v.weights[si, bl, t] = Float32(w)
            end
        end
    end
    v
end

"""
    tsys_weights!(v, ds; syspower=SDM.load_syspower(ds)) → v

Replace `v.weights` with switched-power Tsys weights (see file header):
`Tsys = Tcal·(Psum−Pdif)/(2·Pdif)` from the SysPower table, `Tcal` from
CalDevice `noiseCal` (1.0 K where absent, making that antenna/spw relative),
averaged over the receptors with a positive switched-power difference and
matched to data timestamps by nearest neighbour.
"""
function tsys_weights!(v::VisibilityDataset, ds::SDM.SDMDataset;
                       syspower::AbstractVector=SDM.load_syspower(ds))
    spw_by_id = Dict(s.raw["spectralWindowId"] => i for (i, s) in enumerate(ds.spws))
    ant_by_number = Dict(a.number => i for (i, a) in enumerate(ds.antennas))
    antidx(idstr) = get(ant_by_number, id_index(idstr), 0)

    tcal = Dict{Tuple{Int,Int},Float64}()
    for cdv in ds.calDevices
        a = antidx(cdv.antennaID)
        s = get(spw_by_id, cdv.spectralWindowID, 0)
        (a == 0 || s == 0) && continue
        nc = cdv.noiseCal
        nc isa AbstractVector && !isempty(nc) && isfinite(nc[1]) && nc[1] > 0 ||
            continue
        tcal[(a, s)] = Float64(nc[1])
    end

    series = Dict{Tuple{Int,Int},Tuple{Vector{Int64},Vector{Float64}}}()
    for r in syspower
        a = antidx(r.antennaID)
        s = get(spw_by_id, r.spwID, 0)
        (a == 0 || s == 0) && continue
        d, m = r.switchedDiff, r.switchedSum
        (d isa AbstractVector && m isa AbstractVector) || continue
        acc, n = 0.0, 0
        for k in eachindex(d)
            k <= length(m) || break
            df, sm = Float64(d[k]), Float64(m[k])
            (isfinite(df) && df > 0 && isfinite(sm) && sm > df) || continue
            acc += (sm - df) / (2 * df)
            n += 1
        end
        n == 0 && continue
        ts, vs = get!(series, (a, s), (Int64[], Float64[]))
        push!(ts, Int64(r.time))
        push!(vs, get(tcal, (a, s), 1.0) * acc / n)
    end
    for (ts, vs) in values(series)
        issorted(ts) || (p = sortperm(ts); permute!(ts, p); permute!(vs, p))
    end
    tsys_weights!(v, series)
end
