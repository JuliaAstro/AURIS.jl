struct GainTable
    kind      :: Symbol # :B (bandpass), :G (gain)
    antennas  :: Vector{Int}
    spw_index :: Vector{Int}
    time_ns   :: Vector{Int64}
    gains     :: Array{ComplexF64,5}
    ok        :: Array{Bool,5}
end

function Base.show(io::IO, g::GainTable)
    nr, nc, ns, na, nt = size(g.gains)
    print(io, "GainTable(:", g.kind, ", ", nr, " rec × ", nc, " chan × ",
          ns, " spw × ", na, " ant × ", nt, " times, ",
          round(100 * count(g.ok) / length(g.ok); digits=1), "% ok)")
end

function product_receptors(npol::Int)
    npol == 4 && return [(1, 1), (1, 2), (2, 1), (2, 2)] # RR RL LR LL
    npol == 2 && return [(1, 1), (2, 2)] # RR LL
    npol == 1 && return [(1, 1)] # RR
    error("Unsupported polarization product count: $npol")
end

# Nearest solution slot in time with a valid solution,
# 0 if the antenna never solved.
function nearest_ok(g::GainTable, r::Int, ci::Int, si::Int, a::Int, t_ns::Int64)
    best, bestd = 0, typemax(Int64)
    for ti in eachindex(g.time_ns)
        g.ok[r, ci, si, a, ti] || continue
        d = abs(g.time_ns[ti] - t_ns)
        d < bestd && ((best, bestd) = (ti, d))
    end
    best
end

"""
    applycal!(v::VisibilityDataset, tables...) → v

Correct `v.vis` in place by each `GainTable` in turn:
`V_ij ← V_ij / (g_i(p) conj(g_j(q)))` per polarization product with receptors
(p,q); `:B` tables apply per channel, `:G` tables use the nearest *valid*
solution in time per antenna.  Samples with no valid solution at all are flagged.  
Weights are scaled by |g_i g_j|^2

"""
function applycal!(v::VisibilityDataset, tables::GainTable...)
    npol = Data.n_pol(v)
    recs = product_receptors(npol)
    nrec_v = maximum(maximum.(recs))
    for g in tables
        ant_slot = Dict(a => k for (k, a) in enumerate(g.antennas))
        spw_slot = Dict(s => k for (k, s) in enumerate(g.spw_index))
        nc_sol   = size(g.gains, 2)
        nt_sol   = length(g.time_ns)
        Threads.@threads for t in 1:Data.n_time(v)
            ti1 = Vector{Int}(undef, nrec_v)   
            ti2 = Vector{Int}(undef, nrec_v)
            for bl in 1:Data.n_baseline(v)
                i = get(ant_slot, v.antenna1[bl], 0)
                j = get(ant_slot, v.antenna2[bl], 0)
                for (s, sidx) in enumerate(v.spw_index)
                    si = get(spw_slot, sidx, 0)
                    if i == 0 || j == 0 || si == 0
                        v.flags[:, :, s, bl, t] .= true
                        continue
                    end
                    for r in 1:nrec_v
                        if nt_sol == 1
                            ti1[r] = ti2[r] = 1
                        else
                            ti1[r] = nearest_ok(g, r, 1, si, i, v.time_ns[t])
                            ti2[r] = nearest_ok(g, r, 1, si, j, v.time_ns[t])
                        end
                    end
                    wsum, wn = 0.0, 0
                    for c in 1:Data.n_chan(v)
                        ci = nc_sol == 1 ? 1 : c
                        for (p, (ri, rj)) in enumerate(recs)
                            ta, tb = ti1[ri], ti2[rj]
                            if ta == 0 || tb == 0 ||
                               !(g.ok[ri, ci, si, i, ta] && g.ok[rj, ci, si, j, tb])
                                v.flags[p, c, s, bl, t] = true
                                continue
                            end
                            gg = g.gains[ri, ci, si, i, ta] *
                                 conj(g.gains[rj, ci, si, j, tb])
                            z = ComplexF64(v.vis[p, c, s, bl, t]) / gg
                            zs = convert(eltype(v.vis), z)
                            if isfinite(real(zs)) && isfinite(imag(zs))
                                v.vis[p, c, s, bl, t] = zs
                            else
                                v.flags[p, c, s, bl, t] = true
                            end
                            if ri == rj   
                                wsum += abs2(gg); wn += 1
                            end
                        end
                    end
                    wn > 0 && (v.weights[s, bl, t] *= Float32(wsum / wn))
                end
            end
        end
    end
    v
end

function merge_spw(tables::GainTable...)
    isempty(tables) && error("merge_spw: no tables")
    length(tables) == 1 && return tables[1]
    t1 = tables[1]
    all(t.kind == t1.kind for t in tables) ||
        error("merge_spw: tables have different kinds")
    all(t.antennas == t1.antennas for t in tables) ||
        error("merge_spw: tables have different antenna axes")
    nrec, nchan = size(t1.gains, 1), size(t1.gains, 2)
    all(size(t.gains)[1:2] == (nrec, nchan) for t in tables) ||
        error("merge_spw: tables have different receptor/channel axes")
    spw_index = reduce(vcat, [t.spw_index for t in tables])
    allunique(spw_index) || error("merge_spw: spw sets overlap")

    if t1.kind == :B
        all(t -> length(t.time_ns) == 1, tables) ||
            error("merge_spw: :B tables must have a single time slot")
        return GainTable(:B, copy(t1.antennas), spw_index,
                         [round(Int64, mean(only(t.time_ns) for t in tables))],
                         cat((t.gains for t in tables)...; dims = 3),
                         cat((t.ok    for t in tables)...; dims = 3))
    end

    time_ns = sort!(unique!(reduce(vcat, [t.time_ns for t in tables])))
    tslot   = Dict(tm => k for (k, tm) in enumerate(time_ns))
    gains   = zeros(ComplexF64, nrec, nchan, length(spw_index),
                    length(t1.antennas), length(time_ns))
    ok      = fill(false, size(gains))
    s0 = 0
    for t in tables
        ss = s0 .+ (1:length(t.spw_index))
        for (k, tm) in enumerate(t.time_ns)
            gains[:, :, ss, :, tslot[tm]] = t.gains[:, :, :, :, k]
            ok[:, :, ss, :, tslot[tm]]    = t.ok[:, :, :, :, k]
        end
        s0 += length(t.spw_index)
    end
    GainTable(t1.kind, copy(t1.antennas), spw_index, time_ns, gains, ok)
end