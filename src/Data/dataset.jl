"""
    baseline_pairs(n) → Vector{Tuple{Int,Int}}

BDF triangular baseline ordering over `n` antenna slots:
(1,2), (1,3), (2,3), (1,4), (2,4), (3,4), …  (1-based slot indices, ant1 < ant2).
"""
function baseline_pairs(n::Integer)
    pairs = Vector{Tuple{Int,Int}}(undef, n * (n - 1) ÷ 2)
    k = 1
    for j in 2:n, i in 1:j-1
        pairs[k] = (i, j)
        k += 1
    end
    pairs
end

"""
    VisibilityDataset

Calibration/imaging-ready visibilities for one field and one correlator
configuration, concatenated in time over one or more scans.

Axes of `vis`/`flags`: `(pol, chan, spw, baseline, time)` (BDF-native order).
"""
struct VisibilityDataset{T<:Complex}
    field_id   :: String
    field_name :: String
    ra         :: Float64
    dec        :: Float64
    scan       :: Vector{Int}
    time_ns    :: Vector{Int64}
    jd_utc     :: Vector{Float64}
    interval   :: Vector{Float64}
    antenna1   :: Vector{Int}
    antenna2   :: Vector{Int}
    spw_index  :: Vector{Int}
    freqs      :: Matrix{Float64}
    chan_width :: Vector{Float64}
    vis        :: Array{T,5}
    flags      :: Array{Bool,5}
    weights    :: Array{Float32,3}
    uvw        :: Array{Float64,3}
end

n_time(v::VisibilityDataset)     = length(v.time_ns)
n_baseline(v::VisibilityDataset) = length(v.antenna1)
n_spw(v::VisibilityDataset)      = length(v.spw_index)
n_chan(v::VisibilityDataset)     = size(v.freqs, 1)
n_pol(v::VisibilityDataset)      = size(v.vis, 1)

function Base.show(io::IO, v::VisibilityDataset)
    print(io, "VisibilityDataset(", v.field_name, ": ",
          n_pol(v), " pol × ", n_chan(v), " chan × ", n_spw(v), " spw × ",
          n_baseline(v), " bl × ", n_time(v), " times, ",
          round(100 * count(v.flags) / length(v.flags); digits=1), "% flagged)")
end

# Selection 

scans_with_intent(ds::SDM.SDMDataset, intent::AbstractString) =
    filter(s -> any(occursin(intent, i) for i in s.scanIntent), ds.scans)

mains_for_scans(ds::SDM.SDMDataset, scan_numbers) =
    filter(m -> m.scanNumber in scan_numbers, ds.mains)

mains_for_field(ds::SDM.SDMDataset, field_id::AbstractString) =
    filter(m -> m.fieldID == field_id, ds.mains)

field_by_id(ds::SDM.SDMDataset, field_id::AbstractString) =
    ds.fields[findfirst(f -> f.raw["fieldId"] == field_id, ds.fields)]

# Trailing integer of an ASDM id, e.g. "SpectralWindow_7" → 7.
id_index(s::AbstractString) =
    (m = match(r"_(\d+)$", s); m === nothing ? 0 : parse(Int, m.captures[1]))

# ConfigDescription row referenced by a Main row.
function config_for(ds::SDM.SDMDataset, main::Main)
    i = findfirst(c -> c.configDescriptionID == main.configDescriptionID,
                  ds.configDescriptions)
    i === nothing && error("Main references unknown $(main.configDescriptionID)")
    ds.configDescriptions[i]
end

# Antenna slots: BDF antenna order is the ConfigDescription antennaId order.
# Returns indices into ds.antennas, one per slot.
function antenna_slots(ds::SDM.SDMDataset, cd)
    by_number = Dict(a.number => i for (i, a) in enumerate(ds.antennas))
    [by_number[id_index(aid)] for aid in cd.antennaID]
end

# BDF spw slots: the ConfigDescription dataDescriptionId order defines the spw
# axis.  Returns indices into ds.spws, one per slot.
function spw_slots(ds::SDM.SDMDataset, cd, n_spw::Int)
    dd_by_id  = Dict(d.raw["dataDescriptionId"] => d for d in ds.datadescs)
    spw_by_id = Dict(s.raw["spectralWindowId"] => i for (i, s) in enumerate(ds.spws))
    slots = [spw_by_id[dd_by_id[ddid].spwID] for ddid in cd.dataDescriptionID]
    length(slots) == n_spw ? slots :
        (@warn "ConfigDescription lists $(length(slots)) spws, BDF has $n_spw; using 1:$n_spw";
         collect(1:n_spw))
end

# Online flags (Flag.xml)
# Antenna numbers (id_index of antennaId) flagged by a Flag row.
flag_antenna_numbers(f::Flag) =
    haskey(f.raw, "antennaId") ? id_index.(SDM.asdm_array(f.raw["antennaId"])) : Int[]

function apply_online_flags!(flags::Array{Bool,5}, ds::SDM.SDMDataset,
                             slots::Vector{Int}, pairs::Vector{Tuple{Int,Int}},
                             time_ns::Vector{Int64}, interval_ns::Vector{Int64})
    slot_of = Dict(ds.antennas[s].number => k for (k, s) in enumerate(slots))
    nflagged = 0
    for f in ds.flags
        t0, t1 = f.time, f.time + f.interval
        bad_slots = Set(slot_of[n] for n in flag_antenna_numbers(f) if haskey(slot_of, n))
        isempty(bad_slots) && continue
        for (t, (tm, dt)) in enumerate(zip(time_ns, interval_ns))
            (tm + dt ÷ 2 > t0 && tm - dt ÷ 2 < t1) || continue
            for (k, (i, j)) in enumerate(pairs)
                if i in bad_slots || j in bad_slots
                    flags[:, :, :, k, t] .= true
                    nflagged += 1
                end
            end
        end
    end
    nflagged
end

"""
    load_visibility_dataset(ds, mains; dut1=0.0, eop=nothing) → VisibilityDataset
    load_visibility_dataset(ds, main; ...)
    load_visibility_dataset(ds; field=..., intent=..., scans=..., ...)

Load the BDF binaries for the given `Main` rows (all must share one field and
one ConfigDescription) into a `VisibilityDataset`: complex visibilities, flags
initialised from `Flag.xml` plus exact-zero samples, initial weights
Δν_chan × Δt, and UVWs from the full ITRF→J2000 rotation (see `uvw!`;
`dut1`/`eop` refine Earth orientation).  `T` selects the visibility storage
type (default `Complex{Float16}`.
"""
function load_visibility_dataset(ds::SDM.SDMDataset, mains::Vector{Main};
                                 T::Type{<:Complex}=Complex{Float16},
                                 dut1::Real=0.0, eop=nothing)
    isempty(mains) && error("No Main rows selected")
    allequal(m.fieldID for m in mains) ||
        error("Main rows span several fields; load them separately")
    allequal(m.configDescriptionID for m in mains) ||
        error("Main rows span several ConfigDescriptions; load them separately")

    field = field_by_id(ds, mains[1].fieldID)
    length(field.phaseDir) >= 2 || error("Field $(field.name) has no phase direction")
    ra, dec = field.phaseDir[1], field.phaseDir[2]

    cd    = config_for(ds, mains[1])
    slots = antenna_slots(ds, cd)

    # Per-integration metadata across all Main rows (one BDF file each).
    bdfs = [SDM.open_bdf(ds, m) for m in mains]
    try
        nbl, nspw, nchan, npol =
            bdfs[1].n_baseline, bdfs[1].n_spw, bdfs[1].n_chan, bdfs[1].n_cross_pol
        length(slots) == bdfs[1].n_antenna ||
            error("ConfigDescription lists $(length(slots)) antennas, BDF has $(bdfs[1].n_antenna)")
        spws  = spw_slots(ds, cd, nspw)
        pairs = baseline_pairs(length(slots))

        freqs = Matrix{Float64}(undef, nchan, nspw)
        chanw = Vector{Float64}(undef, nspw)
        for (k, si) in enumerate(spws)
            sw = ds.spws[si]
            sw.numChan == nchan ||
                error("spw $(sw.name) has $(sw.numChan) channels, BDF has $nchan")
            freqs[:, k] = sw.chanFreq
            chanw[k]    = abs(sw.chanWidth[1])
        end

        ntot        = sum(length, bdfs)
        scan        = Vector{Int}(undef, ntot)
        time_ns     = Vector{Int64}(undef, ntot)
        interval_ns = Vector{Int64}(undef, ntot)
        vis         = Array{T,5}(undef, npol, nchan, nspw, nbl, ntot)

        t = 0
        for (m, bdf) in zip(mains, bdfs)
            for k in 1:length(bdf)
                t += 1
                meta           = bdf.integrations[k]
                scan[t]        = m.scanNumber
                time_ns[t]     = meta.time_ns
                interval_ns[t] = meta.interval_ns
                copyto!(view(vis, :, :, :, :, t), SDM.cross_view(bdf, k))
            end
        end

        # Flags: exact-zero samples (unwritten/dropped data) + online Flag.xml.
        flags = map(iszero, vis)
        apply_online_flags!(flags, ds, slots, pairs, time_ns, interval_ns)

        # Initial weights ∝ 1/σ²: channel width × integration time.
        interval_s = interval_ns .* 1e-9
        weights = Array{Float32,3}(undef, nspw, nbl, ntot)
        for t in 1:ntot, bl in 1:nbl, s in 1:nspw
            weights[s, bl, t] = Float32(chanw[s] * interval_s[t])
        end

        # UVWs: one ITRF→J2000 rotation per integration, projected per baseline.
        positions = reduce(hcat, ds.antennas[s].position for s in slots)  # 3×n_ant
        B = Matrix{Float64}(undef, 3, nbl)
        for (k, (i, j)) in enumerate(pairs)
            B[:, k] = positions[:, j] .- positions[:, i]   # r(ant2) − r(ant1)
        end
        jd_utc = tai_ns_to_jd_utc.(time_ns)
        uvw = Array{Float64,3}(undef, 3, nbl, ntot)
        for t in 1:ntot
            uvw!(view(uvw, :, :, t), B, jd_utc[t], ra, dec; dut1, eop)
        end

        ant1 = [slots[i] for (i, _) in pairs]
        ant2 = [slots[j] for (_, j) in pairs]

        return VisibilityDataset(
            mains[1].fieldID, field.name, ra, dec,
            scan, time_ns, jd_utc, interval_s,
            ant1, ant2, spws, freqs, chanw,
            vis, flags, weights, uvw)
    finally
        foreach(close, bdfs)
    end
end

load_visibility_dataset(ds::SDM.SDMDataset, main::Main; kwargs...) =
    load_visibility_dataset(ds, [main]; kwargs...)

function load_visibility_dataset(ds::SDM.SDMDataset;
                                 field=nothing, intent=nothing, scans=nothing,
                                 kwargs...)
    mains = ds.mains
    field  === nothing || (mains = filter(m -> m.fieldID == field, mains))
    scans  === nothing || (mains = filter(m -> m.scanNumber in scans, mains))
    if intent !== nothing
        ok = Set(s.scanNumber for s in scans_with_intent(ds, intent))
        mains = filter(m -> m.scanNumber in ok, mains)
    end
    load_visibility_dataset(ds, mains; kwargs...)
end
