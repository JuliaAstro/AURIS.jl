const PERLEY_BUTLER_2017 = Dict{String,Vector{Float64}}(
    "3C48"  => [1.3253, -0.7553, -0.1914, 0.0498],
    "3C123" => [1.8017, -0.7884, -0.1035, -0.0248, 0.0090],
    "3C138" => [1.0088, -0.4981, -0.1552, -0.0102, 0.0223],
    "3C147" => [1.4516, -0.6961, -0.2007, 0.0640, -0.0464, 0.0289],
    "3C196" => [1.2872, -0.8530, -0.1534, -0.0200, 0.0201],
    "3C286" => [1.2481, -0.4507, -0.1798, 0.0357],
    "3C295" => [1.4701, -0.7658, -0.2780, -0.0347, 0.0399],
)

function setjy_flux(name::AbstractString, freq_ghz::Real)
    for (cal, a) in PERLEY_BUTLER_2017
        occursin(cal, name) || continue
        lf = log10(freq_ghz)
        return exp10(sum(a[i] * lf^(i - 1) for i in eachindex(a)))
    end
    NaN
end

function setjy_model(v::VisibilityDataset)
    S = [setjy_flux(v.field_name, f / 1e9) for f in v.freqs]
    if any(isnan, S)
        @warn "Field $(v.field_name) is not a standard flux calibrator; using a unit model"
        return ones(size(v.freqs))
    end
    S
end

function fluxscale(G_sec::GainTable, G_prim::GainTable)
    nrec, _, nspw, nant, ntime = size(G_sec.gains)
    ant_p = Dict(a => k for (k, a) in enumerate(G_prim.antennas))
    spw_p = Dict(s => k for (k, s) in enumerate(G_prim.spw_index))

    # Instrument gain amplitude from the primary: median over its times.
    flux = fill(NaN, nrec, nspw)
    for s in 1:nspw, r in 1:nrec
        sp = get(spw_p, G_sec.spw_index[s], 0)
        sp == 0 && continue
        ratios = Float64[]
        for a in 1:nant
            ap = get(ant_p, G_sec.antennas[a], 0)
            ap == 0 && continue
            gp = [abs(G_prim.gains[r, 1, sp, ap, t]) for t in 1:length(G_prim.time_ns)
                  if G_prim.ok[r, 1, sp, ap, t]]
            isempty(gp) && continue
            g0 = median(gp)
            g0 > 0 || continue
            for t in 1:ntime
                G_sec.ok[r, 1, s, a, t] || continue
                push!(ratios, abs2(G_sec.gains[r, 1, s, a, t]) / g0^2)
            end
        end
        isempty(ratios) || (flux[r, s] = median(ratios))
    end

    gains = copy(G_sec.gains)
    ok    = copy(G_sec.ok)
    for s in 1:nspw, r in 1:nrec
        if isnan(flux[r, s]) || flux[r, s] <= 0
            ok[r, :, s, :, :] .= false
        else
            gains[r, :, s, :, :] ./= sqrt(flux[r, s])
        end
    end
    (flux = flux,
     table = GainTable(G_sec.kind, copy(G_sec.antennas), copy(G_sec.spw_index),
                       copy(G_sec.time_ns), gains, ok))
end