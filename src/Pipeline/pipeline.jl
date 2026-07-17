struct StageLog
    name    :: String
    elapsed :: Float64
    info    :: Vector{Pair{String,String}}
end

struct TargetImage
    field_name :: String
    scans      :: Vector{Int}
    cell       :: Float64
    beam       :: NTuple{3,Float64}
    peak       :: Float64
    rms        :: Float64
    niter      :: Int
    model_flux :: Float64
    restored   :: Matrix{Float64}
    residual   :: Matrix{Float64}
    model      :: Matrix{Float64}
    path       :: String
end

struct PipelineResult
    refant      :: Int
    refant_name :: String
    ranking     :: Vector{Int}
    bandpass    :: GainTable
    gains_flux  :: GainTable
    gains_phase :: Dict{String,GainTable}
    flux        :: Dict{String,NamedTuple{(:value, :scatter),NTuple{2,Float64}}}
    targets     :: Vector{TargetImage}
    stages      :: Vector{StageLog}
end

okfrac(g::GainTable) = count(g.ok) / max(length(g.ok), 1)
med_amp(g::GainTable) = (a = abs.(g.gains[g.ok]); isempty(a) ? NaN : median(a))
fmt(x; digits=2) = !isfinite(x) ? "n/a" :
    digits == 0 ? string(round(Int, x)) : string(round(x; digits))

sanitize(name::AbstractString) =
    strip(replace(name, r"[^A-Za-z0-9+\-.]+" => "_"), '_')

function band_label(ds::SDM.SDMDataset, cid)
    fallback = "config" * string(SDM.id_index(cid))
    ci = findfirst(c -> c.configDescriptionID == cid, ds.configDescriptions)
    ci === nothing && return fallback
    ddids = ds.configDescriptions[ci].dataDescriptionID
    isempty(ddids) && return fallback
    di = findfirst(d -> d.raw["dataDescriptionId"] == ddids[1], ds.datadescs)
    di === nothing && return fallback
    si = findfirst(s -> s.raw["spectralWindowId"] == ds.datadescs[di].spwID, ds.spws)
    si === nothing && return fallback
    m = match(r"^EVLA_([^#]+)#", string(ds.spws[si].name))
    m === nothing ? fallback : String(m.captures[1])
end

resid_rms(res) = 1.4826 * median(abs.(res))

run_pipeline(sdm_path::AbstractString; kwargs...) =
    run_pipeline(SDM.open_sdm(sdm_path); kwargs...)

function run_pipeline(ds::SDM.SDMDataset;
                      target_scans = nothing,
                      refant::Int = 0,
                      solint_ns::Integer = 30_000_000_000,
                      npix::Int = 512,
                      cellfrac::Real = 4,
                      robust = nothing,
                      niter::Int = 3000,
                      nmajor::Int = 8,
                      major_frac::Real = 0.2,
                      do_statwt::Bool = true,
                      selfcal::Bool = false,
                      combine_spw::Bool = true,
                      outdir = nothing,
                      verbose::Bool = true)
    stages = StageLog[]
    function stage!(name, t0, info::Vector{Pair{String,String}})
        push!(stages, StageLog(name, time() - t0, info))
        verbose && @info join(["pipeline: $name ($(round(time() - t0; digits=1)) s)";
                               ["  $k: $v" for (k, v) in info]], "\n")
    end
    outdir === nothing || mkpath(outdir)

    # Scan selection by intent 
    t0 = time()
    scanset(intent) = Set(s.scanNumber for s in scans_with_intent(ds, intent))
    flux_scans, phase_scans, tgt_scans =
        scanset("CALIBRATE_FLUX"), scanset("CALIBRATE_PHASE"), scanset("OBSERVE_TARGET")
    have_bdf(m) = isfile(SDM.bdf_path(ds, m))

    m_flux_all = [m for m in ds.mains if m.scanNumber in flux_scans]
    isempty(m_flux_all) && error("No scan has CALIBRATE_FLUX intent")
    if !any(have_bdf, m_flux_all)
        bdir = joinpath(ds.path, "ASDMBinary")
        error("No CALIBRATE_FLUX scan has BDF data: e.g. scan " *
              "$(m_flux_all[1].scanNumber) expects $(SDM.bdf_path(ds, m_flux_all[1]))" *
              (isdir(bdir) ? " ($(length(readdir(bdir))) entries in ASDMBinary)" :
                             " — $bdir does not exist"))
    end
    filter!(have_bdf, m_flux_all)

    cal_configs = sort!(unique([m.configDescriptionID for m in m_flux_all]);
                        by = SDM.id_index)
    calset    = Set(cal_configs)
    multiband = length(cal_configs) > 1
    band      = Dict(c => band_label(ds, c) for c in cal_configs)
    tag(name, cid) = multiband ? "$name ($(band[cid]))" : name

    function usable_integrations(m)
        w0, w1 = m.time - m.interval ÷ 2, m.time + m.interval ÷ 2
        lost = 0.0
        for f in ds.flags
            ov = min(w1, f.time + f.interval) - max(w0, f.time)
            ov > 0 && (lost += ov * length(Data.flag_antenna_numbers(f)))
        end
        span = Float64(m.interval) * length(ds.antennas)
        m.numIntegration * max(1.0 - lost / max(span, 1.0), 0.0)
    end

    m_flux = Dict(cid => (cands = [m for m in m_flux_all
                                   if m.configDescriptionID == cid];
                          cands[argmax(usable_integrations.(cands))])
                  for cid in cal_configs)

    wanted(m) = target_scans === nothing || m.scanNumber in target_scans
    m_phase = [m for m in ds.mains if m.scanNumber in phase_scans &&
               !(m.scanNumber in flux_scans) && have_bdf(m)]
    m_target = [m for m in ds.mains if m.scanNumber in tgt_scans &&
                have_bdf(m) && wanted(m)]
    n_uncal = count(m -> !(m.configDescriptionID in calset),
                    Iterators.flatten((m_phase, m_target)))
    n_uncal > 0 &&
        @warn "Dropping $n_uncal scans whose ConfigDescription has no flux calibrator"
    filter!(m -> m.configDescriptionID in calset, m_phase)
    filter!(m -> m.configDescriptionID in calset, m_target)

    pgroups = Dict{Tuple{String,String},Vector{CMain}}()   
    for m in m_phase
        push!(get!(pgroups, (m.fieldID, m.configDescriptionID), CMain[]), m)
    end
    tgroups = Dict{Tuple{String,String},Vector{CMain}}()
    for m in m_target
        push!(get!(tgroups, (m.fieldID, m.configDescriptionID), CMain[]), m)
    end
    isempty(tgroups) && error("No OBSERVE_TARGET scans selected")

    fname(fid) = field_by_id(ds, fid).name
    info = Pair{String,String}[
        "flux scan" => join(["$(m_flux[c].scanNumber) (" *
                             tag(fname(m_flux[c].fieldID), c) * ")"
                             for c in cal_configs], ", "),
        "phase calibrators" =>
            join(sort!(unique([tag(fname(f), c) for (f, c) in keys(pgroups)])), ", "),
        "target fields" =>
            join(sort!(unique([tag(fname(f), c) for (f, c) in keys(tgroups)])), ", "),
        "target scans" => join(sort!([m.scanNumber for m in m_target]), ", ")]
    multiband && pushfirst!(info,
        "bands" => join([band[c] for c in cal_configs], ", "))
    stage!("scan selection", t0, info)

    # Flux calibrator: flag, refant, bandpass, Jy-scale gains per config
    Bs  = Dict{String,GainTable}()
    Gfs = Dict{String,GainTable}()
    auto_ref = refant == 0
    refants  = Dict{String,Int}()
    rankings = Dict{String,Vector{Int}}()
    for cid in cal_configs
        t0 = time()
        vf = load_visibility_dataset(ds, [m_flux[cid]])
        flag_edges!(vf)
        rfi_f = flag_rfi!(vf)
        rankings[cid] = auto_ref ? rank_refants(vf, ds) : Int[refant]
        ref = refants[cid] = rankings[cid][1]
        refinfo = auto_ref ?
            "$(ds.antennas[ref].name) (ranked: " *
            join([ds.antennas[a].name for a in rankings[cid][1:min(end, 5)]],
                 " > ") * " …)" :
            "$(ds.antennas[ref].name) (user-specified)"
        Sf = setjy_model(vf)
        B  = solve_bandpass(vf; model = Sf, refant = ref)
        applycal!(vf, B)
        Gf = solve_gains(vf; solint_ns, model = Sf, refant = ref)
        Bs[cid], Gfs[cid] = B, Gf
        stage!(tag("flux calibrator", cid), t0, [
            "field" => vf.field_name,
            "RFI flagged" => fmt(100 * rfi_f) * " %",
            "refant" => refinfo,
            "model flux" => fmt(minimum(Sf)) * "–" * fmt(maximum(Sf)) * " Jy (Perley–Butler 2017)",
            "bandpass ok" => fmt(100 * okfrac(B); digits=1) * " %",
            "gain ok" => fmt(100 * okfrac(Gf); digits=1) * " %",
            "gain amp (≈1 expected)" => fmt(med_amp(Gf); digits=3)])
    end

    deadset = Set(cid for cid in cal_configs if count(Bs[cid].ok) == 0)
    for cid in deadset
        @warn "No valid bandpass solutions for band $(band[cid]); dropping its scans"
    end
    if !isempty(deadset)
        filter!(c -> !(c in deadset), cal_configs)
        isempty(cal_configs) && error("No band produced valid bandpass solutions")
        filter!(p -> !(p.first[2] in deadset), pgroups)
        filter!(p -> !(p.first[2] in deadset), tgroups)
        isempty(tgroups) &&
            error("No OBSERVE_TARGET scans left after dropping uncalibratable bands")
    end
    refant  = refants[cal_configs[1]]      
    ranking = rankings[cal_configs[1]]
    refname = ds.antennas[refant].name
    B  = merge_spw((Bs[c]  for c in cal_configs)...)
    Gf = merge_spw((Gfs[c] for c in cal_configs)...)

    # Phase calibrators: gains + fluxscale bootstrap per field × config
    phase_tables = Dict{Tuple{String,String},GainTable}()   
    gains_phase  = Dict{String,GainTable}()                 
    flux = Dict{String,NamedTuple{(:value, :scatter),NTuple{2,Float64}}}()
    for ((fid, cid), ms) in sort!(collect(pgroups);
                                  by = p -> (SDM.id_index(p.first[2]),
                                             minimum(m.time for m in p.second)))
        t0 = time()
        vp = load_visibility_dataset(ds, ms)
        flag_edges!(vp)
        rfi_p = flag_rfi!(vp)
        applycal!(vp, Bs[cid])
        Gp = solve_gains(vp; solint_ns, refant = refants[cid])
        fs = fluxscale(Gp, Gfs[cid])
        good = filter(x -> isfinite(x) && x > 0, vec(fs.flux))
        S  = isempty(good) ? NaN : median(good)
        dS = isempty(good) ? NaN : 1.4826 * median(abs.(good .- S))
        pname = tag(vp.field_name, cid)
        phase_tables[(fid, cid)] = fs.table
        gains_phase[pname] = fs.table
        flux[pname] = (value = S, scatter = dS)
        stage!("phase calibrator $pname", t0, [
            "scans" => join(sort!([m.scanNumber for m in ms]), ", "),
            "RFI flagged" => fmt(100 * rfi_p) * " %",
            "gain ok" => fmt(100 * okfrac(fs.table); digits=1) * " %",
            "bootstrapped flux" => fmt(S; digits=3) * " ± " * fmt(dS; digits=3) * " Jy"])
    end

    function table_for(ms)
        cid  = ms[1].configDescriptionID
        keys_c = [k for k in keys(phase_tables) if k[2] == cid]
        isempty(keys_c) && return Gfs[cid]
        tmid = median([m.time for m in ms])
        _, key = findmin(Dict(k => minimum(abs(m.time - tmid) for m in pgroups[k])
                              for k in keys_c))
        phase_tables[key]
    end

    # Targets: apply, flag, weight, CLEAN, optional self-cal, FITS
    targets = TargetImage[]
    for ((fid, cid), ms) in sort!(collect(tgroups);
                                  by = p -> minimum(m.time for m in p.second))
        t0 = time()
        vt = load_visibility_dataset(ds, ms)
        tname = tag(vt.field_name, cid)
        flag_edges!(vt)
        applycal!(vt, Bs[cid], table_for(ms))
        rfi_t = flag_rfi!(vt)
        if all(vt.flags)
            @warn "All data flagged for target $tname after calibration; skipping"
            continue
        end
        do_statwt && statwt!(vt)

        function image()
            uv, vis, wt = uv_samples(vt)
            cell = nyquist_cell(uv) / cellfrac
            robust === nothing || (wt = briggs_weights(uv, wt; npix, cell, robust))
            r = cs_clean(uv, vis, wt; npix, cell, gain = 0.1, niter,
                         nmajor, major_frac)
            r, cell
        end
        r, cell = image()
        info = Pair{String,String}["scans" => join(sort!([m.scanNumber for m in ms]), ", "),
                                   "RFI flagged" => fmt(100 * rfi_t) * " %",
                                   "weighting" => (do_statwt ? "statwt" : "natural") *
                                                  (robust === nothing ? "" : ", Briggs robust=$(robust)")]
        if selfcal && sum(r.model) > 0
            mv = predict_dataset_vis(vt, r.model, cell)
            Gs = solve_selfcal(vt, mv; mode = :phase, solint_ns,
                               refant = refants[cid], combine_spw)
            if count(Gs.ok) > 0
                applycal!(vt, Gs)
                dr0 = maximum(r.restored) / resid_rms(r.residual)
                r, cell = image()
                ph = median(abs.(rad2deg.(angle.(Gs.gains[Gs.ok]))))
                push!(info, "self-cal" => "phase-only, ok " *
                      fmt(100 * okfrac(Gs); digits=1) * " %, median |φ| " *
                      fmt(ph; digits=1) * "°, DR ×" *
                      fmt(maximum(r.restored) / resid_rms(r.residual) / dr0))
            else
                push!(info, "self-cal" => "no valid solutions; kept transferred calibration")
            end
        end

        peak = maximum(r.restored)
        rms  = resid_rms(r.residual)
        path = ""
        if outdir !== nothing
            path = joinpath(outdir, sanitize(tname) * ".image.fits")
            f0, f1 = extrema(vt.freqs)
            write_fits(path, r.restored;
                       ra = vt.ra, dec = vt.dec, cell,
                       freq = (f0 + f1) / 2, bandwidth = f1 - f0 + median(vt.chan_width),
                       beam = r.beam, object = vt.field_name,
                       date_obs = Dates.format(Dates.julian2datetime(vt.jd_utc[1]),
                                               dateformat"yyyy-mm-ddTHH:MM:SS.sss"))
        end
        casec = cell * 180 / π * 3600
        append!(info, ["CLEAN" => "$(r.iters) components, model flux " *
                                  fmt(sum(r.model); digits=3) * " Jy",
                       "beam" => fmt(r.beam[1] * 206265; digits=1) * "″ × " *
                                 fmt(r.beam[2] * 206265; digits=1) * "″ (cell " *
                                 fmt(casec; digits=1) * "″)",
                       "peak / rms / DR" => fmt(1e3 * peak) * " mJy / " *
                                            fmt(1e3 * rms) * " mJy / " *
                                            fmt(peak / rms; digits=0)])
        isempty(path) || push!(info, "FITS" => path)
        stage!("target $tname", t0, info)
        push!(targets, TargetImage(tname,
                                   sort!([m.scanNumber for m in ms]), cell,
                                   r.beam, peak, rms, r.iters, sum(r.model),
                                   r.restored, r.residual, r.model, path))
    end

    result = PipelineResult(refant, refname, ranking, B, Gf, gains_phase, flux,
                            targets, stages)

    if outdir !== nothing
        anames = [a.name for a in ds.antennas]
        export_gains(joinpath(outdir, "bandpass.csv"), B; names = anames)
        export_gains(joinpath(outdir, "gains_flux.csv"), Gf; names = anames)
        for (name, g) in gains_phase
            export_gains(joinpath(outdir, "gains_" * sanitize(name) * ".csv"), g;
                         names = anames)
        end
        export_summary(joinpath(outdir, "summary.csv"), result)
        open(joinpath(outdir, "report.md"), "w") do io
            print(io, pipeline_report(result;
                dataset = get(ds.observation.raw, "execBlockUID", "")))
        end
        verbose && @info "pipeline: outputs written" outdir
    end
    result
end
