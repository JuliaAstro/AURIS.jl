parallel_products(npol::Int) = npol >= 2 ? [1, npol] : [1]

function solint_groups(v::VisibilityDataset, solint_ns::Integer)
    solint_ns <= 0 && return [[t] for t in 1:Data.n_time(v)]
    gs = Vector{Vector{Int}}()
    t0 = v.time_ns[1]
    for t in 1:Data.n_time(v)
        if isempty(gs) || v.time_ns[t] - t0 >= solint_ns
            push!(gs, Int[])
            t0 = v.time_ns[t]
        end
        push!(gs[end], t)
    end
    gs
end

function solve_antennas(v::VisibilityDataset)
    ants = sort(unique(vcat(v.antenna1, v.antenna2)))
    slot = Dict(a => k for (k, a) in enumerate(ants))
    i1   = [slot[a] for a in v.antenna1]
    i2   = [slot[a] for a in v.antenna2]
    ants, i1, i2
end

model_at(model::Number, c, s) = Float64(model)
model_at(model::AbstractMatrix, c, s) = Float64(model[c, s])

function gain_ok(g::Vector{ComplexF64}, snr::Vector{Float64},
                 frac::Vector{Float64};
                 minsnr::Real=3.0, amp_window::Real=10.0, minfrac::Real=0.3)
    amps = [abs(x) for x in g if isfinite(abs(x)) && abs(x) > 0]
    isempty(amps) && return falses(length(g))
    med = median(amps)
    [isfinite(abs(g[a])) && med / amp_window < abs(g[a]) < med * amp_window &&
     snr[a] >= minsnr && frac[a] >= minfrac for a in eachindex(g)]
end

function amp_screen!(gains, ok; amp_window::Real=10.0)
    amps = [abs(gains[i]) for i in eachindex(ok)
            if ok[i] && isfinite(abs(gains[i]))]
    isempty(amps) && return
    med = median(amps)
    for i in eachindex(ok)
        ok[i] || continue
        a = abs(gains[i])
        if !(med / amp_window < a < med * amp_window)
            ok[i] = false
            gains[i] = complex(NaN, NaN)
        end
    end
end

function accumulate_matrix(v::VisibilityDataset, i1, i2, nant::Int,
                           p::Int, s::Int, cr, tr, model)
    Vs = zeros(ComplexF64, nant, nant)
    Ws = zeros(Float64, nant, nant)
    nused = zeros(Int, nant)
    act = fill(false, length(cr))
    nactive = 0
    for t in tr
        fill!(act, false)
        for bl in 1:Data.n_baseline(v)
            i, j = i1[bl], i2[bl]
            w = Float64(v.weights[s, bl, t])
            w > 0 || continue
            for (ci, c) in enumerate(cr)
                v.flags[p, c, s, bl, t] && continue
                m = model_at(model, c, s)
                m > 0 || continue
                wm = w * m^2
                z  = ComplexF64(v.vis[p, c, s, bl, t]) / m
                Vs[i, j] += wm * z
                Ws[i, j] += wm
                nused[i] += 1
                nused[j] += 1
                act[ci] = true
            end
        end
        nactive += count(act)
    end
    for j in 1:nant, i in 1:j-1     
        if Ws[i, j] > 0
            Vs[i, j] /= Ws[i, j]
            Vs[j, i] = conj(Vs[i, j])
            Ws[j, i] = Ws[i, j]
        end
    end
    npossible = (nant - 1) * nactive
    Vs, Ws, nused ./ max(npossible, 1)
end

"""
    solve_bandpass(v; refant=1, model=1.0) → GainTable

Per-channel antenna gains (`:B`) from a calibrator `VisibilityDataset`,
averaging over the full time range per channel Phases are referenced to `refant` 
(an index into the SDM antenna vector; defaults to the first antenna present).

The independent (spw, receptor, channel) solves run in parallel when Julia is
started with threads.
"""
function solve_bandpass(v::VisibilityDataset; refant::Int=0,
                        model::Union{Real,AbstractMatrix}=1.0,
                        minsnr::Real=15.0, amp_window::Real=10.0,
                        minfrac::Real=0.3)
    ants, i1, i2 = solve_antennas(v)
    nant  = length(ants)
    npol  = Data.n_pol(v)
    prods = parallel_products(npol)
    nrec  = length(prods)
    nchan, nspw = Data.n_chan(v), Data.n_spw(v)
    ref   = refant == 0 ? 1 : something(findfirst(==(refant), ants), 1)

    gains = fill(complex(NaN, NaN), nrec, nchan, nspw, nant, 1)
    ok   = fill(false, nrec, nchan, nspw, nant, 1)
    jobs = [(s, r, p, c) for s in 1:nspw for (r, p) in enumerate(prods)
                         for c in 1:nchan]
    Threads.@threads for job in jobs
        s, r, p, c = job
        Vs, Ws, frac = accumulate_matrix(v, i1, i2, nant, p, s, c:c,
                                         1:Data.n_time(v), model)
        g, conv, snr = stefcal(Vs, Ws; refant=ref)
        conv || continue
        good = gain_ok(g, snr, frac; minsnr, minfrac)
        for a in 1:nant
            if good[a]
                gains[r, c, s, a, 1] = g[a]
                ok[r, c, s, a, 1] = true
            end
        end
    end
    amp_screen!(gains, ok; amp_window)
    lo, hi = extrema(v.time_ns)
    mid = lo + (hi - lo) ÷ 2
    GainTable(:B, ants, copy(v.spw_index), [mid], gains, ok)
end

"""
    solve_gains(v; solint_ns=0, refant=1, model=1.0) → GainTable

Time-resolved antenna gains (`:G`) from a calibrator `VisibilityDataset`,
averaging over all channels per solution interval. `solint_ns = 0` solves per
integration otherwise integrations are grouped into intervals of `solint_ns`
Apply the bandpass first (`applycal!(v, B)`) so the channel
average is coherent.

The independent (interval, spw, receptor) solves run in parallel when Julia is
started with threads.
"""
function solve_gains(v::VisibilityDataset; solint_ns::Integer=0, refant::Int=0,
                     model::Union{Real,AbstractMatrix}=1.0, minsnr::Real=15.0,
                     amp_window::Real=3.0, minfrac::Real=0.3)
    ants, i1, i2 = solve_antennas(v)
    nant  = length(ants)
    npol  = Data.n_pol(v)
    prods = parallel_products(npol)
    nrec  = length(prods)
    nspw  = Data.n_spw(v)
    ref   = refant == 0 ? 1 : something(findfirst(==(refant), ants), 1)

    groups = solint_groups(v, solint_ns)
    nsol  = length(groups)
    gains = fill(complex(NaN, NaN), nrec, 1, nspw, nant, nsol)
    ok    = fill(false, nrec, 1, nspw, nant, nsol)
    tmid  = Vector{Int64}(undef, nsol)
    for (k, tr) in enumerate(groups)
        tmid[k] = v.time_ns[tr[1]] + (v.time_ns[tr[end]] - v.time_ns[tr[1]]) ÷ 2  
    end
    jobs = [(k, s, r, p) for k in 1:nsol for s in 1:nspw
                         for (r, p) in enumerate(prods)]
    Threads.@threads for job in jobs
        k, s, r, p = job
        Vs, Ws, frac = accumulate_matrix(v, i1, i2, nant, p, s,
                                         1:Data.n_chan(v), groups[k], model)
        g, conv, snr = stefcal(Vs, Ws; refant=ref)
        conv || continue
        good = gain_ok(g, snr, frac; minsnr, minfrac)
        for a in 1:nant
            if good[a]
                gains[r, 1, s, a, k] = g[a]
                ok[r, 1, s, a, k] = true
            end
        end
    end
    amp_screen!(gains, ok; amp_window)
    GainTable(:G, ants, copy(v.spw_index), tmid, gains, ok)
end

function accumulate_selfcal!(Vs, Ws, v::VisibilityDataset, model_vis,
                             i1, i2, p::Int, ss, tr)
    nant = size(Vs, 1)
    nused = zeros(Int, nant)
    act = fill(false, Data.n_chan(v))
    nactive = 0
    for s in ss, t in tr
        fill!(act, false)
        for bl in 1:Data.n_baseline(v)
            i, j = i1[bl], i2[bl]
            w = Float64(v.weights[s, bl, t])
            w > 0 || continue
            for c in 1:Data.n_chan(v)
                v.flags[p, c, s, bl, t] && continue
                m = ComplexF64(model_vis[c, s, bl, t])
                am2 = abs2(m)
                am2 > 0 || continue
                wm = w * am2
                Vs[i, j] += wm * (ComplexF64(v.vis[p, c, s, bl, t]) / m)
                Ws[i, j] += wm
                nused[i] += 1
                nused[j] += 1
                act[c] = true
            end
        end
        nactive += count(act)
    end
    npossible = (nant - 1) * nactive
    nused ./ max(npossible, 1)
end

"""
    solve_selfcal(v, model_vis; mode=:phase, solint_ns=0, refant=0, minsnr=3.0,
                  combine_spw=false) → GainTable

Self-calibration: time-resolved antenna gains (`:G`) solved against a sky
model rather than a point source. `model_vis` holds the model visibility of
every sample, shaped `(chan, spw, baseline, time)` produced from a CLEAN
component image by `Imaging.predict_dataset_vis`. Each sample is divided by
its model (weights scaled by |M|^2).
"""
function solve_selfcal(v::VisibilityDataset, model_vis::AbstractArray{<:Complex,4};
                       mode::Symbol=:phase, solint_ns::Integer=0, refant::Int=0,
                       minsnr::Real=3.0, combine_spw::Bool=false,
                       amp_window::Real=3.0, minfrac::Real=0.3)
    mode in (:phase, :amphase) ||
        throw(ArgumentError("mode must be :phase or :amphase"))
    size(model_vis) == (Data.n_chan(v), Data.n_spw(v), Data.n_baseline(v),
                        Data.n_time(v)) ||
        throw(DimensionMismatch("model_vis must be (chan, spw, baseline, time)"))

    ants, i1, i2 = solve_antennas(v)
    nant  = length(ants)
    npol  = Data.n_pol(v)
    prods = parallel_products(npol)
    nrec  = length(prods)
    nspw  = Data.n_spw(v)
    ref   = refant == 0 ? 1 : something(findfirst(==(refant), ants), 1)

    groups  = solint_groups(v, solint_ns)
    spwsets = combine_spw ? [1:nspw] : [s:s for s in 1:nspw]
    nsol  = length(groups)
    gains = fill(complex(NaN, NaN), nrec, 1, nspw, nant, nsol)
    ok    = falses(nrec, 1, nspw, nant, nsol)
    tmid  = Vector{Int64}(undef, nsol)
    Vs = Matrix{ComplexF64}(undef, nant, nant)
    Ws = Matrix{Float64}(undef, nant, nant)
    for (k, tr) in enumerate(groups)
        tmid[k] = v.time_ns[tr[1]] + (v.time_ns[tr[end]] - v.time_ns[tr[1]]) ÷ 2  
        for ss in spwsets, (r, p) in enumerate(prods)
            fill!(Vs, 0); fill!(Ws, 0)
            frac = accumulate_selfcal!(Vs, Ws, v, model_vis, i1, i2, p, ss, tr)
            for j in 1:nant, i in 1:j-1    
                if Ws[i, j] > 0
                    Vs[i, j] /= Ws[i, j]
                    Vs[j, i] = conj(Vs[i, j])
                    Ws[j, i] = Ws[i, j]
                end
            end
            g, conv, snr = stefcal(Vs, Ws; refant=ref)
            conv || continue
            good = gain_ok(g, snr, frac; minsnr, minfrac)
            for a in 1:nant
                good[a] || continue
                for s in ss
                    gains[r, 1, s, a, k] = g[a]
                    ok[r, 1, s, a, k] = true
                end
            end
        end
    end
    amp_screen!(gains, ok; amp_window)
    if mode == :phase
        for i in eachindex(gains)
            ok[i] && (gains[i] /= abs(gains[i]))
        end
    end
    GainTable(:G, ants, copy(v.spw_index), tmid, gains, collect(ok))
end