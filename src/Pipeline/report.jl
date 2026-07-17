
function pipeline_report(result::PipelineResult; dataset::AbstractString="")
    io = IOBuffer()
    println(io, "# AURIS pipeline run report")
    println(io)
    isempty(dataset) || println(io, "- dataset: `", dataset, "`")
    println(io, "- date: ", Dates.format(Dates.now(), dateformat"yyyy-mm-dd HH:MM"))
    println(io, "- reference antenna: **", result.refant_name, "**")
    println(io, "- total time: ",
            round(sum(s.elapsed for s in result.stages); digits=1), " s")
    println(io)

    println(io, "## Stages")
    for s in result.stages
        println(io)
        println(io, "### ", s.name, "  (", round(s.elapsed; digits=1), " s)")
        println(io)
        for (k, v) in s.info
            println(io, "- ", k, ": ", v)
        end
    end
    println(io)

    println(io, "| quantity | AURIS value |")
    println(io, "|---|---|")
    println(io, "| bandpass solutions ok | ",
            fmt(100 * okfrac(result.bandpass); digits=1), " % |")
    println(io, "| flux-cal gain amplitude (≈1) | ",
            fmt(med_amp(result.gains_flux); digits=3), " |")
    for (name, f) in sort!(collect(result.flux); by=first)
        println(io, "| ", name, " bootstrapped flux | ",
                fmt(f.value; digits=3), " ± ", fmt(f.scatter; digits=3), " Jy |")
    end
    for t in result.targets
        println(io, "| ", t.field_name, " peak / rms / DR | ",
                fmt(1e3 * t.peak), " mJy / ", fmt(1e3 * t.rms), " mJy / ",
                fmt(t.peak / t.rms; digits=0), " |")
        println(io, "| ", t.field_name, " restoring beam | ",
                fmt(t.beam[1] * 206265; digits=1), "″ × ",
                fmt(t.beam[2] * 206265; digits=1), "″, PA ",
                fmt(t.beam[3] * 180 / π; digits=0), "° |")
    end
    println(io)
    String(take!(io))
end

function export_gains(path::AbstractString, g::GainTable; names=nothing)
    nr, nc, ns, na, nt = size(g.gains)
    open(path, "w") do io
        println(io, "kind,receptor,chan,spw,antenna,antenna_name,time_ns,amp,phase_deg,ok")
        for t in 1:nt, a in 1:na, s in 1:ns, c in 1:nc, r in 1:nr
            z = g.gains[r, c, s, a, t]
            println(io, g.kind, ",", r, ",", c, ",", g.spw_index[s], ",",
                    g.antennas[a], ",",
                    names === nothing ? "" : names[g.antennas[a]], ",",
                    g.time_ns[t], ",",
                    abs(z), ",", rad2deg(angle(z)), ",",
                    g.ok[r, c, s, a, t] ? 1 : 0)
        end
    end
    path
end

function export_summary(path::AbstractString, result::PipelineResult)
    open(path, "w") do io
        println(io, "quantity,value")
        for (name, f) in sort!(collect(result.flux); by=first)
            println(io, "flux:", name, ",", f.value)
        end
        for t in result.targets
            println(io, "peak:", t.field_name, ",", t.peak)
            println(io, "rms:", t.field_name, ",", t.rms)
            println(io, "dr:", t.field_name, ",", t.peak / t.rms)
            println(io, "beam_maj_arcsec:", t.field_name, ",", t.beam[1] * 206265)
            println(io, "beam_min_arcsec:", t.field_name, ",", t.beam[2] * 206265)
            println(io, "beam_pa_deg:", t.field_name, ",", t.beam[3] * 180 / π)
        end
    end
    path
end
