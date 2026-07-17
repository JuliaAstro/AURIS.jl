using AURIS

args = copy(ARGS)

outdir = joinpath(@__DIR__, "pipeline_out")

i = findfirst(==("--outdir"), args)
if i !== nothing
    i == length(args) && error("--outdir requires a directory")
    outdir = args[i + 1]
    deleteat!(args, i:i+1)
end

sdm_path = length(args) >= 1 ? args[1] : get(ENV, "AURIS_SDM_PATH", "")
isempty(sdm_path) && error(
    "usage: julia --project=. $(PROGRAM_FILE) [--outdir DIR] <path/to/ASDM> [target_scan ...]"
)

tscans = length(args) >= 2 ? parse.(Int, args[2:end]) : nothing

result = run_pipeline(sdm_path; target_scans=tscans, outdir)

println()
println(pipeline_report(result))
println("Outputs in: $outdir")
