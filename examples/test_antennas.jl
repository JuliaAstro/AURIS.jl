using AURIS
using Plots
using Statistics: mean

# Add path to SDM here
ds = SDM.open_sdm("")
@info "Loaded SDM" telescope=ds.observation.telescopeName n_antenna=ds.n_antenna n_spw=ds.n_spw

antennas = ds.antennas
positions = [collect(Float64, a.position) for a in antennas]   
names = [a.name for a in antennas]

centroid = mean(positions)
lambda = atan(centroid[2], centroid[1])                     
phi = atan(centroid[3], hypot(centroid[1], centroid[2])) 

function east_north(p)
    d = p .- centroid
    east  = -sin(lambda)*d[1] + cos(lambda)*d[2]
    north = -sin(phi)*cos(lambda)*d[1] - sin(phi)*sin(lambda)*d[2] + cos(phi)*d[3]
    east, north
end

en = east_north.(positions)
east  = first.(en)
north = last.(en)


config = get(ds.observation.raw, "configName", "")
title  = "$(ds.observation.telescopeName)" * " ($config-array)" * " — $(ds.n_antenna) antennas"

gr()
plt = scatter(east, north;
    marker=(:circle, 6, :steelblue),
    label="",
    xlabel="East (m)", ylabel="North (m)",
    title=title, aspect_ratio=:equal, framestyle=:box, size=(760, 720))

span = maximum(north) - minimum(north)
for (e, n, nm) in zip(east, north, names)
    annotate!(plt, e, n + 0.02*span, text(nm, 7, :black, :bottom))
end

outfile = joinpath(@__DIR__, "test_antennas.png")
savefig(plt, outfile)
