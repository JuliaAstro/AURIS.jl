module Components

const components_dir = joinpath(@__DIR__, "Components")

include(joinpath(components_dir, "antenna.jl"))
include(joinpath(components_dir, "geometry.jl"))
include(joinpath(components_dir, "receiver.jl"))
include(joinpath(components_dir, "atmosphere.jl"))
include(joinpath(components_dir, "caldata.jl"))
include(joinpath(components_dir, "caldevice.jl"))
include(joinpath(components_dir, "calreduction.jl"))
include(joinpath(components_dir, "configdescription.jl"))
include(joinpath(components_dir, "correlatormode.jl"))
include(joinpath(components_dir, "datadescrip.jl"))
include(joinpath(components_dir, "feed.jl"))
include(joinpath(components_dir, "field.jl"))
include(joinpath(components_dir, "flags.jl"))
include(joinpath(components_dir, "history.jl"))
include(joinpath(components_dir, "main.jl"))
include(joinpath(components_dir, "observation.jl"))
include(joinpath(components_dir, "pointing.jl"))
include(joinpath(components_dir, "pointingmodel.jl"))
include(joinpath(components_dir, "polarization.jl"))
include(joinpath(components_dir, "processor.jl"))
include(joinpath(components_dir, "sbsummary.jl"))
include(joinpath(components_dir, "scan.jl"))
include(joinpath(components_dir, "source.jl"))
include(joinpath(components_dir, "specwindow.jl"))
include(joinpath(components_dir, "state.jl"))
include(joinpath(components_dir, "station.jl"))
include(joinpath(components_dir, "subscan.jl"))
include(joinpath(components_dir, "switchcycle.jl"))
include(joinpath(components_dir, "syscal.jl"))
include(joinpath(components_dir, "syspower.jl"))
include(joinpath(components_dir, "weather.jl"))

export Antenna, Atmosphere, CalData, CalDevice, CalReduction, Circular,
       ConfigDescription, CorrelatorMode, DataDescrip, Feed, Field, Flag,
       Geometry, History, Linear, Main, Observation, Pointing, PointingModel,
       Polarization, Processor, Receiver, SBSummary, Scan, Source, specWindow,
       State, Station, Subscan, SwitchCycle, SysCal, SysPower, Weather

end
