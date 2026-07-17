# SDM XML tables → Components structs

struct SDMDataset
    path               :: String
    observation        :: Observation                
    mains              :: Vector{Main}                 
    scans              :: Vector{Scan}                 
    subscans           :: Vector{Subscan}              
    antennas           :: Vector{Antenna}              
    stations           :: Vector{Station}              
    spws               :: Vector{specWindow}           
    polarizations      :: Vector{Polarization}         
    datadescs          :: Vector{DataDescrip}          
    fields             :: Vector{Field}               
    sources            :: Vector{Source}               
    configDescriptions :: Vector{ConfigDescription}    
    correlatorModes    :: Vector{CorrelatorMode}       
    pointingModels     :: Vector{PointingModel}        
    sbSummaries        :: Vector{SBSummary}            
    switchCycles       :: Vector{SwitchCycle}          
    calData            :: Vector{CalData}              
    calReductions      :: Vector{CalReduction}         
    calDevices         :: Vector{CalDevice}            
    feeds              :: Vector{Feed}                
    flags              :: Vector{Flag}                 
    processors         :: Vector{Processor}            
    states             :: Vector{State}                
    receivers          :: Vector{Receiver}             
    weather            :: Vector{Weather}              
    n_antenna          :: Int
    n_spw              :: Int
    n_chan             :: Int                          
    n_cross_pol        :: Int                          
    n_auto_pol         :: Int                          
end

# Timestamp conversion 
const MJD_TO_UNIX_NS = Int64(40587) * Int64(86400) * Int64(1_000_000_000)
asdm_time_to_unix(ns::Int64) = (ns - MJD_TO_UNIX_NS) * 1e-9

# UID → BDF path
# uid:///evla/bdf/1444569572563  →  uid____evla_bdf_1444569572563
uid_to_filename(uid::String) = replace(uid, r"[^a-zA-Z0-9]" => "_")

function bdf_path(dataset::SDMDataset, main::Main)
    joinpath(dataset.path, "ASDMBinary", uid_to_filename(String(main.dataUID)))
end

# XML helpers 
function xml_values(text::AbstractString, tag::AbstractString)
    pat = Regex("<$(tag)>([^<]*)</$(tag)>")
    String[String(m.captures[1]) for m in eachmatch(pat, text)]
end

function xml_first(text::AbstractString, tag::AbstractString, default::String="")
    vs = xml_values(text, tag)
    isempty(vs) ? default : vs[1]
end

function xml_entityref(text::AbstractString, tag::AbstractString)
    m = match(Regex("(?s)<$(tag)>.*?entityId=\"([^\"]+)\""), text)
    m === nothing ? "" : String(m.captures[1])
end

# Attempt to get everything, in case something is missing from test SDM
function xml_all(row::AbstractString)
    d = Dict{String,String}()
    for m in eachmatch(r"<([A-Za-z][A-Za-z0-9]*)>([^<]*)</\1>", row)
        d[String(m.captures[1])] = strip(String(m.captures[2]))
    end
    for m in eachmatch(r"(?s)<([A-Za-z][A-Za-z0-9]*)>\s*<EntityRef entityId=\"([^\"]+)\"", row)
        d[String(m.captures[1])] = String(m.captures[2])
    end
    d
end

function asdm_array(s::AbstractString)
    parts = split(strip(s))
    isempty(parts) && return String[]
    ndims = parse(Int, parts[1])
    ndims == 1 || error("Only 1-D ASDM arrays handled here (got ndims=$ndims)")
    n = parse(Int, parts[2])
    String.(parts[3:2+n])
end

function asdm_floats(s::AbstractString)
    parts = split(strip(s))
    isempty(parts) && return Float64[]
    ndims = parse(Int, parts[1])
    dims  = parse.(Int, parts[2:1+ndims])   
    n     = prod(dims)
    parse.(Float64, parts[2+ndims:1+ndims+n])
end

function asdm_floats_all(s::AbstractString)
    parts = split(strip(s))
    isempty(parts) && return Float64[]
    ndims = parse(Int, parts[1])
    start = 2 + ndims
    start > length(parts) ? Float64[] : parse.(Float64, parts[start:end])
end

function asdm_interval(s::AbstractString)
    p = split(strip(s))
    t = length(p) >= 1 ? parse(Int64, p[1]) : Int64(0)
    d = length(p) >= 2 ? parse(Int64, p[2]) : Int64(0)
    t, d
end

xml_opt(text::AbstractString, tag::AbstractString) =
    occursin("<$tag>", text) ? xml_first(text, tag) : nothing

xml_rows(text::AbstractString) =
    String[String(m.match) for m in eachmatch(r"(?s)<row>(.*?)</row>", text)]

function table_rows(path::String, file::String)
    fp = joinpath(path, file)
    isfile(fp) ? xml_rows(read(fp, String)) : String[]
end

to_float(s, default=0.0) = (v = tryparse(Float64, s); v === nothing ? default : v)
to_int(s, default=0)     = (v = tryparse(Int, s);     v === nothing ? default : v)
to_int64(s, default=0)   = parse(Int64, isempty(s) ? string(default) : s)

# ExecBlock → Components.Observation 
function parse_observation(path::String)
    rows = table_rows(path, "ExecBlock.xml")
    isempty(rows) && error("No rows in ExecBlock.xml")
    r = rows[1]

    obs = Observation(
        nothing,                              
        xml_first(r, "observingLog"),         
        xml_first(r, "observerName"),         
        xml_entityref(r, "projectUID"),       
        nothing,                              
        nothing,                              
        xml_first(r, "schedulerMode"),        
        xml_first(r, "telescopeName"),        
        (to_int64(xml_first(r, "startTime")), to_int64(xml_first(r, "endTime"))),  
        xml_all(r),                           
    )
    nant = to_int(xml_first(r, "numAntenna", "0"))
    obs, nant
end

# Main.xml → Components.Main 
function parse_mains(path::String)
    mains = Main[]
    for r in table_rows(path, "Main.xml")
        push!(mains, Main(
            to_int64(xml_first(r, "time")),                
            to_int(xml_first(r, "numAntenna", "0")),      
            xml_first(r, "timeSampling"),              
            to_int64(xml_first(r, "interval")),            
            to_int(xml_first(r, "numIntegration", "0")),  
            to_int(xml_first(r, "scanNumber", "0")),      
            to_int(xml_first(r, "subscanNumber", "0")),   
            to_int64(xml_first(r, "dataSize")),            
            xml_entityref(r, "dataUID"),               
            xml_first(r, "configDescriptionId"),       
            xml_first(r, "execBlockId"),               
            xml_first(r, "fieldId"),                   
            asdm_array(xml_first(r, "stateId")),       
            xml_all(r),                                
        ))
    end
    sort!(mains, by = m -> (m.scanNumber, m.subscanNumber))
    mains
end

# Scan.xml → Components.Scan 
function parse_scans(path::String)
    scans = Scan[]
    for r in table_rows(path, "Scan.xml")
        push!(scans, Scan(
            to_int(xml_first(r, "scanNumber", "0")),
            to_int64(xml_first(r, "startTime")),
            to_int64(xml_first(r, "endTime")),
            to_int(xml_first(r, "numIntent", "0")),
            to_int(xml_first(r, "numSubscan", "0")),
            asdm_array(xml_first(r, "scanIntent")),
            asdm_array(xml_first(r, "calDataType")),
            asdm_array(xml_first(r, "calibrationOnLine")),
            xml_first(r, "sourceName"),
            xml_first(r, "execBlockId"),
            xml_all(r),                                
        ))
    end
    sort!(scans, by = s -> s.scanNumber)
    scans
end

# Subscan.xml → Components.Subscan 
function parse_subscans(path::String)
    subs = Subscan[]
    for r in table_rows(path, "Subscan.xml")
        push!(subs, Subscan(
            to_int(xml_first(r, "scanNumber", "0")),
            to_int(xml_first(r, "subscanNumber", "0")),
            to_int64(xml_first(r, "startTime")),
            to_int64(xml_first(r, "endTime")),
            xml_first(r, "fieldName"),
            xml_first(r, "subscanIntent"),
            to_int(xml_first(r, "numIntegration", "0")),
            to_int.(asdm_array(xml_first(r, "numSubintegration"))),
            xml_first(r, "execBlockId"),
            xml_all(r),                                
        ))
    end
    subs
end

# Station.xml → Components.Station 
function parse_stations(path::String)
    stations = Station[]
    for r in table_rows(path, "Station.xml")
        push!(stations, Station(
            xml_first(r, "stationId"),
            xml_first(r, "name"),
            asdm_floats(xml_first(r, "position")),   
            xml_first(r, "type"),
            xml_all(r),                              
        ))
    end
    stations
end

# SpectralWindow.xml → Components.specWindow 
function parse_spectral_windows(path::String)
    spws = specWindow[]
    for r in table_rows(path, "SpectralWindow.xml")
        nc   = to_int(xml_first(r, "numChan", "0"))
        fs   = to_float(xml_first(r, "chanFreqStart", "0"))
        step = to_float(xml_first(r, "chanFreqStep", "0"))
        cw   = to_float(xml_first(r, "chanWidth", "0"))
        chanFreq = fs .+ (0:nc-1) .* step   
        push!(spws, specWindow(
            nothing,                                  
            nothing,                                  
            xml_first(r, "basebandName"),             
            collect(chanFreq),                        
            fill(cw, nc),                             
            to_float(xml_first(r, "effectiveBw", "0")),    
            nothing,                                  
            nothing,                                  
            nothing,                                  
            nothing,                                  
            nothing,                                  
            xml_first(r, "name"),                     
            xml_first(r, "netSideband"),              
            nc,                                       
            to_float(xml_first(r, "refFreq", "0")),        
            to_float(xml_first(r, "resolution", "0")),     
            xml_first(r, "correlationBit"),           
            nothing,                                  
            xml_first(r, "windowFunction"),           
            to_float(xml_first(r, "totBandwidth", "0")),   
            xml_all(r),                               
        ))
    end
    spws
end

# Polarization.xml → Components.Polarization 
function parse_polarizations(path::String)
    pols = Polarization[]
    for r in table_rows(path, "Polarization.xml")
        push!(pols, Polarization(
            split(strip(xml_first(r, "corrProduct"))),  
            asdm_array(xml_first(r, "corrType")),        
            nothing,                                     
            to_int(xml_first(r, "numCorr", "0")),           
            xml_all(r),                                  
        ))
    end
    pols
end

# DataDescription.xml → Components.DataDescrip 
function parse_datadescs(path::String)
    dds = DataDescrip[]
    for r in table_rows(path, "DataDescription.xml")
        push!(dds, DataDescrip(
            nothing,                              
            xml_first(r, "polOrHoloId"),          
            xml_first(r, "spectralWindowId"),     
            xml_all(r),                           
        ))
    end
    dds
end

# Field.xml → Components.Field 
function parse_fields(path::String)
    fields = Field[]
    for r in table_rows(path, "Field.xml")
        push!(fields, Field(
            xml_first(r, "code"),                        
            asdm_floats(xml_first(r, "delayDir")),       
            nothing,                                     
            nothing,                                     
            nothing,                                     
            xml_first(r, "fieldName"),                   
            to_int(xml_first(r, "numPoly", "0")),           
            asdm_floats(xml_first(r, "phaseDir")),       
            nothing,                                     
            asdm_floats(xml_first(r, "referenceDir")),   
            nothing,                                     
            xml_first(r, "sourceId"),                    
            to_int64(xml_first(r, "time")),                  
            xml_all(r),                                  
        ))
    end
    fields
end

# Source.xml → Components.Source 
function parse_sources(path::String)
    sources = Source[]
    for r in table_rows(path, "Source.xml")
        ti = split(strip(xml_first(r, "timeInterval")))   
        t_start = length(ti) >= 1 ? parse(Int64, ti[1]) : Int64(0)
        t_dur   = length(ti) >= 2 ? parse(Int64, ti[2]) : Int64(0)
        push!(sources, Source(
            nothing,                                       
            xml_first(r, "code"),                          
            asdm_floats(xml_first(r, "direction")),        
            t_dur,                                         
            xml_first(r, "sourceName"),                    
            to_int(xml_first(r, "numLines", "0")),            
            nothing,                                       
            asdm_floats(xml_first(r, "properMotion")),     
            asdm_floats(xml_first(r, "restFrequency")),    
            xml_first(r, "sourceId"),                      
            xml_first(r, "spectralWindowId"),              
            asdm_floats(xml_first(r, "sysVel")),           
            t_start,                                       
            nothing,                                       
            xml_all(r),                                     
        ))
    end
    sources
end

#  Antenna.xml + Station.xml → Components.Antenna
id_index(s::AbstractString) = (m = match(r"_(\d+)$", s); m === nothing ? 0 : parse(Int, m.captures[1]))

function parse_antennas(path::String, stations::Vector{Station})
    sta_pos = Dict{String, Vector{Float64}}(s.stationID => s.position for s in stations)

    ants = Antenna[]
    for r in table_rows(path, "Antenna.xml")
        sid    = xml_first(r, "stationId")
        pos    = get(sta_pos, sid, Float64[0.0, 0.0, 0.0])
        offset = asdm_floats(xml_first(r, "offset"))
        diam   = to_float(xml_first(r, "dishDiameter", "0"))
        push!(ants, Antenna(
            id_index(xml_first(r, "antennaId")),   
            xml_first(r, "name"),                   
            sid,                                    
            pos,                                    
            offset,                                 
            nothing,                                
            xml_first(r, "antennaType"),            
            π * (diam / 2)^2,                       
            xml_all(r),                             
        ))
    end
    ants
end

# ConfigDescription.xml → Components.ConfigDescription 
function parse_configdescriptions(path::String)
    cds = ConfigDescription[]
    for r in table_rows(path, "ConfigDescription.xml")
        push!(cds, ConfigDescription(
            xml_first(r, "configDescriptionId"),
            to_int(xml_first(r, "numAntenna", "0")),
            to_int(xml_first(r, "numDataDescription", "0")),
            to_int(xml_first(r, "numFeed", "0")),
            xml_first(r, "correlationMode"),
            to_int(xml_first(r, "numAtmPhaseCorrection", "0")),
            asdm_array(xml_first(r, "atmPhaseCorrection")),
            xml_first(r, "processorType"),
            xml_first(r, "spectralType"),
            asdm_array(xml_first(r, "antennaId")),
            asdm_array(xml_first(r, "dataDescriptionId")),
            asdm_array(xml_first(r, "feedId")),
            xml_first(r, "processorId"),
            asdm_array(xml_first(r, "switchCycleId")),
            xml_all(r),                           # raw
        ))
    end
    cds
end

# CorrelatorMode.xml → Components.CorrelatorMode 
function parse_correlatormodes(path::String)
    cms = CorrelatorMode[]
    for r in table_rows(path, "CorrelatorMode.xml")
        push!(cms, CorrelatorMode(
            xml_first(r, "correlatorModeId"),
            to_int(xml_first(r, "numBaseband", "0")),
            asdm_array(xml_first(r, "basebandNames")),
            asdm_array(xml_first(r, "basebandConfig")),
            xml_first(r, "accumMode"),
            to_int(xml_first(r, "binMode", "0")),
            to_int(xml_first(r, "numAxes", "0")),
            asdm_array(xml_first(r, "axesOrderArray")),
            asdm_array(xml_first(r, "filterMode")),
            xml_first(r, "correlatorName"),
            xml_all(r),                           
        ))
    end
    cms
end

# PointingModel.xml → Components.PointingModel 
function parse_pointingmodels(path::String)
    pms = PointingModel[]
    for r in table_rows(path, "PointingModel.xml")
        push!(pms, PointingModel(
            xml_first(r, "pointingModelId"),
            to_int(xml_first(r, "numCoeff", "0")),
            asdm_array(xml_first(r, "coeffName")),
            asdm_floats(xml_first(r, "coeffVal")),
            xml_first(r, "polarizationType"),
            xml_first(r, "receiverBand"),
            xml_first(r, "assocNature"),
            xml_first(r, "antennaId"),
            xml_first(r, "assocPointingModelId"),
            xml_all(r),                           
        ))
    end
    pms
end

# SBSummary.xml → Components.SBSummary 
function parse_sbsummaries(path::String)
    sbs = SBSummary[]
    for r in table_rows(path, "SBSummary.xml")
        push!(sbs, SBSummary(
            xml_first(r, "sBSummaryId"),
            xml_entityref(r, "sbSummaryUID"),
            xml_entityref(r, "projectUID"),
            xml_entityref(r, "obsUnitSetUID"),
            to_float(xml_first(r, "frequency", "0")),
            xml_first(r, "frequencyBand"),
            xml_first(r, "sbType"),
            to_int64(xml_first(r, "sbDuration")),
            to_int(xml_first(r, "numObservingMode", "0")),
            asdm_array(xml_first(r, "observingMode")),
            to_int(xml_first(r, "numberRepeats", "0")),
            to_int(xml_first(r, "numScienceGoal", "0")),
            asdm_array(xml_first(r, "scienceGoal")),
            to_int(xml_first(r, "numWeatherConstraint", "0")),
            asdm_array(xml_first(r, "weatherConstraint")),
            asdm_floats(xml_first(r, "centerDirection")),
            xml_all(r),                           
        ))
    end
    sbs
end

# SwitchCycle.xml → Components.SwitchCycle 
function parse_switchcycles(path::String)
    scs = SwitchCycle[]
    for r in table_rows(path, "SwitchCycle.xml")
        push!(scs, SwitchCycle(
            xml_first(r, "switchCycleId"),
            to_int(xml_first(r, "numStep", "0")),
            asdm_floats(xml_first(r, "weightArray")),
            asdm_floats(xml_first(r, "dirOffsetArray")),
            asdm_floats(xml_first(r, "freqOffsetArray")),
            asdm_array(xml_first(r, "stepDurationArray")),   
            xml_all(r),                                      
        ))
    end
    scs
end

# CalData.xml → Components.CalData
function parse_caldata(path::String)
    cds = CalData[]
    for r in table_rows(path, "CalData.xml")
        push!(cds, CalData(
            xml_first(r, "calDataId"),
            to_int64(xml_first(r, "startTimeObserved")),
            to_int64(xml_first(r, "endTimeObserved")),
            xml_entityref(r, "execBlockUID"),
            xml_first(r, "calDataType"),
            xml_first(r, "calType"),
            to_int(xml_first(r, "numScan", "0")),
            to_int.(asdm_array(xml_first(r, "scanSet"))),
            xml_all(r),                          
        ))
    end
    cds
end

# CalReduction.xml → Components.CalReduction 
function parse_calreductions(path::String)
    crs = CalReduction[]
    for r in table_rows(path, "CalReduction.xml")
        push!(crs, CalReduction(
            xml_first(r, "calReductionId"),
            to_int(xml_first(r, "numApplied", "0")),
            asdm_array(xml_first(r, "appliedCalibrations")),
            to_int(xml_first(r, "numParam", "0")),
            xml_first(r, "paramSet"),        
            to_int(xml_first(r, "numInvalidConditions", "0")),
            asdm_array(xml_first(r, "invalidConditions")),
            to_int64(xml_first(r, "timeReduced")),
            xml_first(r, "messages"),
            xml_first(r, "software"),
            xml_first(r, "softwareVersion"),
            xml_all(r),                      
        ))
    end
    crs
end

# Processor.xml → Components.Processor 
function parse_processors(path::String)
    ps = Processor[]
    for r in table_rows(path, "Processor.xml")
        push!(ps, Processor(
            nothing,                              
            xml_first(r, "modeId"),               
            xml_first(r, "processorSubType"),     
            xml_first(r, "processorType"),        
            xml_first(r, "processorId"),          
            xml_all(r),                           
        ))
    end
    ps
end

# State.xml → Components.State 
function parse_states(path::String)
    sts = State[]
    for r in table_rows(path, "State.xml")
        push!(sts, State(
            xml_first(r, "calDeviceName"),        
            nothing,                              
            nothing,                              
            nothing,                              
            xml_first(r, "ref") == "true",        
            xml_first(r, "sig") == "true",        
            nothing,                              
            xml_all(r),                           
        ))
    end
    sts
end

# Flag.xml → Components.Flag 
function parse_flags(path::String)
    flags = Flag[]
    for r in table_rows(path, "Flag.xml")
        t0 = to_int64(xml_first(r, "startTime"))
        t1 = to_int64(xml_first(r, "endTime"))
        push!(flags, Flag(
            nothing,                              
            nothing,                              
            t1 - t0,                              
            nothing,                              
            xml_first(r, "reason"),               
            nothing,                              
            t0,                                   
            nothing,                              
            xml_all(r),                          
        ))
    end
    flags
end

# Receiver.xml → Components.Receiver 
function parse_receivers(path::String)
    recs = Receiver[]
    for r in table_rows(path, "Receiver.xml")
        push!(recs, Receiver(
            xml_first(r, "name"),                       
            asdm_floats(xml_first(r, "freqLO")),        
            xml_first(r, "frequencyBand"),              
            to_int(xml_first(r, "numLO", "0")),            
            xml_first(r, "receiverId"),                 
            xml_first(r, "receiverSideband"),           
            asdm_array(xml_first(r, "sidebandLO")),     
            xml_first(r, "spectralWindowId"),           
            asdm_interval(xml_first(r, "timeInterval")),
            xml_all(r),                                 
        ))
    end
    recs
end

# Weather.xml → Components.Weather 
function parse_weather(path::String, stations::Vector{Station})
    sta_pos = Dict{String, Vector{Float64}}(s.stationID => s.position for s in stations)
    ws = Weather[]
    for r in table_rows(path, "Weather.xml")
        t0, dur = asdm_interval(xml_first(r, "timeInterval"))
        sid     = xml_first(r, "stationId")
        push!(ws, Weather(
            nothing,                              
            to_float(xml_first(r, "dewPoint", "0")),   
            nothing,                              
            dur,                                  
            sid,                                  
            get(sta_pos, sid, nothing),           
            to_float(xml_first(r, "pressure", "0")),   
            nothing,                              
            to_float(xml_first(r, "relHumidity", "0")),
            nothing,                              
            to_float(xml_first(r, "temperature", "0")),
            nothing,                              
            t0,                                   
            to_float(xml_first(r, "windDirection", "0")), 
            nothing,                              
            to_float(xml_first(r, "windSpeed", "0")),  
            nothing,                              
            xml_all(r),                           
        ))
    end
    ws
end

# CalDevice.xml → Components.CalDevice 
function parse_caldevices(path::String)
    cds = CalDevice[]
    for r in table_rows(path, "CalDevice.xml")
        t0, dur = asdm_interval(xml_first(r, "timeInterval"))
        calEff  = xml_opt(r, "calEff")
        tLoad   = xml_opt(r, "temperatureLoad")
        push!(cds, CalDevice(
            xml_first(r, "antennaId"),                      
            calEff === nothing ? nothing : asdm_floats_all(calEff),  
            asdm_array(xml_first(r, "calLoadNames")),        
            xml_first(r, "feedId"),                          
            dur,                                             
            asdm_floats_all(xml_first(r, "noiseCal")),       
            to_int(xml_first(r, "numCalload", "0")),            
            to_int(xml_first(r, "numReceptor", "0")),           
            xml_first(r, "spectralWindowId"),                
            tLoad === nothing ? nothing : asdm_floats_all(tLoad),    
            t0,                                              
            xml_all(r),                                      
        ))
    end
    cds
end

# Feed.xml → Components.Feed 
function parse_feeds(path::String)
    feeds = Feed[]
    for r in table_rows(path, "Feed.xml")
        t0, dur = asdm_interval(xml_first(r, "timeInterval"))
        push!(feeds, Feed(
            xml_first(r, "antennaId"),                       
            xml_opt(r, "beamId"),                            
            asdm_floats_all(xml_first(r, "beamOffset")),     
            xml_first(r, "feedId"),                          
            nothing,                                         
            dur,                                             
            to_int(xml_first(r, "numReceptor", "0")),           
            asdm_array(xml_first(r, "polarizationTypes")),   
            asdm_floats_all(xml_first(r, "polResponse")),    
            nothing,                                         
            asdm_floats_all(xml_first(r, "receptorAngle")),  
            xml_first(r, "spectralWindowId"),                
            t0,                                              
            xml_all(r),                                      
        ))
    end
    feeds
end

# Polarization product counts from ConfigDescription + DataDescription + Polarization.
function pol_counts(cds::Vector{ConfigDescription}, dds::Vector{DataDescrip},
                    pols::Vector{Polarization})
    isempty(cds) && return (4, 3)
    dd_ids = cds[1].dataDescriptionID
    isempty(dd_ids) && return (4, 3)
    first_dd = dd_ids[1]
    dd_idx = id_index(first_dd) + 1
    (1 <= dd_idx <= length(dds)) || return (4, 3)
    pol_id = dds[dd_idx].polarID
    pol_idx = id_index(pol_id) + 1
    (1 <= pol_idx <= length(pols)) || return (4, 3)
    n_cross = pols[pol_idx].numCorr
    n_auto = n_cross == 4 ? 3 : n_cross
    n_cross, n_auto
end

# Public API 
function open_sdm(path::String)
    isdir(path) || error("Not a directory: $path")

    observation, n_ant_hdr = parse_observation(path)
    mains       = parse_mains(path)
    scans       = parse_scans(path)
    subscans    = parse_subscans(path)
    stations    = parse_stations(path)
    antennas    = parse_antennas(path, stations)
    spws        = parse_spectral_windows(path)
    pols        = parse_polarizations(path)
    datadescs   = parse_datadescs(path)
    fields      = parse_fields(path)
    sources     = parse_sources(path)
    cds         = parse_configdescriptions(path)
    cms         = parse_correlatormodes(path)
    pms         = parse_pointingmodels(path)
    sbs         = parse_sbsummaries(path)
    scs         = parse_switchcycles(path)
    caldata     = parse_caldata(path)
    calreds     = parse_calreductions(path)
    caldevices  = parse_caldevices(path)
    feeds       = parse_feeds(path)
    flags       = parse_flags(path)
    procs       = parse_processors(path)
    states      = parse_states(path)
    receivers   = parse_receivers(path)
    weather     = parse_weather(path, stations)

    n_cross_pol, n_auto_pol = pol_counts(cds, datadescs, pols)
    n_ant  = isempty(antennas) ? n_ant_hdr : length(antennas)
    n_chan = isempty(spws) ? 0 : spws[1].numChan

    SDMDataset(
        path, observation, mains, scans, subscans, antennas, stations,
        spws, pols, datadescs, fields, sources, cds, cms, pms, sbs, scs,
        caldata, calreds, caldevices, feeds, flags, procs, states, receivers,
        weather,
        n_ant, length(spws), n_chan, n_cross_pol, n_auto_pol,
    )
end
