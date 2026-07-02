# Binary Data Format (BDF) visibility reader

# Data structures
struct IntegrationMeta
    time_ns        :: Int64    
    interval_ns    :: Int64    
    cross_offset   :: Int64    
    auto_offset    :: Int64    
end

struct BDFFile
    path         :: String
    iostream     :: IOStream   # kept open to hold the mmap
    mmap         :: Vector{UInt8}
    n_antenna    :: Int
    n_baseline   :: Int
    n_spw        :: Int
    n_chan       :: Int
    n_cross_pol  :: Int        
    n_auto_pol   :: Int        
    cross_size   :: Int        
    auto_size    :: Int        
    integrations :: Vector{IntegrationMeta}
end

Base.length(b::BDFFile) = length(b.integrations)

function Base.close(b::BDFFile)
    finalize(b.mmap)
    close(b.iostream)
end
const close_bdf = close

# Size calculations
cross_data_size(n_bl, n_spw, n_chan, n_pol) = n_bl * n_spw * n_chan * n_pol * 8

auto_data_size(n_ant, n_spw, n_chan, n_auto_pol) =
    n_ant * n_spw * n_chan * ((n_auto_pol - 1) * 4 + 8)

# BDF header parsing 
function parse_subheader_time(line::String)
    tm = match(r"<time>(\d+)</time>", line)
    iv = match(r"<interval>(\d+)</interval>", line)
    t  = tm === nothing ? Int64(0) : parse(Int64, tm.captures[1])
    i  = iv === nothing ? Int64(0) : parse(Int64, iv.captures[1])
    t, i
end

# MIME structure scan 
function scan_bdf_structure(path::String, cross_size::Int, auto_size::Int)
    integrations = IntegrationMeta[]

    open(path, "r") do io
        time_ns     = Int64(0)
        interval_ns = Int64(0)
        cross_off   = Int64(0)
        auto_off    = Int64(0)
        got_cross   = false

        state = :seek_global_header

        while !eof(io)
            line = readline(io; keep=false)

            if state == :seek_global_header
                if startswith(line, "--MIME_boundary-1")
                    state = :in_outer_header
                end

            elseif state == :in_outer_header
                if isempty(line)
                    state = :in_outer_body
                end

            elseif state == :in_outer_body
                if startswith(line, "--MIME_boundary-1")
                    state = :in_outer_header
                elseif contains(line, "sdmDataHeader")
                elseif contains(line, "multipart/related")
                    state = :in_integration
                elseif startswith(line, "--MIME_boundary-2")
                    got_cross = false
                    state = :in_inner_header
                elseif contains(line, "sdmDataSubsetHeader")
                    time_ns, interval_ns = parse_subheader_time(line)
                end

            elseif state == :in_integration
                if startswith(line, "--MIME_boundary-2")
                    state = :in_inner_header
                end

            elseif state == :in_inner_header
                if contains(line, "crossData.bin")
                    state = :after_cross_header
                elseif contains(line, "autoData.bin")
                    state = :after_auto_header
                elseif contains(line, "desc.xml")
                    state = :in_desc_header
                elseif isempty(line)
                    state = :in_outer_body
                end

            elseif state == :in_desc_header
                if isempty(line)
                    state = :in_desc_body
                end

            elseif state == :in_desc_body
                if contains(line, "sdmDataSubsetHeader")
                    time_ns, interval_ns = parse_subheader_time(line)
                elseif startswith(line, "--MIME_boundary-2")
                    state = :in_inner_header
                end

            elseif state == :after_cross_header
                if isempty(line)
                    cross_off = position(io)
                    skip(io, cross_size)
                    eof(io) || readline(io; keep=false)
                    got_cross = true
                    state = :in_integration
                end

            elseif state == :after_auto_header
                if isempty(line)
                    auto_off = position(io)
                    skip(io, auto_size)
                    eof(io) || readline(io; keep=false)
                    if got_cross
                        push!(integrations, IntegrationMeta(
                            time_ns, interval_ns, cross_off, auto_off))
                        got_cross = false
                    end
                    state = :in_integration
                end
            end
        end
    end

    integrations
end

# Public API 
function open_bdf(dataset::SDMDataset, main::Main)
    path = bdf_path(dataset, main)
    isfile(path) || error("BDF file not found: $path")

    n_ant  = dataset.n_antenna
    n_bl   = n_ant * (n_ant - 1) ÷ 2
    n_spw  = dataset.n_spw
    n_chan = dataset.n_chan
    n_xpol = dataset.n_cross_pol
    n_apol = dataset.n_auto_pol

    xsize = cross_data_size(n_bl, n_spw, n_chan, n_xpol)
    asize = auto_data_size(n_ant, n_spw, n_chan, n_apol)

    integrations = scan_bdf_structure(path, xsize, asize)

    if length(integrations) != main.numIntegration
        @warn "Expected $(main.numIntegration) integrations, found $(length(integrations)) in $(basename(path))"
    end

    io   = open(path, "r")
    mmap = Mmap.mmap(io, Vector{UInt8})

    BDFFile(path, io, mmap, n_ant, n_bl, n_spw, n_chan, n_xpol, n_apol,
            xsize, asize, integrations)
end

# Data access
const MJD_TO_UNIX_NS = Int64(40587) * Int64(86400) * Int64(1_000_000_000)

function integration_times(bdf::BDFFile)
    [(m.time_ns - MJD_TO_UNIX_NS) * 1e-9 for m in bdf.integrations]
end

function cross_view(bdf::BDFFile, t::Int)
    m   = bdf.integrations[t]
    raw = view(bdf.mmap, (m.cross_offset + 1):(m.cross_offset + bdf.cross_size))
    cf  = reinterpret(ComplexF32, raw)
    reshape(cf, bdf.n_cross_pol, bdf.n_chan, bdf.n_spw, bdf.n_baseline)
end

const POL_IDX = Dict{Symbol,Int}(:RR => 1, :RL => 2, :LR => 3, :LL => -1)

function load_amplitudes(bdf::BDFFile; stokes::Symbol=:I)
    n_t  = length(bdf.integrations)
    n_f  = bdf.n_spw * bdf.n_chan
    n_bl = bdf.n_baseline
    n_p  = bdf.n_cross_pol

    out = Array{Float16}(undef, n_t, n_f, n_bl)

    ll_idx  = n_p   
    pol_sel = stokes == :I  ? 0 :
              stokes == :LL ? ll_idx :
              get(POL_IDX, stokes, 1)

    for t in 1:n_t
        block = cross_view(bdf, t)  
        for bl in 1:n_bl
            fi = 1
            for spw in 1:bdf.n_spw, c in 1:bdf.n_chan
                amp = if pol_sel == 0
                    0.5f0 * (abs(block[1, c, spw, bl]) + abs(block[ll_idx, c, spw, bl]))
                else
                    abs(block[pol_sel, c, spw, bl])
                end
                out[t, fi, bl] = Float16(amp)
                fi += 1
            end
        end
    end

    out
end
