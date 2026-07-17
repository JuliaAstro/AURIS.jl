# SDM SysPower binary table reader

read_be(io::IO, ::Type{T}) where {T} = ntoh(read(io, T))

function read_string(io::IO)
    n = read_be(io, Int32)
    String(read(io, n))
end

function read_optional(io::IO)
    read(io, UInt8) == 0x00 && return nothing
    n = read_be(io, Int32)
    Float32[read_be(io, Float32) for _ in 1:n]
end

function find_bytes(data::Vector{UInt8}, pattern::Vector{UInt8}, from::Integer)
    n, m = length(data), length(pattern)
    @inbounds for i in from:(n - m + 1)
        matched = true
        for j in 1:m
            if data[i + j - 1] != pattern[j]
                matched = false
                break
            end
        end
        matched && return i
    end
    0
end

# Public API

load_syspower(dataset::SDMDataset; kwargs...) = load_syspower(dataset.path; kwargs...)

function load_syspower(path::AbstractString; limit::Integer=typemax(Int))
    file = joinpath(path, "SysPower.bin")
    isfile(file) || error("SysPower.bin not found in $path")
    data = read(file)

    cid = find_bytes(data, Vector{UInt8}("<content.bin>"), 1)
    cid == 0 && error("SysPower.bin: <content.bin> part not found")
    blank = find_bytes(data, UInt8[0x0a, 0x0a], cid)          
    blank == 0 && error("SysPower.bin: malformed content headers")
    start = blank + 2
    boundary = find_bytes(data, Vector{UInt8}("--MIME_boundary"), start)
    stop = (boundary == 0 ? length(data) + 1 : boundary) - 1  
    while stop >= start && (data[stop] == 0x0a || data[stop] == 0x0d)
        stop -= 1                                             
    end

    io = IOBuffer(view(data, start:stop))
    for _ in 1:10                                             
        read_string(io)
    end
    read(io, UInt32)                                          

    rows = SysPower[]
    while !eof(io) && length(rows) < limit
        antenna       = read_string(io)
        spw           = read_string(io)
        feed          = Int(read_be(io, Int32))
        timestamp     = read_be(io, Int64)
        interval      = read_be(io, Int64)
        nreceptor     = Int(read_be(io, Int32))
        switched_diff = read_optional(io)
        switched_sum  = read_optional(io)
        requant_gain  = read_optional(io)

        push!(rows, SysPower(antenna, feed, interval, requant_gain, spw,
                             switched_diff, switched_sum, timestamp,
                             (numReceptor = nreceptor,)))
    end
    rows
end
