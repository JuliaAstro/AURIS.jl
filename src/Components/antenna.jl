struct Antenna{I<:Integer, S<:AbstractString}
    number::I
    name::S
    station
    position
    offset
    mount
    type
    area
    raw   # complete original row (lossless catch-all)
end
