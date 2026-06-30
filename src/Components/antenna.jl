struct Antenna{I,S} where {I<:Integer, S<:AbstractString}
    number::I
    name::S
    station
    position
    offset
    mount
    type
    area
end
