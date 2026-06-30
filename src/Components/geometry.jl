struct Geometry{} where {I<:Integer, S<:AbstractString}
    n::I
    name::S
    station::Tuple
    deriv::Tuple
    mount
    antenna::Tuple
end
