function write_fits(path::AbstractString, image::AbstractMatrix{<:Real};
                    ra::Real, dec::Real, cell::Real, freq::Real,
                    bandwidth::Real=0.0, beam=nothing,
                    bunit::AbstractString="JY/BEAM", object::AbstractString="",
                    telescope::AbstractString="EVLA",
                    date_obs::AbstractString="",
                    origin::AbstractString="AURIS.jl")
    nx, ny = size(image)
    celldeg = cell * 180 / π
    scard(key, value) = Card(key, value; slash = 72)
    fcard(key, value, comment="") =
        Card(key, round(Float64(value); sigdigits = 12), comment)
    data = Array{Float32,4}(undef, nx, ny, 1, 1)
    data[:, :, 1, 1] .= Float32.(reverse(image, dims=1))

    cards = [Card("SIMPLE", true, "conforms to FITS standard"),
             Card("BITPIX", -32, "IEEE 32-bit float"),
             Card("NAXIS", 4),
             Card("NAXIS1", nx), Card("NAXIS2", ny),
             Card("NAXIS3", 1), Card("NAXIS4", 1),
             scard("BUNIT", bunit),
             scard("BTYPE", "Intensity"),
             scard("OBJECT", object),
             scard("TELESCOP", telescope),
             scard("ORIGIN", origin),
             Card("EQUINOX", 2000.0),
             scard("RADESYS", "FK5"),
             scard("SPECSYS", "TOPOCENT"),
             scard("CTYPE1", "RA---SIN"),
             fcard("CRVAL1", mod(ra * 180 / π, 360.0), "deg"),
             fcard("CDELT1", -celldeg, "deg"),
             fcard("CRPIX1", nx - nx ÷ 2),
             scard("CUNIT1", "deg"),
             scard("CTYPE2", "DEC--SIN"),
             fcard("CRVAL2", dec * 180 / π, "deg"),
             fcard("CDELT2", celldeg, "deg"),
             fcard("CRPIX2", ny ÷ 2 + 1),
             scard("CUNIT2", "deg"),
             scard("CTYPE3", "FREQ"),
             fcard("CRVAL3", freq, "Hz"),
             fcard("CDELT3", bandwidth > 0 ? bandwidth : 1.0, "Hz"),
             fcard("CRPIX3", 1.0),
             scard("CUNIT3", "Hz"),
             scard("CTYPE4", "STOKES"),
             fcard("CRVAL4", 1.0, "Stokes I"),
             fcard("CDELT4", 1.0),
             fcard("CRPIX4", 1.0)]
    isempty(date_obs) || push!(cards, scard("DATE-OBS", date_obs))
    if beam !== nothing
        push!(cards, fcard("BMAJ", beam[1] * 180 / π, "deg, restoring beam"))
        push!(cards, fcard("BMIN", beam[2] * 180 / π, "deg"))
        push!(cards, fcard("BPA", beam[3] * 180 / π, "deg"))
    end

    write(path, HDU[HDU{FITSFiles.Primary}(cards, data)])
    path
end
