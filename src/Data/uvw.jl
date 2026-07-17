const JD_MJD0 = 2400000.5   # Julian date of MJD = 0

"""
    tai_ns_to_jd_utc(t_ns) → Float64

Convert an ASDM timestamp (ns since MJD=0, TAI) to a UTC Julian date, applying
the leap-second table (ΔAT = TAI − UTC).
"""
function tai_ns_to_jd_utc(t_ns::Integer)
    jd_tai = JD_MJD0 + t_ns / 86400e9
    jd_tai - get_Δat(jd_tai) / 86400
end

"""
    itrf_to_j2000(jd_utc; dut1=0.0, eop=nothing) → SMatrix{3,3}

Rotation matrix taking ITRF/ECEF vectors to the J2000 celestial frame at the
given UTC Julian date.
"""
function itrf_to_j2000(jd_utc::Real; dut1::Real=0.0, eop=nothing)
    D = eop === nothing ?
        r_ecef_to_eci(PEF(), J2000(), jd_utc + dut1 / 86400) :
        r_ecef_to_eci(ITRF(), GCRF(), jd_utc, eop)
    SMatrix{3,3,Float64}(D)
end

"""
    uvw_basis(ra, dec) → (û, v̂, ŵ)

Orthonormal (u,v,w) basis vectors in the J2000 frame for a phase centre at
J2000 `ra`/`dec` (radians): ŵ towards the source, û east, v̂ north.
"""
function uvw_basis(ra::Real, dec::Real)
    w = SVector(cos(dec) * cos(ra), cos(dec) * sin(ra), sin(dec))
    u = SVector(-sin(ra), cos(ra), 0.0)          
    v = w × u
    u, v, w
end

"""
    uvw!(out, baselines_itrf, jd_utc, ra, dec; dut1=0.0, eop=nothing)

Fill `out` (3 × n_baseline, metres) with (u,v,w) for each ITRF baseline vector
(columns of `baselines_itrf`, metres) at UTC Julian date `jd_utc`, phase centre
J2000 `ra`/`dec`.

Baseline sign convention: pass B = r(ant2) − r(ant1) for the pair (ant1,ant2).
"""
function uvw!(out::AbstractMatrix{<:AbstractFloat},
              baselines_itrf::AbstractMatrix{<:Real},
              jd_utc::Real, ra::Real, dec::Real; dut1::Real=0.0, eop=nothing)
    size(out) == size(baselines_itrf) ||
        throw(DimensionMismatch("out and baselines_itrf must both be 3×n"))
    R = itrf_to_j2000(jd_utc; dut1, eop)
    û, v̂, ŵ = uvw_basis(ra, dec)
    for k in axes(baselines_itrf, 2)
        B = R * SVector{3,Float64}(baselines_itrf[1, k], baselines_itrf[2, k],
                                   baselines_itrf[3, k])
        out[1, k] = û ⋅ B
        out[2, k] = v̂ ⋅ B
        out[3, k] = ŵ ⋅ B
    end
    out
end
