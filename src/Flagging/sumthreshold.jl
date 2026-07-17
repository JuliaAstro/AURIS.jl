# To be replaced by Paul's implementation
function sumthreshold!(mask::AbstractMatrix{Bool}, plane::AbstractMatrix{<:Real},
                       t1::Real; seqs=(1, 2, 4, 8, 16, 32, 64), rho::Real=1.5)
    size(mask) == size(plane) || throw(DimensionMismatch("mask/plane size mismatch"))
    for M in seqs
        tM = t1 / rho^log2(M)
        _sumthreshold_axis!(mask, plane, M, tM, Val(1))
        _sumthreshold_axis!(mask, plane, M, tM, Val(2))
    end
    mask
end

function _sumthreshold_axis!(mask, plane, M::Int, tM::Real, ::Val{ax}) where {ax}
    nlines = size(plane, ax == 1 ? 2 : 1)
    nsamp  = size(plane, ax)
    nsamp >= M || return
    newmask = falses(nsamp)
    line  = Vector{Float64}(undef, nsamp)
    lmask = Vector{Bool}(undef, nsamp)
    for k in 1:nlines
        for i in 1:nsamp
            line[i]  = ax == 1 ? plane[i, k] : plane[k, i]
            lmask[i] = ax == 1 ? mask[i, k]  : mask[k, i]
        end
        fill!(newmask, false)
        s = 0.0
        for i in 1:M
            s += lmask[i] ? tM : line[i]
        end
        for i0 in 1:(nsamp - M + 1)
            if s > tM * M
                @views newmask[i0:i0+M-1] .= true
            end
            if i0 + M <= nsamp
                s -= lmask[i0] ? tM : line[i0]
                s += lmask[i0+M] ? tM : line[i0+M]
            end
        end
        for i in 1:nsamp
            newmask[i] || continue
            if ax == 1
                mask[i, k] = true
            else
                mask[k, i] = true
            end
        end
    end
end

function mad_sigma(plane::AbstractMatrix{<:Real}, mask::AbstractMatrix{Bool})
    vals = [Float64(plane[i]) for i in eachindex(plane) if !mask[i]]
    isempty(vals) && return 0.0, 0.0
    med = median(vals)
    1.4826 * median(abs.(vals .- med)), med
end

function flag_rfi!(v::VisibilityDataset; threshold::Real=6.0,
                   seqs=(1, 2, 4, 8, 16, 32, 64), rho::Real=1.5, niter::Int=2)
    npol  = Data.n_pol(v)
    rr, ll = 1, npol
    nchan, nspw = Data.n_chan(v), Data.n_spw(v)
    nbl, nt = Data.n_baseline(v), Data.n_time(v)

    nbefore = count(v.flags)
    work = [(bl, s) for s in 1:nspw for bl in 1:nbl]
    chunks = collect(Iterators.partition(work,
                     cld(length(work), Threads.nthreads())))
    Threads.@threads for chunk in chunks
        plane = Matrix{Float64}(undef, nchan, nt)
        mask  = Matrix{Bool}(undef, nchan, nt)
        for (bl, s) in chunk
            for t in 1:nt, c in 1:nchan
                plane[c, t] = 0.5 * (abs(ComplexF64(v.vis[rr, c, s, bl, t])) +
                                     abs(ComplexF64(v.vis[ll, c, s, bl, t])))
                mask[c, t]  = v.flags[rr, c, s, bl, t] || v.flags[ll, c, s, bl, t]
            end
            all(mask) && continue
            for c in 1:nchan
                vals = [plane[c, t] for t in 1:nt if !mask[c, t]]
                m = isempty(vals) ? 1.0 : median(vals)
                m = m > 0 ? m : 1.0
                @views plane[c, :] ./= m
            end
            for _ in 1:niter
                σ, med = mad_sigma(plane, mask)
                σ > 0 || break
                centred = plane .- med
                sumthreshold!(mask, centred, threshold * σ; seqs, rho)
            end
            for t in 1:nt, c in 1:nchan
                if mask[c, t]
                    @views v.flags[:, c, s, bl, t] .= true
                end
            end
        end
    end
    (count(v.flags) - nbefore) / max(length(v.flags) - nbefore, 1)
end

function flag_edges!(v::VisibilityDataset, nedge::Int=3)
    nchan = Data.n_chan(v)
    nedge >= 1 || return v
    v.flags[:, 1:min(nedge, nchan), :, :, :] .= true
    v.flags[:, max(1, nchan - nedge + 1):nchan, :, :, :] .= true
    v
end