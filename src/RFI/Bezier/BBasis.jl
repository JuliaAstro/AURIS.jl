module BBasis

using Statistics

export bezier_basis, get_tile_indices, get_sliding_window_indices, reflect_pad_2d, safe_median, hann1d, tile_window

function bezier_basis(t, n=3)
    [binomial(n, i) * (1 - t)^(n - i) * t^i for i in 0:n]
end

function get_tile_indices(npoints::Int, ntiles::Int, overlap_frac::Float64)
    tile_size = max(1, ceil(Int, npoints / ntiles))
    step = max(1, ceil(Int, tile_size * (1 - overlap_frac)))
    starts = collect(1:step:max(1, npoints - 1))
    ends = map(s -> min(s + tile_size - 1, npoints), starts)
    collect(zip(starts, ends))
end

function get_sliding_window_indices(npoints::Int, window_size::Int)
    starts = collect(1:(npoints-window_size+1))
    ends = starts .+ (window_size - 1)
    collect(zip(starts, ends))
end

function reflect_pad_2d(y::AbstractMatrix{T}, pad_t::Int=6, pad_f::Int=6) where T
    nt, nf = size(y)
    pad_t = min(pad_t, nt - 1)
    pad_f = min(pad_f, nf - 1)

    top = y[pad_t:-1:1, :]
    bottom = y[end:-1:end-pad_t+1, :]
    ytmp = vcat(top, y, bottom)

    left = ytmp[:, pad_f:-1:1]
    right = ytmp[:, end:-1:end-pad_f+1]

    hcat(left, ytmp, right)
end

function safe_median(v)
    isempty(v) ? 0.0 : median(v)
end

function hann1d(n::Int)
    n == 1 ? ones(1) : 0.5 .* (1 .- cos.(2π .* (0:(n-1)) ./ (n - 1)))
end

function tile_window(nt::Int, nf::Int; γ_freq::Float64=0.5)
    w_t = hann1d(nt)
    w_f = hann1d(nf)
    w_f .= w_f .^ γ_freq
    return w_t * w_f'
end

end