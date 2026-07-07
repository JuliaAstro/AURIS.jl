using StatsBase

abstract type Threshold end
abstract type FlagMethod end

#=
function setweights(mset::MSet)

    Npol, Nchn, Ncor, Nexp, Nspw = size(mset.data)

    weight = deepcopy(mset.weight)
    for p in 1:Npol
        for (j::UInt32, s::UInt32) in collect(enumerate(mset.spw))
            datam = Array{Float64}(undef,Nchn)
            Threads.@threads for k in 1:Nchn
                I = findall(!iszero, mset.weight[p,k,:,:,j])
                datam[k] = mean(abs.(mset.data[p,k,:,:,j][I]))
            end
            datad = abs.(datam .- median(datam))
            datae = median(datad)
            weight[p,:,:,:,j] .*= map(x -> x<10. ? 1. : 0.,
                                      datae != 0. ? datad./datae :
                                      fill(100.0,size(datad)))
        end
    end
    return MSet(mset.antenna, mset.frequency, mset.spw, mset.field, mset.scan,
                mset.time, mset.uvw, weight, mset.data)
end
=#

const MAD_NORMAL = 1.4826

#=
  Assume data array axes are: polarization, channels, correlations
=#

"""
    winsorized_mean_and_std(mask, data)

Calculate the mean and standard deviation of the data.

Primarily used for data having Gaussian statistics.
Outliers are excluding by using only the middle 80% of the data. 
"""
function winsorized_mean_and_std(mask::AbstractArray{Bool}, data::AbstractArray{<:AbstractFloat})

    unmaskdata = vec(abs.(data[.!mask .&& isfinite.(data)]))
    unmasklen = length(unmaskdata)
    if unmasklen > 0
        #  Get values at 10% and 90% indices
        minindex, maxindex = floor(Int, 0.1*unmasklen)+1, ceil(Int, 0.9*unmasklen)
        minvalue, maxvalue = partialsort(unmaskdata, minindex), partialsort(unmaskdata, maxindex)
        #  Clip outliers to min and max values
        unmaskdata[unmaskdata .< minvalue] .= minvalue
        unmaskdata[unmaskdata .> maxvalue] .= maxvalue
        #  Calculate mean and std dev
        res = StatsBase.mean_and_std(unmaskdata)
    else
        res = (0., 0.)
    end
end

"""
    winsorized_mode(mask, data)

Calculate the mode of the data.

Primarily used for data having Rayleigh (or 2D chi-squared) statistics.
Outliers are excluding by using only the lower 90% of the data.
"""
function winsorized_mode(mask::AbstractArray{Bool}, data::AbstractArray{<:AbstractFloat})

    unmaskdata = vec(abs.(data[.!mask .&& isfinite.(data)]))
    unmasklen = length(unmaskdata)
    if unmasklen > 0
        #  Get values at 90% index
        maxindex = floor(Int, 0.8*unmasklen)
        maxvalue = partialsort(unmaskdata, maxindex)
        #  Clip outliers to max value
        unmaskdata[unmaskdata .> maxvalue] .= maxvalue
        #  Calculate mode
        res = StatsBase.mode(unmaskdata)
    else
        res = 0
    end
    res
end

"""
    gaussian_scaling(mask, data, threshold, timesen, freqsen; verbose=false)

Calculate the initial time and frequency scales for data having Gaussian noise.

Primarily used for real-valued data.
"""
function gaussian_scaling(mask::AbstractArray{Bool}, data::AbstractArray{<:AbstractFloat}, threshold::AbstractFloat,
                          timesen::AbstractFloat, freqsen::AbstractFloat; verbose=false)
    
    mean, stddev = winsorized_mean_and_std(mask, data)
    timescale, freqscale = (stddev == 0.0 ? 1.0 : stddev) .* (timesen, freqsen)
    if verbose
        println("Stddev=$stddev first time-direction threshold=$(threshold*timescale)")
    end
    timescale, freqscale
end

"""
    rayleigh_scaling(mask, data, threshold, timesen, freqsen; verbose=false)

Calculate the initial time and frequency scales for data having Rayleigh noise.

Primarily used for complex-valued data.
"""
function rayleigh_scaling(mask::AbstractArray{Bool}, data::AbstractArray{<:AbstractFloat}, threshold::AbstractFloat,
                          timesen::AbstractFloat, freqsen::AbstractFloat; verbose=false)

    mode = winsorized_mode(mask, data)
    timescale, freqscale = (mode == 0.0 ? 1.0 : mode) .* (timesen, freqsen)
    if verbose
        mean, stddev = winsorized_mean_and_std(mask, data)
        println("Mode=$mode first time-direction threshold=$(threshold*timescale)")
        println("Stddev=$stddev")
    end
    timescale, freqscale
end

function median_time_scaling(mask::AbstractArray{Bool}, data::AbstractArray{<:AbstractFloat}, threshold::AbstractFloat,
                             timesen::AbstractFloat, freqsen::AbstractFloat; verbose=false)
    median_ = median(abs.(data[.!mask]), dims=2)
    
    timescale, freqscale = median_ .* (timesen, freqsen)
    if verbose
        stddev = median_*MAD_NORMAL
        println("Median=$median_ first time-direction threshold=$(threshold*timescale)")
        println("Stddev=$stddev")
    end
    timescale, freqscale
end

"""
    constant_scaling(mask, data, threshold, timesen, freqsen, verbose=false)

Calculate the initial time and frequency scales for data.

Default method.
"""
function contant_scaling(mask::AbstractArray{Bool}, data::AbstractArray{<:AbstractFloat}, threshold::AbstractFloat,
                         timesen::AbstractFloat, freqsen::AbstractFloat; verbose=false)
    if verbose
        println("Stddev=0.0 first time-direction threshold=$timesen")
    end
    timesen, freqsen
end

#=
function filtersamples(mask, minsamples, connected)
    ylength, xlength size(mask)
    for j=1:ylength
end
=#

struct FreqThreshold <: Threshold
    length::Int
    value::Float64
end

struct TimeThreshold <: Threshold
    length::Int
    value::Float64
end

"""
    thresholds(method, axis, scaling)

Generate an array of lengths and thresholds for the given method, axis, and scaling.

Allowd methods are SumThreshold and VarThreshold and allowed axes are TimeThreshold and FreqThreshold.
"""
function thresholds(method::FlagMethod, axis::typeof(Threshold), scaling::AbstractFloat)
    # [axis(2^(j-1), scaling * method.initthresh / method.expfactor^log2(2^(j-1)) / 2^(j-1))
    # MeerKAT tricolour
    [axis(2^(j-1), scaling * method.initthresh / method.expfactor^log2(2^(j-1)))
     for j=1:(0<=method.iterations<=9 ? method.iterations : 1)]
end

#  Var-Threshold Algorithm

struct VarThreshold <: FlagMethod
    iterators::Int
    minsample::Int
    initthresh::Float64
    expfactor::Float64
end

"""
    VarThreshold(iterations, minsample, [initthresh=1.0, expfactor=1.2])

Configuration values for the var-threshold flagging algorithm.
"""
function (::VarThreshold)(mask::AbstractArray{Bool}, data::AbstractArray{<:Number}, threshold::TimeThreshold;
                          verbose::Bool=false)

    if size(mask) != size(data)
        DimensionMismatch("mask and data array do not have same dimensions.")
    end
    if verbose
        println("Performing VarThreshold with length $(threshold.length), threshold $(threshold.value)...")
    end
    nrows, ncols = size(data)
    if threshold.length <= ncols
        for k=1:ncols - threshold.length + 1
            for j=1:nrows
                slice = k:k + threshold.length - 1
                mask[j,slice] .= any(-threshold.value .<= data[j,slice] .<= threshold.value) ? false : true
            end
        end
    end
    mask
end

VarThreshold(iters, minsamp) = VarThreshold(iters, minsamp, 1.0, 1.2)

function (::VarThreshold)(mask::AbstractArray{Bool}, data::AbstractArray{<:Number}, threshold::FreqThreshold;
                          verbose::Bool=false)

    if size(mask) != size(data)
        DimensionMismatch("mask and data array do not have same dimensions.")
    end
    nrows, ncols = size(data)
    if threshold.length <= nrows
        for k=1:ncols
            for j=1:nrows - threshold.length + 1
                slice = j:j + threshold.length - 1
                mask[slice,k] .= any(-threshold.value .<= data[slice,k] .<= threshold.value) ? false : true
            end
        end
    end
    mask
end

#  Sum-Threshold Algorithm

"""
    SumThreshold(iterations, minsample, [initthresh=1.0, expfactor=1.5])

Configuration values for the sum-threshold flagging algorithm.
"""
struct SumThreshold <: FlagMethod
    iterations::Int
    minsample::Int
    initthresh::Float64
    expfactor::Float64
end

# MeerKAT tricolour uses rho = 1.3
SumThreshold(iters, minsamp) = SumThreshold(iters, minsamp, 1.0, 1.5)

"""
    (SumThreshold)(mask, data, threshold; verbose)

Functor for the given method.

Arguments are the mask to be modifed, the data to be flagged, and the threshold to be used.

Allowed thresholds are TimeThreshold and FreqThreshold. The mask and data size are identical.
"""
function (::SumThreshold)(mask::AbstractArray{Bool}, data::AbstractArray{<:Number}, threshold::TimeThreshold;
                          verbose::Bool=false)

    if size(mask) != size(data)
        DimensionMismatch("mask and data array do not have same dimensions.")
    end
    if verbose
        println("Performing SumThreshold with length $(threshold.length), threshold $(threshold.value)...")
    end
    nrows, ncols = size(data)
    if threshold.length <= ncols
        for k=1:ncols - threshold.length + 1
            for j=1:nrows
                slice = k:k+threshold.length-1
                sflag = sum(abs.(data[j,slice][.!mask[j,slice]]))
                nflag = count(.!mask[j,slice])
                if nflag > 0 && abs(sflag/nflag) > threshold.value
                    mask[j,slice] .= true
                end
            end
        end
    end
    mask
end

function (::SumThreshold)(mask::AbstractArray{Bool}, data::AbstractArray{<:Number},
                          threshold::FreqThreshold; verbose::Bool=false)

    if size(mask) != size(data)
        DimensionMismatch("mask and data array do not have same dimensions.")
    end
    nrows, ncols = size(data)
    if threshold.length <= nrows
        for k=1:ncols
            for j=1:nrows - threshold.length + 1
                slice = j:j+threshold.length-1
                sflag = sum(abs.(data[slice,k][.!mask[slice,k]]))
                nflag = count(.!mask[slice,k])
                if nflag > 0 && abs(sflag/nflag) > threshold.value
                    mask[slice,k] .= true
                end
            end
        end
    end
    mask
end

function flagdata!(mask::AbstractArray{Bool}, data::AbstractArray{<:Number},
                   method::FlagMethod=SumThreshold, statistic=rayleigh_scaling;
                   freqsen::AbstractFloat=1., timesen::AbstractFloat=1.,
                   additive::Bool=true, verbose::Bool=false)

    rdata = abs.(data)
    timescaling, freqscaling = statistic(mask, rdata, timesen, freqsen, method.initthresh;
                                         verbose=verbose)

    if !additive mask .= false end

    time, freq = thresholds(method, TimeThreshold, timescaling), thresholds(method, FreqThreshold, freqscaling)

    for j=1:max(length(time), length(freq))
        if j < length(time) method(mask, rdata, time[j], verbose=verbose) end
        if j < length(freq) method(mask, rdata, freq[j], verbose=verbose) end
    end

    if method.minsample > 1 filtersamples(mask, method.minsample) end
    mask
end
