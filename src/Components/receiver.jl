struct Receiver
    name
    freqLO
    freqBand
    numLO
    receiverID
    receiverSideband
    sidebandLO
    spectralWindowID
    timeInterval
    raw   # complete original row (lossless catch-all)
end

struct Circular{F<:AbstractFloat}
    L::F
    R::F
end

struct Linear{F<:AbstractFloat}
    X::F
    Y::F
end
