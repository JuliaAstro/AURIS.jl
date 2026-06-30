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
end

struct Circular{F} where F<:AbstractFloat
    L::F
    R::F
end

struct Linear{F} where F<:AbstractFloat
    X::F
    Y::F
end
