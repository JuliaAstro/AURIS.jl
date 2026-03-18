module SBasis

export basis, bases, knot_vector

function basis(t, k, p, x)
    n = length(t) - p - 1
    if k == n && x == t[end]
        return 1.0
    end
    if x < t[k] || x ≥ t[k+p+1]
        return 0.0
    end
    if p == 0
        return 1.0
    end

    d1 = t[k+p] - t[k]
    d2 = t[k+p+1] - t[k+1]

    a = d1 == 0 ? 0.0 : (x - t[k]) / d1 * basis(t, k, p - 1, x)
    b = d2 == 0 ? 0.0 : (t[k+p+1] - x) / d2 * basis(t, k + 1, p - 1, x)

    a + b
end

function bases(kx, ky, p, ξ, η)
    nbx = length(kx) - p - 1
    nby = length(ky) - p - 1

    Bx = [basis(kx, a, p, ξi) for ξi in ξ, a in 1:nbx]
    By = [basis(ky, b, p, ηj) for ηj in η, b in 1:nby]

    Bx, By
end

function knot_vector(a, b, n, p)
    n_inner = n - p - 1

    if n_inner < 0
        error("Invalid spline: n=$n < p+1=$(p+1)")
    elseif n_inner == 0
        return vcat(fill(a, p + 1), fill(b, p + 1))
    end

    inner = range(a, b, length=n_inner + 2)[2:end-1]
    vcat(fill(a, p + 1), inner, fill(b, p + 1))
end

end