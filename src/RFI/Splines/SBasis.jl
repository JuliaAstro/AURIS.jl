module SBasis

export bspline_basis, basis_matrices, clamped_knots_inclusive, decode_knots

### Cox-de Boor recursion
function bspline_basis(t, k, p, x)
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

    a = d1 == 0 ? 0.0 : (x - t[k]) / d1 * bspline_basis(t, k, p - 1, x)
    b = d2 == 0 ? 0.0 : (t[k+p+1] - x) / d2 * bspline_basis(t, k + 1, p - 1, x)

    a + b
end

function basis_matrices(kx, ky, p, ξ, η)
    nbx = length(kx) - p - 1
    nby = length(ky) - p - 1

    Bx = [bspline_basis(kx, a, p, ξi) for ξi in ξ, a in 1:nbx]
    By = [bspline_basis(ky, b, p, ηj) for ηj in η, b in 1:nby]

    Bx, By
end

### Knots

function clamped_knots_inclusive(a, b, n, p)
    n_inner = n - p - 1

    if n_inner < 0
        error("Invalid spline: n=$n < p+1=$(p+1)")
    elseif n_inner == 0
        return vcat(fill(a, p + 1), fill(b, p + 1))
    end

    inner = range(a, b, length=n_inner + 2)[2:end-1]
    vcat(fill(a, p + 1), inner, fill(b, p + 1))
end


function decode_knots(v, n, p; δ=1e-3)
    n_inner = n - p - 1
    inner = sort(clamp.(v[1:n_inner], δ, 1 - δ))
    for i in 2:length(inner)
        inner[i] = max(inner[i], inner[i-1] + δ)
    end
    vcat(fill(0.0, p + 1), inner, fill(1.0, p + 1))
end

end