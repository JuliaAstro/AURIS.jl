module SFitting

using LinearAlgebra, Statistics
using ..SProbability: prob_good

export fit_ls, fit_weighted, fit_masked, fit_em, evaluate

function fit_ls(Bx, By, Z; λ=1e-6)
    nkx, nky = size(Bx, 2), size(By, 2)
    n = nkx * nky
    AtA = kron(By' * By, Bx' * Bx)
    Atz = vec(Bx' * Z * By)
    c = (AtA + λ * I(n)) \ Atz
    reshape(c, nkx, nky)
end

function fit_weighted(Bx, By, Z, w; λ=1e-6)
    nkx, nky = size(Bx, 2), size(By, 2)
    n = nkx * nky
    W = reshape(w, size(Bx, 1), size(By, 1))

    AtWA = zeros(n, n)
    for j in axes(By, 1)
        Qj = Bx' * (W[:, j] .* Bx)
        byj = By[j, :]
        for d in 1:nky, b in 1:nky
            ib = (b-1)*nkx+1:b*nkx
            id = (d-1)*nkx+1:d*nkx
            @views AtWA[ib, id] .+= (byj[b] * byj[d]) .* Qj
        end
    end

    AtWz = vec(Bx' * (W .* Z) * By)
    c = (AtWA + λ * I(n)) \ AtWz
    reshape(c, nkx, nky)
end

function fit_masked(Bx, By, Z, mask; λ=1e-6)
    fit_weighted(Bx, By, Z, Float64.(mask); λ=λ)
end

function fit_em(
    Bx, By, Z;
    λ=1e-6,
    β::Float64=0.01,
    γ::Float64=10.0,
    maxiters::Int=30,
    tol=1e-6
)
    C = fit_ls(Bx, By, Z; λ=λ)

    for _ in 1:maxiters
        Ẑ = evaluate(Bx, By, C)
        resid = Z .- Ẑ
        med_r = median(resid)
        σ = max(1.4826 * median(abs.(resid .- med_r)), 1e-10)
        Z_fit_adj = max.(Ẑ .+ med_r, 1e-10)
        w = prob_good(Z, Z_fit_adj, σ; β=β, γ=γ)
        w[resid.≤0] .= 1.0
        C_new = fit_weighted(Bx, By, Z, w; λ=λ)

        Δ = norm(vec(C_new - C))
        C = C_new
        Δ < tol * (norm(vec(C)) + 1e-10) && break
    end

    return C
end

function evaluate(Bx, By, C)
    Bx * C * By'
end

end