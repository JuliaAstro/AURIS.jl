module SOptimisation

using ..SBasis: clamped_knots_inclusive, basis_matrices, decode_knots
using ..SFitting: fit_bspline_ls, eval_surface, fit_bspline_robust, fit_bspline_ls_weighted, fit_bspline_mask_outliers, fit_bspline_ls_masked
using Random, LinearAlgebra, Statistics

export pso, model_loss, knot_loss

function pso(f; ndim, nparticles=30, maxiters=40, bounds, ω=0.5, c1=1.5, c2=1.5)
    low, high = bounds
    X = low .+ (high - low) .* rand(nparticles, ndim)
    V = zeros(nparticles, ndim)
    pbest = copy(X)
    pbest_val = [f(X[i, :]) for i in 1:nparticles]
    gidx = argmin(pbest_val)
    gbest = copy(pbest[gidx, :])
    gbest_val = pbest_val[gidx]

    for _ in 1:maxiters
        for i in 1:nparticles
            V[i, :] .= ω .* V[i, :] .+ c1 .* rand(ndim) .* (pbest[i, :] .- X[i, :]) .+ c2 .* rand(ndim) .* (gbest .- X[i, :])
            X[i, :] .+= V[i, :]
            clamp!(X[i, :], low, high)
            val = f(X[i, :])
            if val < pbest_val[i]
                pbest[i, :] .= X[i, :]
                pbest_val[i] = val
                if val < gbest_val
                    gbest .= X[i, :]
                    gbest_val = val
                end
            end
        end
    end
    gbest, gbest_val
end

### Knot number loss
function model_loss(v, Z, ξ, η, p; β=0.01, γ=10.0, robust=true)
    nbx = clamp(round(Int, v[1]), p + 1, 25)
    nby = clamp(round(Int, v[2]), p + 1, 25)
    kx = clamped_knots_inclusive(0.0, 1.0, nbx, p)
    ky = clamped_knots_inclusive(0.0, 1.0, nby, p)
    Bx, By = basis_matrices(kx, ky, p, ξ, η)
    C = robust ? fit_bspline_robust(Bx, By, Z; β=β, γ=γ) : fit_bspline_ls(Bx, By, Z)
    Z_fit = eval_surface(Bx, By, C)
    rmse = sqrt(mean((Z_fit .- Z) .^ 2))
    complexity = 1e-3 * nbx * nby
    rmse + complexity
end

### Knot positions loss
function knot_loss(v, Z, ξ, η, p, nbx_opt, nby_opt; β=0.01, γ=10.0, robust=true)
    mx = nbx_opt - p - 1
    kx = decode_knots(v[1:mx], nbx_opt, p)
    ky = decode_knots(v[mx+1:end], nby_opt, p)
    Bx, By = basis_matrices(kx, ky, p, ξ, η)
    C = robust ? fit_bspline_robust(Bx, By, Z; β=β, γ=γ) : fit_bspline_ls(Bx, By, Z)
    Z_fit = eval_surface(Bx, By, C)
    sqrt(mean((Z_fit .- Z) .^ 2))
end

end
