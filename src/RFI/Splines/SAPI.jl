module SAPI

using LinearAlgebra, Statistics, Random

using ..SBasis: basis_matrices, clamped_knots_inclusive, decode_knots
using ..SFitting: fit_bspline_ls, eval_surface, fit_bspline_robust, fit_bspline_ls_weighted, fit_bspline_mask_outliers, fit_bspline_ls_masked
using ..SOptimisation: pso, model_loss, knot_loss

export fit_bspline_surface_pso

function fit_bspline_surface_pso(
    Z::AbstractMatrix;
    p::Int=3,
    knot_particles=40,
    knot_iters=50,
    knotpos_particles=60,
    knotpos_iters=80,
    knotpos_restarts::Int=3,
    λ=1e-6,
    β=0.01,
    γ=10.0,
    prob_threshold=0.5,
    pso_downsample::Int=4,
)

    nx, ny = size(Z)

    ε = eps(Float64)
    ξ = [clamp(xi, ε, 1 - ε) for xi in range(0.0, 1.0, length=nx)]
    η = [clamp(yj, ε, 1 - ε) for yj in range(0.0, 1.0, length=ny)]

    ds = max(1, pso_downsample)
    Z_ds = Z[1:ds:end, 1:ds:end]
    nx_ds, ny_ds = size(Z_ds)
    ξ_ds = [clamp(xi, ε, 1 - ε) for xi in range(0.0, 1.0, length=nx_ds)]
    η_ds = [clamp(yj, ε, 1 - ε) for yj in range(0.0, 1.0, length=ny_ds)]

    knot_max = max(p + 1, min(nx, ny) ÷ 4)

    best_v, _ = pso(
        v -> model_loss(v, Z_ds, ξ_ds, η_ds, p; β=β, γ=γ, robust=false),
        ndim=2,
        bounds=(Float64(p + 1), Float64(knot_max)),
        nparticles=knot_particles,
        maxiters=knot_iters
    )

    nbx_opt = clamp(round(Int, best_v[1]), p + 1, knot_max)
    nby_opt = clamp(round(Int, best_v[2]), p + 1, knot_max)

    ndim_knots = (nbx_opt - p - 1) + (nby_opt - p - 1)

    best_knots = nothing
    best_loss = Inf
    for _ in 1:knotpos_restarts
        knots, _ = pso(
            v -> knot_loss(v, Z_ds, ξ_ds, η_ds, p, nbx_opt, nby_opt; β=β, γ=γ, robust=false),
            ndim=ndim_knots,
            bounds=(0.0, 1.0),
            nparticles=knotpos_particles,
            maxiters=knotpos_iters
        )
        loss = knot_loss(knots, Z, ξ, η, p, nbx_opt, nby_opt; β=β, γ=γ, robust=false)
        if loss < best_loss
            best_loss = loss
            best_knots = knots
        end
    end

    mx = nbx_opt - p - 1
    kx_kopt = decode_knots(best_knots[1:mx], nbx_opt, p)
    ky_kopt = decode_knots(best_knots[mx+1:end], nby_opt, p)

    Bx_kopt, By_kopt = basis_matrices(kx_kopt, ky_kopt, p, ξ, η)

    C_robust = fit_bspline_robust(Bx_kopt, By_kopt, Z; λ=λ, β=β, γ=γ)

    C_kopt, mask = fit_bspline_mask_outliers(
        Bx_kopt, By_kopt, Z;
        λ=λ, β=β, γ=γ, prob_threshold=prob_threshold,
        C_init=C_robust
    )

    Z_kopt = eval_surface(Bx_kopt, By_kopt, C_kopt)

    return Z_kopt, mask
end

end
