module SAPI

using LinearAlgebra, Statistics, Random

using ..SBasis: basis_matrices, clamped_knots_inclusive, decode_knots
using ..SFitting: fit_bspline_ls, eval_surface, fit_bspline_robust, fit_bspline_ls_weighted, fit_bspline_mask_outliers, fit_bspline_ls_masked
using ..SOptimisation: pso, model_loss, knot_loss

export fit_bspline_surface_pso

function fit_bspline_surface_pso(
    Z::AbstractMatrix;
    p::Int=3,
    knot_bounds=(4.0, 25.0),
    knot_particles=40,
    knot_iters=50,
    knotpos_particles=60,
    knotpos_iters=80,
    λ=1e-6,
)

    nx, ny = size(Z)

    x = range(0.0, 1.0, length=nx)
    y = range(0.0, 1.0, length=ny)

    ϵ = eps(Float64)
    ξ = [clamp(xi, ϵ, 1 - ϵ) for xi in x]
    η = [clamp(yj, ϵ, 1 - ϵ) for yj in y]

    knot_max = max(p + 1, min(nx, ny) ÷ 2)

    best_v, _ = pso(
        v -> model_loss(v, Z, ξ, η, p),
        ndim=2,
        bounds=(Float64(p + 1), Float64(knot_max)),
        nparticles=knot_particles,
        maxiters=knot_iters
    )

    #nbx_opt = round(Int, best_v[1])
    #nby_opt = round(Int, best_v[2])

    nbx_opt = clamp(round(Int, best_v[1]), p + 1, knot_max)
    nby_opt = clamp(round(Int, best_v[2]), p + 1, knot_max)

    kx = clamped_knots_inclusive(0.0, 1.0, nbx_opt, p)
    ky = clamped_knots_inclusive(0.0, 1.0, nby_opt, p)

    Bx, By = basis_matrices(kx, ky, p, ξ, η)
    C = fit_bspline_ls(Bx, By, Z; λ=λ)

    ndim_knots = (nbx_opt - p - 1) + (nby_opt - p - 1)

    best_knots, _ = pso(
        v -> knot_loss(v, Z, ξ, η, p, nbx_opt, nby_opt),
        ndim=ndim_knots,
        bounds=(0.0, 1.0),
        nparticles=knotpos_particles,
        maxiters=knotpos_iters
    )

    mx = nbx_opt - p - 1
    kx_kopt = decode_knots(best_knots[1:mx], nbx_opt, p)
    ky_kopt = decode_knots(best_knots[mx+1:end], nby_opt, p)

    Bx_kopt, By_kopt = basis_matrices(kx_kopt, ky_kopt, p, ξ, η)
    ### Original
    #C_kopt = fit_bspline_ls(Bx_kopt, By_kopt, Z; λ=λ)
    ### Test feature 1/2
    #C_kopt = fit_bspline_robust(Bx_kopt, By_kopt, Z; λ=λ)
    ### Test feature 3/4
    C_kopt, mask = fit_bspline_mask_outliers(Bx_kopt, By_kopt, Z; λ=λ)

    Z_kopt = eval_surface(Bx_kopt, By_kopt, C_kopt)

    return Z_kopt, mask#, (kx_kopt, ky_kopt), C_kopt
end

end