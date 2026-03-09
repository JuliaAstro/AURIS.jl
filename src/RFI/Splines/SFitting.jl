module SFitting

using LinearAlgebra, Statistics
using ..BProbability: bayesian_prob_good

function fit_bspline_ls(Bx, By, Z; λ=1e-6)
    A = kron(By, Bx)
    z = vec(Z)
    c = (A' * A + λ * I(size(A, 2))) \ (A' * z)
    reshape(c, size(Bx, 2), size(By, 2))
end

function fit_bspline_ls_weighted(Bx, By, Z, w; λ=1e-6)
    A = kron(By, Bx)
    z = vec(Z)
    wv = vec(w)

    W = Diagonal(wv)

    c = (A' * W * A + λ * I(size(A, 2))) \ (A' * W * z)
    reshape(c, size(Bx, 2), size(By, 2))
end

function fit_bspline_robust(Bx, By, Z;
    λ=1e-6,
    β=0.01,
    γ=10.0,
    maxiters=10,
    tol=1e-4)

    nx, ny = size(Z)
    w = ones(nx, ny)

    C = fit_bspline_ls(Bx, By, Z; λ=λ)

    for _ in 1:maxiters
        Ẑ = eval_surface(Bx, By, C)
        resid = Z .- Ẑ

        med_r = median(resid)
        σ = 1.4826 * median(abs.(resid .- med_r))

        w_new = bayesian_prob_good(resid .- med_r, σ; β=β, γ=γ)
        w_new[resid.≤0] .= 1.0

        C_new = fit_bspline_ls_weighted(Bx, By, Z, w_new; λ=λ)

        if norm(C_new - C) < tol
            C = C_new
            break
        end

        C = C_new
        w = w_new
    end

    C
end

function fit_bspline_ls_masked(Bx, By, Z, mask; λ=1e-6)
    A = kron(By, Bx)
    z = vec(Z)
    mv = vec(mask)

    A_masked = A[mv, :]
    z_masked = z[mv]

    c = (A_masked' * A_masked + λ * I(size(A, 2))) \ (A_masked' * z_masked)

    reshape(c, size(Bx, 2), size(By, 2))
end

function fit_bspline_mask_outliers(
    Bx, By, Z;
    λ=1e-6,
    β=0.01,
    γ=10.0,
    prob_threshold=0.5,
    maxiters=10,
    C_init=nothing
)

    nx, ny = size(Z)
    mask = trues(nx, ny)

    C = isnothing(C_init) ? fit_bspline_ls(Bx, By, Z; λ=λ) : C_init

    resid_init = vec(Z .- eval_surface(Bx, By, C))
    med_init = median(resid_init)
    σ_init = 1.4826 * median(abs.(resid_init .- med_init))
    tail_frac = mean(abs.(resid_init .- med_init) .> 3.0 * σ_init)
    β_eff = clamp(tail_frac, β / 10, 0.2)

    for _ in 1:maxiters
        Ẑ = eval_surface(Bx, By, C)
        resid = Z .- Ẑ

        resid_in = resid[mask]
        med_r = median(resid_in)
        σ = 1.4826 * median(abs.(resid_in .- med_r))

        p_good = bayesian_prob_good(resid .- med_r, σ; β=β_eff, γ=γ)
        p_good[resid.≤0] .= 1.0

        new_mask = p_good .> prob_threshold

        if new_mask == mask
            break
        end

        mask = new_mask
        C = fit_bspline_ls_masked(Bx, By, Z, mask; λ=λ)
    end

    return C, mask
end


function eval_surface(Bx, By, C)
    Bx * C * By'
end

end
