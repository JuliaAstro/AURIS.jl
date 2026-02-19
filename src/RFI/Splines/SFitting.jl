module SFitting

using LinearAlgebra, Statistics
using ..BProbability: bayesian_prob_good

function fit_bspline_ls(Bx, By, Z; λ=1e-6)
    A = kron(By, Bx)
    z = vec(Z)
    c = (A' * A + λ * I(size(A, 2))) \ (A' * z)
    reshape(c, size(Bx, 2), size(By, 2))
end

### Test feature
function fit_bspline_ls_weighted(Bx, By, Z, w; λ=1e-6)
    A = kron(By, Bx)
    z = vec(Z)
    wv = vec(w)

    W = Diagonal(wv)

    c = (A' * W * A + λ * I(size(A, 2))) \ (A' * W * z)
    reshape(c, size(Bx, 2), size(By, 2))
end

### Test feature 2
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
        Ẑ = eval_surface(Bx, By, C)
        resid = Z .- Ẑ

        σ = 1.4826 * median(abs.(resid))

        w_new = bayesian_prob_good(resid, σ; β=β, γ=γ)

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

### Test feature 3
function fit_bspline_ls_masked(Bx, By, Z, mask; λ=1e-6)
    A = kron(By, Bx)
    z = vec(Z)
    mv = vec(mask)

    # Keep only inliers
    A_masked = A[mv, :]
    z_masked = z[mv]

    c = (A_masked' * A_masked + λ * I(size(A, 2))) \ (A_masked' * z_masked)

    reshape(c, size(Bx, 2), size(By, 2))
end
### Test feature 4
function fit_bspline_mask_outliers(
    Bx, By, Z;
    λ=1e-6,
    β=0.01,
    γ=10.0,
    prob_threshold=0.5,
    maxiters=10
)

    nx, ny = size(Z)
    mask = trues(nx, ny)

    C = fit_bspline_ls(Bx, By, Z; λ=λ)

    for _ in 1:maxiters
        Ẑ = eval_surface(Bx, By, C)
        resid = Z .- Ẑ

        # robust sigma estimate
        σ = 1.4826 * median(abs.(resid[mask]))

        p_good = bayesian_prob_good(resid, σ; β=β, γ=γ)

        new_mask = p_good .> prob_threshold

        # Stop if mask stable
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