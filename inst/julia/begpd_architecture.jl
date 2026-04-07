# Neural network architectures for neural-Bayes EGPD estimation
# Loaded by egpd R package during .init_julia_begpd()

function _begpd_expert_summary(Z)
    log.(samplesize(Z)) .- 7.7f0
end

function _begpd_residual_block(dim::Integer)
    layer = Chain(Dense(dim, dim), LayerNorm(dim), relu)
    SkipConnection(layer, +)
end

function initialize_point_network(n::Integer, d::Integer; residuals::Bool = true)
    S = _begpd_expert_summary
    w = 128
    w_phi = w + 1

    if residuals
        psi = Chain(
            Dense(n, w), LayerNorm(w), relu,
            _begpd_residual_block(w),
            _begpd_residual_block(w)
        )

        phi = Chain(
            _begpd_residual_block(w_phi),
            _begpd_residual_block(w_phi),
            LayerNorm(w_phi), relu,
            Dense(w_phi, d)
        )
    else
        psi = Chain(Dense(n, w, relu), Dense(w, w, relu), Dense(w, w, relu))
        phi = Chain(Dense(w_phi, w, relu), Dense(w, w, relu), Dense(w, d))
    end

    DeepSet(psi, phi; S = S)
end

function initialize_ratio_network(n::Integer, d::Integer; residuals::Bool = true)
    S = _begpd_expert_summary
    w = 128
    w_phi = w + 1 + d

    if residuals
        psi = Chain(
            Dense(n, w), LayerNorm(w), relu,
            _begpd_residual_block(w),
            _begpd_residual_block(w)
        )

        phi = Chain(
            _begpd_residual_block(w_phi),
            _begpd_residual_block(w_phi),
            LayerNorm(w_phi), relu,
            Dense(w_phi, 1)
        )
    else
        psi = Chain(Dense(n, w, relu), Dense(w, w, relu), Dense(w, w, relu))
        phi = Chain(Dense(w_phi, w, relu), Dense(w, w, relu), Dense(w, 1))
    end

    DeepSet(psi, phi; S = S)
end

# Backwards-compatible alias used by the training helpers for point estimators.
initializenetwork(n::Integer, d::Integer; residuals::Bool = true) =
    initialize_point_network(n, d; residuals = residuals)
