# Neural network architectures and cross-version helpers for neural-Bayes EGPD estimation
# Loaded by egpd R package during .init_julia_begpd() (after `using NeuralEstimators, Flux, Random, BSON`)
# Executed with NeuralEstimators.jl v0.1.0 (Flux 0.14); RatioEstimator call methods of v0.1.3/v0.1.4
# checked verbatim; written against the v0.2.1+ and v0.3.x APIs (the `grid` keyword of
# logratio() first appears in v0.2.1; training through R NeuralEstimators needs v0.3+)

using Statistics: mean

function _begpd_expert_summary(Z)
    log.(samplesize(Z)) .- 7.7f0
end

function _begpd_residual_block(dim::Integer)
    layer = Chain(Dense(dim, dim), LayerNorm(dim), relu)
    SkipConnection(layer, +)
end

function _make_deepset(psi, phi; S = _begpd_expert_summary)
    if hasmethod(NeuralEstimators.DeepSet, Tuple{typeof(psi), typeof(phi), typeof(mean)}, (:condition_on_sample_size,))
        return NeuralEstimators.DeepSet(psi, phi, mean; condition_on_sample_size = !isnothing(S))
    else
        return NeuralEstimators.DeepSet(psi, phi, mean; S = S)
    end
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

    _make_deepset(psi, phi; S = S)
end

function initialize_ratio_network(n::Integer, d::Integer; residuals::Bool = true)
    summary = initialize_point_network(n, 128; residuals = residuals)
    if hasmethod(NeuralEstimators.RatioEstimator, Tuple{typeof(d), typeof(summary)}, (:num_summaries,))
        return NeuralEstimators.RatioEstimator(d, summary; num_summaries = 128)
    elseif hasmethod(NeuralEstimators.RatioEstimator, Tuple{typeof(summary), typeof(d)}, (:num_summaries,))
        return NeuralEstimators.RatioEstimator(summary, d; num_summaries = 128)
    else
        w = 128
        w_phi = w + 1
        w_in = w_phi + d
        if residuals
            psi = Chain(
                Dense(n, w), LayerNorm(w), relu,
                _begpd_residual_block(w),
                _begpd_residual_block(w)
            )
            phi = Chain(
                Dense(w_in, w), LayerNorm(w), relu,
                _begpd_residual_block(w),
                _begpd_residual_block(w),
                Dense(w, 1)
            )
        else
            psi = Chain(Dense(n, w, relu), Dense(w, w, relu), Dense(w, w, relu))
            phi = Chain(Dense(w_in, w, relu), Dense(w, w, relu), Dense(w, 1))
        end
        deepset = _make_deepset(psi, phi; S = _begpd_expert_summary)
        return NeuralEstimators.RatioEstimator(deepset)
    end
end

function egpd_estimate(estimator, Z; batchsize::Integer = 32, use_gpu::Bool = false)
    if isdefined(NeuralEstimators, :estimate)
        return NeuralEstimators.estimate(estimator, Z; batchsize = Int(batchsize), use_gpu = use_gpu)
    elseif isdefined(NeuralEstimators, :estimateinbatches)
        return NeuralEstimators.estimateinbatches(estimator, Z; batchsize = Int(batchsize), use_gpu = use_gpu)
    else
        return estimator(Float32.(Z))
    end
end

# Log likelihood-to-evidence ratio log r(Z, theta) for one data set Z and each
# column of theta, evaluated in batches of columns. Uses the network output on
# the logit scale (logratio() in v0.2.1+, the inner network of the v0.1.x
# RatioEstimator otherwise) so that saturated classifier probabilities never
# round to log r = +/-Inf.
function egpd_logratio(estimator, Z, theta::AbstractMatrix; batchsize::Integer = 4096)
    Zf = Float32.(Z)
    G = size(theta, 2)
    out = Vector{Float64}(undef, G)
    for first in 1:Int(batchsize):G
        last = min(first + Int(batchsize) - 1, G)
        out[first:last] = _egpd_logratio_batch(estimator, Zf, Float32.(theta[:, first:last]))
    end
    return out
end

function _egpd_logratio_batch(estimator, Z, theta)
    if isdefined(NeuralEstimators, :logratio)
        # v0.2.1+: K x G matrix for K = 1 data set
        lr = NeuralEstimators.logratio(estimator, [Z]; grid = theta)
    else
        # v0.1.x: the wrapped network outputs logit(c) = log r before any
        # sigmoid; the field is `deepset` in v0.1.0 and `network` in v0.1.3+
        # (whose call methods return sigmoid(.) or, in v0.1.4, the logit).
        net = if hasproperty(estimator, :deepset)
            estimator.deepset
        elseif hasproperty(estimator, :network)
            estimator.network
        else
            error("Unsupported NeuralEstimators.jl RatioEstimator: no logratio() ",
                  "and no `deepset` or `network` field")
        end
        lr = net((Z, theta))
        if lr isa AbstractVector && !(eltype(lr) <: Number)
            lr = reduce(vcat, vec.(lr))
        end
    end
    v = vec(Float64.(lr))
    length(v) == size(theta, 2) || error("log-ratio has length $(length(v)), expected $(size(theta, 2))")
    return v
end

# Flux state saved/loaded with BSON under the key `model_state`, the format
# written by NeuralEstimators::savestate() for Flux estimators.
function egpd_savestate(estimator, filename::AbstractString)
    BSON.bson(filename, Dict(:model_state => Flux.state(estimator)))
    return filename
end

function egpd_loadstate(estimator, filename::AbstractString)
    model_state = BSON.load(filename, @__MODULE__)[:model_state]
    Flux.loadmodel!(estimator, model_state)
    return estimator
end

function egpd_backend_versions()
    v(m) = try
        string(pkgversion(m))
    catch
        "unknown"
    end
    return [string(VERSION), v(NeuralEstimators), v(Flux)]
end

# Backwards-compatible alias used by the training helpers for point estimators.
initializenetwork(n::Integer, d::Integer; residuals::Bool = true) =
    initialize_point_network(n, d; residuals = residuals)
