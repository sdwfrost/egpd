# Neural network architecture for bivariate BEGPD estimation
# Loaded by egpd R package during .init_julia_begpd()
#
# Defines initializenetwork(n, d) which builds a DeepSet network with
# ResidualBlocks, LayerNorm, and expert summary statistics.

import NeuralEstimators: ResidualBlock

function initializenetwork(n::Integer, d::Integer; residuals::Bool = true)

    # Expert summary statistics: log sample size (centered)
    S(Z) = log.(samplesize(Z)) .- 7.7f0

    # Width of each hidden layer
    w = 128

    # Input dimension to phi (includes expert summary statistic)
    w_phi = w + !isnothing(S)

    if residuals
        # Architecture using residual connections and layer normalization
        function ResidualBlock(dim::Integer)
            layer = Chain(Dense(dim, dim), LayerNorm(dim), relu)
            return SkipConnection(layer, +)
        end

        psi = Chain(
            Dense(n, w), LayerNorm(w), relu,
            ResidualBlock(w),
            ResidualBlock(w)
        )

        phi = Chain(
            ResidualBlock(w_phi),
            ResidualBlock(w_phi),
            LayerNorm(w_phi), relu,
            Dense(w_phi, d)
        )
    else
        # Simple architecture without residual connections
        psi = Chain(Dense(n, w, relu), Dense(w, w, relu), Dense(w, w, relu))
        phi = Chain(Dense(w_phi, w, relu), Dense(w, w, relu), Dense(w, d))
    end

    network = DeepSet(psi, phi; S = S)

    return network
end
