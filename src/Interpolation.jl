"""
    barycentric_interp_1d(x_nodes, f_nodes, x_new; weights=nothing)

Evaluate the barycentric Lagrange interpolant defined by the nodes
`x_nodes` and values `f_nodes` at the query points `x_new`.
Optional pre-computed barycentric `weights` can be provided to avoid
recomputing them across calls.
"""
function barycentric_interp_1d(x_nodes::AbstractVector,
                               f_nodes::AbstractVector,
                               x_new::AbstractVector;
                               weights=nothing)
    length(x_nodes) == length(f_nodes) ||
        throw(ArgumentError("x_nodes and f_nodes must have the same length"))
    nodes = x_nodes isa Vector{Float64} ? x_nodes : Vector{Float64}(x_nodes)
    vals = Vector{Float64}(f_nodes)
    queries = x_new isa Vector{Float64} ? x_new : Vector{Float64}(x_new)
    w = weights === nothing ? barycentric_weights(nodes) :
        (weights isa Vector{Float64} ? weights : Vector{Float64}(weights))
    length(w) == length(nodes) ||
        throw(ArgumentError("weights must match the number of nodes"))
    result = Vector{Float64}(undef, length(queries))
    tol = 50 * eps(Float64)
    for (iq, xq) in enumerate(queries)
        numerator = 0.0
        denominator = 0.0
        exact = false
        exact_val = 0.0
        for j in eachindex(nodes)
            diff = xq - nodes[j]
            if abs(diff) <= tol * max(1.0, abs(nodes[j]))
                exact = true
                exact_val = vals[j]
                break
            end
            wj = w[j] / diff
            numerator += wj * vals[j]
            denominator += wj
        end
        result[iq] = exact ? exact_val : numerator / denominator
    end
    return result
end

barycentric_interp_1d(x_nodes::AbstractVector,
                      f_nodes::AbstractVector,
                      x_new::Number;
                      weights=nothing) =
    barycentric_interp_1d(x_nodes, f_nodes, [x_new]; weights=weights)[1]

"""
    interp2D_spectral(x_nodes, y_nodes, F, x_new, y_new)

Interpolate the 2D field `F` defined on the tensor grid given by
`x_nodes` and `y_nodes` onto the new tensor grid `x_new × y_new` using
successive barycentric interpolations.
"""
function interp2D_spectral(x_nodes::AbstractVector,
                           y_nodes::AbstractVector,
                           F::AbstractMatrix,
                           x_new::AbstractVector,
                           y_new::AbstractVector)
    size(F, 1) == length(y_nodes) ||
        throw(ArgumentError("First dimension of F must match length(y_nodes)"))
    size(F, 2) == length(x_nodes) ||
        throw(ArgumentError("Second dimension of F must match length(x_nodes)"))
    xq = Vector{Float64}(x_new)
    yq = Vector{Float64}(y_new)
    wx = barycentric_weights(Vector{Float64}(x_nodes))
    wy = barycentric_weights(Vector{Float64}(y_nodes))
    tmp = Array{Float64}(undef, length(y_nodes), length(xq))
    for j in 1:length(y_nodes)
        row = view(F, j, :)
        tmp[j, :] = barycentric_interp_1d(x_nodes, row, xq; weights=wx)
    end
    result = Array{Float64}(undef, length(yq), length(xq))
    for k in 1:length(xq)
        column = view(tmp, :, k)
        result[:, k] = barycentric_interp_1d(y_nodes, column, yq; weights=wy)
    end
    return result
end

function barycentric_weights(nodes::AbstractVector)
    n = length(nodes)
    w = ones(Float64, n)
    for j in 1:n
        prod_val = 1.0
        xj = nodes[j]
        for k in 1:n
            k == j && continue
            prod_val *= (xj - nodes[k])
        end
        w[j] = 1.0 / prod_val
    end
    return w
end
