"""
    cheb_lobatto_nodes(N)

Return the `N + 1` Chebyshev–Lobatto points on ``[-1, 1]`` in descending order.
"""
function cheb_lobatto_nodes(N::Integer)
    N < 1 && throw(ArgumentError("N must be at least 1"))
    nodes = cos.(pi .* (0:N) ./ N)
    return Vector{Float64}(nodes)
end

"""
    cheb_D_matrices(N)

Construct the first- and second-order Chebyshev spectral differentiation
matrices for ``N + 1`` Chebyshev–Lobatto nodes.
"""
function cheb_D_matrices(N::Integer)
    x = cheb_lobatto_nodes(N)
    n = length(x)
    c = ones(Float64, n)
    c[1] = 2.0
    c[end] = 2.0
    c .*= (-1.0).^(0:n-1)

    xcol = reshape(x, n, 1)
    X = repeat(xcol, 1, n)
    dX = X .- X'
    D = (c * (1.0 ./ c)') ./ (dX + I)
    D .-= Diagonal(vec(sum(D, dims = 2)))
    D2 = D * D
    return Matrix{Float64}(D), Matrix{Float64}(D2)
end
