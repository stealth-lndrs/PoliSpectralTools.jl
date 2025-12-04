"""
    grid2D(x, y)

Construct a tensor-product grid from the 1D nodes `x` and `y`.
Returns matrices `(X, Y)` where `X[i, j]` stores the ``x``-coordinate and
`Y[i, j]` the ``y``-coordinate of the grid point `(x_j, y_i)`.
"""
function grid2D(x::AbstractVector, y::AbstractVector)
    nx = length(x)
    ny = length(y)
    Tx = promote_type(eltype(x), eltype(y))
    xrow = reshape(Tx.(x), 1, nx)
    ycol = reshape(Tx.(y), ny, 1)
    X = repeat(xrow, ny, 1)
    Y = repeat(ycol, 1, nx)
    return X, Y
end

"""
    laplacian_plus_I_operator(D2x, D2y)

Assemble the operator ``\nabla^2 + I`` on a tensor-product grid given the
1D second-derivative matrices `D2x` and `D2y`.
"""
function laplacian_plus_I_operator(D2x::AbstractMatrix, D2y::AbstractMatrix)
    size(D2x, 1) == size(D2x, 2) || throw(ArgumentError("D2x must be square"))
    size(D2y, 1) == size(D2y, 2) || throw(ArgumentError("D2y must be square"))
    nx = size(D2x, 1)
    ny = size(D2y, 1)
    Ix = Matrix{Float64}(I, nx, nx)
    Iy = Matrix{Float64}(I, ny, ny)
    A = kron(Ix, D2y) + kron(D2x, Iy) + I(nx * ny)
    return sparse(A)
end

"""
    apply_dirichlet!(A, b, mask_internal, boundary_values)

Apply Dirichlet boundary conditions to the linear system defined by `A` and
`b`. The Boolean matrix `mask_internal` marks the interior points while
`boundary_values` holds the known solution along the boundary.
The function returns `(A_reduced, b_reduced, boundary_vec)`.
"""
function apply_dirichlet!(A::AbstractMatrix,
                          b::AbstractVector,
                          mask_internal::AbstractMatrix{Bool},
                          boundary_values::AbstractMatrix)
    ensure_same_shape("mask_internal", mask_internal, "boundary_values", boundary_values)
    total_pts = length(mask_internal)
    length(b) == total_pts ||
        throw(ArgumentError("Length of b must match the number of grid points"))
    mask_vec, internal_idx, boundary_idx = flatten_mask(mask_internal)
    boundary_vec = vec(boundary_values)
    b_internal = b[internal_idx]
    if !isempty(boundary_idx)
        b_internal .-= A[internal_idx, boundary_idx] * boundary_vec[boundary_idx]
    end
    A_internal = A[internal_idx, internal_idx]
    return sparse(A_internal), b_internal, boundary_vec
end

"""
    solve_poisson_like(A, b)

Solve the linear system defined by the operator `A` and the right-hand side
`b` using the default Julia linear solver.
"""
solve_poisson_like(A::AbstractMatrix, b::AbstractVector) = A \ b

"""
    rebuild_solution(F_internal, boundary_values, mask)

Rebuild the full 2D field inserting the computed interior solution
`F_internal` inside `boundary_values` according to the Boolean `mask`.
"""
function rebuild_solution(F_internal::AbstractVector,
                          boundary_values::AbstractMatrix,
                          mask::AbstractMatrix{Bool})
    ensure_same_shape("mask", mask, "boundary_values", boundary_values)
    count(identity, mask) == length(F_internal) ||
        throw(ArgumentError("Number of interior values does not match mask"))
    result = copy(boundary_values)
    result[mask] .= F_internal
    return result
end
