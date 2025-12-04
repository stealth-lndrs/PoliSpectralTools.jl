"""
    unvec(v, m, n)

Reshape the vector `v` into an `m × n` matrix assuming column-major ordering.
"""
function unvec(v::AbstractVector, m::Integer, n::Integer)
    length(v) == m * n || throw(ArgumentError("vector length does not match target dimensions"))
    return reshape(v, m, n)
end

"""
    ensure_same_shape(name_a, A, name_b, B)

Throw an error if `A` and `B` do not share the same size.
"""
function ensure_same_shape(name_a::AbstractString, A::AbstractArray, name_b::AbstractString, B::AbstractArray)
    size(A) == size(B) ||
        throw(ArgumentError("$name_a and $name_b must have identical dimensions"))
    return nothing
end

"""
    flatten_mask(mask)

Return a tuple `(mask_vec, internal_idx, boundary_idx)` for a Boolean `mask`.
"""
function flatten_mask(mask::AbstractArray{Bool})
    mask_vec = vec(mask)
    internal_idx = findall(mask_vec)
    boundary_idx = findall(x -> !x, mask_vec)
    return mask_vec, internal_idx, boundary_idx
end
