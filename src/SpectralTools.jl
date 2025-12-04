module SpectralTools

using LinearAlgebra
using SparseArrays
using Kronecker

include("Utils.jl")
include("Chebyshev.jl")
include("Operators2D.jl")
include("Interpolation.jl")

export cheb_lobatto_nodes,
       cheb_D_matrices,
       grid2D,
       laplacian_plus_I_operator,
       apply_dirichlet!,
       solve_poisson_like,
       rebuild_solution,
       barycentric_interp_1d,
       interp2D_spectral,
       unvec

end
