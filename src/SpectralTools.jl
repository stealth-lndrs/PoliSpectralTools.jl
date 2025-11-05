module SpectralTools

# 1. Inclui o código dos seus arquivos
include("Chebyshev.jl")
include("Legendre.jl")
include("Fourier.jl")
include("PseudoGen.jl")



# 2. Traz as funções dos submódulos para o escopo de "SpectralTools"
# O ponto (.) antes do nome é crucial, significa "deste módulo"
using .Chebyshev
using .Legendre
using .Fourier
using .PseudoGen

# 3. Exporta as funções que você quer que o *usuário final* veja
# (Isso permite ao usuário escrever `D_Cheby` em vez de `SpectralTools.Chebyshev.D_Cheby`)

export DCT_Chb_Gauss, DCT_Chb_Lob, D_Cheby, ChebyQuadrature
export Legendre_Gauss_Basis, Legendre_Lobatto_Basis, Legendre_Radau_Basis, eval_legendre, GaussQuadTypes, D_Legendre
export D_Fourier, Multi_Diff_Mat, Bary_Trig_Mat, I_Fourier, fourier_quad
export Bary_Interp, Generalized_Diff_Mat


greet() = print("Hello World!")

end # module SpectralTools
