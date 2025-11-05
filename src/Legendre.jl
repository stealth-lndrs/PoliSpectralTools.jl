module Legendre    # Início do módulo

# O que o módulo usa
using LinearAlgebra

# Funções públicas
export Legendre_Gauss_Basis, Legendre_Lobatto_Basis, Legendre_Radau_Basis, eval_legendre, GaussQuadTypes, D_Legendre 


# Base de Legendre nos nodos de Gauss, Lobatto e Radau
# 1. Legendre_Gauss_Basis
# 2. Legendre_Lobatto_Basis
# 3. Legendre_Radau_Basis
#   gauss_nodes_initial
#   eval_legendre



using LinearAlgebra # Necessário para SymTridiagonal e eigen


"""
    gauss_nodes_initial(N1::Int)

Calcula os N1 nós iniciais (xG) e os pesos (wG) de Gauss-Legendre 
usando a Matriz de Jacobi (Golub-Welsch).
"""
function gauss_nodes_initial(N1::Int)
    N = N1 - 1
    
    # Coeficientes 'beta' (Matriz de Jacobi sub/superdiagonal)
    # beta_k = 0.5 / sqrt(1 - (2k)^-2) para k=1:N
    k_vec = collect(1:N)
    beta = 0.5 ./ sqrt.(1.0 .- (2.0 .* k_vec).^(-2))
    
    # Matriz de Jacobi J (Tridiagonal Simétrica)
    J = SymTridiagonal(zeros(N1), beta)
    
    # Encontrar Autovalores e Autovetores
    decomp = eigen(J)
    
    # Os autovalores (λ) são os nós iniciais
    x = decomp.values
    
    # Índice para ordenação (o MATLAB faz sort(x) e usa o índice i)
    # O Julia's sortperm() retorna os índices da ordem crescente
    i = sortperm(x)
    xG = x[i]
    
    # Pesos (wG): w_i = 2 * V[1, i]^2
    # V[1, i] é a primeira linha dos autovetores, V.vectors
    w = 2.0 .* decomp.vectors[1, i].^2
    
    return (xG, w)
end



# 1. Nodos de Gauss

function Legendre_Gauss_Basis(n::Int)
    
    N1 = n + 1 # Número de pontos (roots of P_N1)
    N2 = n + 2 # Ponto de parada na recursão (grau N+1)
    
    (xG, wG_initial) = gauss_nodes_initial(N1)
    
    x0 = 2.0 * ones(N1) 
    epsilon = eps(Float64)
    
    L = zeros(N1, N2)
    Lp = zeros(N1) 
    
    while maximum(abs.(xG .- x0)) > epsilon
        x0 = copy(xG)
        
        L[:, 1] .= 1.0 
        L[:, 2] .= xG  
        
        for k = 2:N1 
            L[:, k+1] = (((2k - 1) .* xG .* L[:, k]) .- ((k - 1) .* L[:, k-1])) ./ k
        end
        
        # --- A CORREÇÃO ESTÁ AQUI ---
        # A derivada P'_{N1} usa o fator N1, não N2.
        Lp = (N1) .* (L[:, N1] .- xG .* L[:, N2]) ./ (1.0 .- xG.^2)
        # --------------------------------
        
        xG = x0 .- L[:, N2] ./ Lp
    end
    
    # --- 3. CONSTRUÇÃO DAS MATRIZES DE SAÍDA ---
    
    BL = L[:, 1:N1] 
    
    # Agora esta fórmula de peso (wG) está correta,
    # porque o 'Lp' que ela usa está correto.
    wG = 2.0 ./ ((1.0 .- xG.^2) .* Lp.^2) 
    
    C = Diagonal((2 .* (0:n) .+ 1) ./ 2.0)
    W = Diagonal(wG)
    BI = C * BL' * W
    
    return BL, BI, xG, wG
end

# 2. Legendre Lobatto 

"""
    Legendre_Lobatto_Basis(n::Int)

Calcula os nós (xLL), pesos (wGL) e a Matriz de Base (P) 
de Legendre-Gauss-Lobatto (LGL).

O método utiliza a relação de recorrência de Legendre 
e o refinamento de Newton-Raphson.

# Argumentos
- `n`: Grau máximo do polinômio (N = n+1 é o número de pontos).

# Retorno
- `P`: Matriz de Base/Vandermonde LGL (N x N).
- `xLL`: Nós de Gauss-Lobatto (Vetor N x 1).
- `wGL`: Pesos de Gauss-Lobatto (Vetor N x 1).
"""
function Legendre_Lobatto_Basis(n::Int)
    # Truncation + 1 (Número de Pontos)
    N1 = n + 1 
    
    # --- 1. Palpite Inicial e Matriz Base ---
    # Use os nós de Chebyshev-Gauss-Lobatto como palpite inicial
    # (0:n) * pi/n
    x = -cos.((0:n) .* π ./ n) # ]
    
    # Matriz de Base Legendre (Vandermonde)
    P = zeros(N1, N1) 

    # --- 2. Refinamento por Newton-Raphson ---
    
    xold = 2.0 * ones(N1) # Palpite inicial inválido para entrar no loop
    epsilon = eps(Float64) # Precisão de máquina (o 'eps' do MATLAB)

    # O laço de Newton-Raphson (continua enquanto a mudança for maior que 'eps')
    while maximum(abs.(x .- xold)) > epsilon
        xold = copy(x) # Salva o valor anterior

        # Preenche a Matriz P (Vandermonde) usando a Recorrência
        P[:, 1] .= 1.0 # P0(x) = 1
        P[:, 2] .= x   # P1(x) = x

        # Loop da Recorrência P_k+1(x) = [(2k-1)*x*P_k(x) - (k-1)*P_k-1(x)] / k
        for k = 2:n
            # Os pontos e a recorrência são vetorizados (Broadcasting)
            P[:, k+1] = (((2k - 1) .* x .* P[:, k]) .- ((k - 1) .* P[:, k-1])) ./ k
        end
        
        # Fórmula de Refinamento de Newton: x_new = x_old - f(x) / f'(x)
        # Onde f(x) = P_n(x) [do MATLAB, P(:,N1) é Pn(x)] e 
        # a derivada é P'_{n}(x) [do MATLAB, N1*P(:,N1)]
        
        # A expressão (x.*P(:,N1)-P(:,n)) / (N1*P(:,N1)) é a correção do Newton:
        # Corr_Term = f(x) / f'(x) [Note que MATLAB/von Winckel usam uma forma específica de P'n]

        # O termo de correção é calculado (vetorizado)
        Corr_Term = (x .* P[:, N1] .- P[:, n]) ./ (N1 .* P[:, N1])

        # Aplica a correção
        x = xold .- Corr_Term
        
        # O MATLAB lida com o caso x=1 e x=-1 na derivada P'n(x) automaticamente,
        # mas aqui confiamos na estabilidade do algoritmo de von Winckel.
    end

    # --- 3. Cálculo de Pesos (wGL) e Saída ---
    
    # Fórmula dos Pesos: wGL = 2 / [n*N1*P_n(x)^2]
    # O P(:,N1) já é o vetor Pn(x) no último passo do Newton.
    wGL = 2.0 ./ (n .* N1 .* P[:, N1].^2)

    xLL = x # Os nós refinados

    # A função original do MATLAB retorna (P, xLL, wGL)
    return P, xLL, wGL
end

# 3. Radau
"""
    Legendre_Radau_Basis(N::Int)

Calcula os nós (xG), pesos (wG), a Matriz de Base (BL) e a Inversa da Base (BI)
de Legendre-Gauss-Radau (LGR), com nó fixo em x = -1.

# Argumentos
- `N`: Número de pontos de quadratura (N = grau máximo + 1).

# Retorno
- `BL`: Matriz de Base (N x N).
- `BI`: Matriz Inversa de Base (N x N).
- `xG`: Nós de Gauss-Radau (Vetor N x 1).
- `wG`: Pesos de Gauss-Radau (Vetor N x 1).
"""



function Legendre_Radau_Basis(N::Int)
    
    N_total = N 
    
# --- 1. Chute Inicial (USANDO SUA SUGESTÃO: NÓS DE CHEBYSHEV-RADAU) ---
    k_indices = 0:N_total-1
    xG = -cos.((2 .* k_indices .* π) ./ (2 * N_total - 1))
    
    # --- 2. Refinamento por Newton-Raphson no Polinômio Radau ---
    
    x0 = 2.0 * ones(N_total)
    epsilon = eps(Float64)
    
    # O loop refina os N-1 nós livres (índices 2 a N), deixando o xG[1] = -1.0 fixo.
    while maximum(abs.(xG[2:end] .- x0[2:end])) > epsilon
        x0 = copy(xG)
        
        # O polinômio de Radau (LGR à esquerda) é: P_N(x) + P_{N-1}(x)
        
        # 1. Calcular P_N(x) e P_{N-1}(x) nos pontos xG
        P_N = eval_legendre(N_total, xG)
        P_N_minus_1 = eval_legendre(N_total - 1, xG)
        
        # 2. O Polinômio f(x) = P_N(x) + P_{N-1}(x)
        f_radau = P_N .+ P_N_minus_1
        
        # 3. A Derivada f'(x)
        P_N_deriv = eval_legendre(N_total, xG, true)[2]
        P_N_minus_1_deriv = eval_legendre(N_total - 1, xG, true)[2]
        f_radau_deriv = P_N_deriv .+ P_N_minus_1_deriv
        
        # 4. Refinamento de Newton-Raphson: Apenas para os N-1 nós livres (índices 2:N)
        
        # Corr_Term = f(x) / f'(x)
        Corr_Term = f_radau[2:end] ./ f_radau_deriv[2:end]
        
        # Aplica a correção APENAS aos nós livres
        xG[2:end] = x0[2:end] .- Corr_Term
        
        # Garante o nó fixo em -1.0
        xG[1] = -1.0 
    end
    
    # --- 3. CÁLCULO DE PESOS E MATRIZES ---
  # --- 3. CÁLCULO DE PESOS (wG) (A CORREÇÃO ESTÁ AQUI) ---
    
    # Recalcula P_{N-1} nos nós finais xG
    P_N_minus_1_final = eval_legendre(N_total - 1, xG)
    
    # Fórmula de Trefethen/Davis-Rabinowitz (Unificada):
    wG = (1.0 .- xG) ./ (N_total^2 .* P_N_minus_1_final.^2)
    
    # --- 4. CONSTRUÇÃO DA MATRIZ DE BASE BL E INVERSA BI ---
    
    # BL[i, j] = P_{j-1}(x_i)
    BL = Matrix{Float64}(undef, N_total, N_total) 
    
    for k = 0:N_total-1 # k é o grau do polinômio, de 0 até N-1
        BL[:, k + 1] = eval_legendre(k, xG) 
    end
    
    # BI (Matriz Inversa) - Usando a fórmula de projeção com pesos.
    C = Diagonal((2 .* (0:N_total-1) .+ 1) ./ 2.0)
    W = Diagonal(wG)
    BI = C * BL' * W

    return BL, BI, xG, wG

end


"""
    eval_legendre(n, x, deriv=false)

Calcula Pn(x) e P'n(x) (se deriv=true) usando a recursão de três termos.
O argumento 'x' pode ser um escalar ou um vetor.
"""
function eval_legendre(n, x, deriv::Bool=false)
    # Garante que x seja Float64 e trata o caso escalar para vetorizar
    X = x isa Number ? [x] : vec(Float64.(x))
    
    # Inicializa P_n-1 e P_n-2
    P_prev2 = ones(length(X)) # P0(x) = 1
    P_prev1 = copy(X)         # P1(x) = x

    if n == 0
        P_n = P_prev2
    elseif n == 1
        P_n = P_prev1
    else
        P_n = P_prev1 # Inicialização para o loop
        # Loop da recorrência P_k+1(x) = [(2k+1)x * P_k(x) - k * P_k-1(x)] / (k+1)
        for k = 1:n-1
            P_n = (((2k + 1) .* X .* P_prev1) .- (k .* P_prev2)) ./ (k + 1)
            P_prev2 = P_prev1
            P_prev1 = P_n
        end
    end

    if deriv
        # Usa a identidade (1-x^2)P'_n = n(P_{n-1} - x*P_n)
        P_deriv = (n .* (P_prev2 .- X .* P_n)) ./ (1.0 .- X.^2)

        # Trata os endpoints (x = +/- 1) onde a divisão por zero ocorreria
        for i in 1:length(X)
            if isapprox(X[i], 1.0, atol=1e-15)
                P_deriv[i] = n * (n + 1) / 2.0
            elseif isapprox(X[i], -1.0, atol=1e-15)
                P_deriv[i] = ((-1.0)^(n+1)) * n * (n + 1) / 2.0
            end
        end
        return P_n, P_deriv
    else
        return P_n
    end
end


"""
tipo: 
1. Gauss           2. Lobatto   
3. Radau [a,b)     4. Radau (a,b]


Example:  a = -1; b = 1; n = 16
z,w = GaussQuadTypes(a, b,n, 2)

"""



"""
    eval_legendre_p(n, x)

Avalia o polinômio de Legendre P_n(x) em x usando a
recorrência de 3 termos. É numericamente estável.
"""
function eval_legendre_p(n::Int, x::Number)
    if n == 0
        return 1.0
    elseif n == 1
        return x
    end
    
    # Inicializa P_0(x) e P_1(x)
    z0 = 1.0
    z1 = x
    
    # z2 irá guardar P_j(x)
    z2 = 0.0
    
    # Itera de j=1 até n-1 para encontrar P_n(x)
    for j = 1:n-1
        z2 = x * z1 * (2*j + 1) / (j + 1) - z0 * j / (j + 1)
        z0 = z1
        z1 = z2
    end
    # z1 agora contém P_n(x)
    return z1
end

# Versão "broadcasted" para aplicar a um vetor de pontos `x_vec`
eval_legendre_p(n::Int, x_vec::AbstractVector) = eval_legendre_p.(n, x_vec)


"""
    GaussQuadTypes(a, b, n, tipo=1)

Calcula os nós (z) e pesos (w) para quadratura de Gauss.
(Versão 2, com correção para Radau)

# Argumentos
- `a`, `b`: Limites do intervalo de integração.
- `n`: Número de pontos de quadratura.
- `tipo`: 1 (Gauss-Legendre), 2 (Gauss-Lobatto), 3 (Gauss-Radau, nó fixo em -1).
"""
function GaussQuadTypes(a::Real, b::Real, n::Int, tipo::Int = 1)
    
    local z, w

    if tipo == 1
        # --- Caso 1: Gauss-Legendre (Correto) ---
        i = 1:n-1
        beta = @. 0.5 * (1 - (2*i)^(-2))^(-0.5)
        alpha = zeros(n)
        t = Tridiagonal(beta, alpha, beta)
        
        F = eigen(t)
        z = F.values
        V = F.vectors
        w = @. 2 * (V[1, :])^2
        
    elseif tipo == 2
        # --- Caso 2: Gauss-Lobatto (Correto) ---
        N = n - 1
        w = zeros(n)
        z = zeros(n)
        z[1] = -1.0
        z[end] = 1.0
        
        if N > 1
            if N == 2
                z[2] = 0.0
            else
                j = 1:N-2
                beta = @. 0.5 * sqrt((j * (j + 2)) / ((j + 0.5) * (j + 1.5)))
                alpha = zeros(N - 1)
                M = Tridiagonal(beta, alpha, beta)
                z[2:N] = eigen(M).values
            end
        end
        
        w[1] = 2.0 / (N * n)
        w[end] = w[1]
        
        # Pesos para os nós internos
        # Usamos a nova função `eval_legendre_p`
        # P_N(x) é P_{n-1}(x)
        y = eval_legendre_p(n - 1, z[2:N])
        w[2:N] = @. 2.0 / (N * n * y^2)
        
    elseif tipo == 3
        # --- Caso 3: Gauss-Radau (CORRIGIDO) ---
        
        # O cálculo dos nós (abscissas z) estava correto
        j = 0:n-2
        alpha = @. 1.0 / ((2*j + 1) * (2*j + 3))
        j = 1:n-2
        beta = @. sqrt(j * (j + 1)) / (2*j + 1)
        
        A = Tridiagonal(beta, alpha, beta)
        z_interior = eigen(A).values
        z = [-1.0; z_interior]
        
        # --- CORREÇÃO NO CÁLCULO DOS PESOS ---
        # Substituímos o problezmático `LegPoly` + `polyval` pela
        # função `eval_legendre_p`, que é estável.
        
        # A fórmula original (correta) depende de P_{n-1}(z)
        y = eval_legendre_p(n - 1, z)
        
        # Calcula os pesos com a fórmula
        w = @. (1 - z) / (n^2 * y^2)
    elseif tipo == 4
        # --- Caso 3: Gauss-Radau (CORRIGIDO) ---
        
        # O cálculo dos nós (abscissas z) estava correto
        j = 0:n-2
        alpha = @. 1.0 / ((2*j + 1) * (2*j + 3))
        j = 1:n-2
        beta = @. sqrt(j * (j + 1)) / (2*j + 1)
        
        A = Tridiagonal(beta, alpha, beta)
        z_interior = eigen(A).values
        z = [-1.0; z_interior]
        z = -reverse(z)
        
        # --- CORREÇÃO NO CÁLCULO DOS PESOS ---
        # Substituímos o problezmático `LegPoly` + `polyval` pela
        # função `eval_legendre_p`, que é estável.
        
        # A fórmula original (correta) depende de P_{n-1}(z)
        y = eval_legendre_p(n - 1, z)
        
        # Calcula os pesos com a fórmula
        w = @. (1 + z) / (n^2 * y^2)
        
        
        
        
    else
        error("Tipo de quadratura inválido: $tipo. Use 1 (Gauss), 2 (Lobatto) ou 3 (Radau).")
    end
    
    # --- Etapa Final: Mapear de [-1, 1] para [a, b] ---
    z_scaled = @. ((b - a) * z + b + a) / 2
    w_scaled = @. w * (b - a) / 2
    
    return (z_scaled, w_scaled)
end


# ====================================================================
# 1. Matriz de Diferenciação Modal de LEGENDRE (D_Legendre)
#    Aplica-se aos coeficientes a[0] a a[N-1]
#    Dimensão: N x N
# ====================================================================
"""
    D_Legendre(N::Int)

Cria a matriz de diferenciação modal de Legendre de dimensão N x N.
A matriz D, de elementos D[i, j], transforma o vetor de coeficientes
a[j] no vetor de coeficientes da derivada a'[i].
Indices (i, j) vão de 1 a N, correspondendo aos índices modais 0 a N-1.
"""
function D_Legendre(n::Int)
    N = n+1;
    D = zeros(Float64, N, N)

    # i e j são os índices baseados em 1 (1 a N).
    # O índice modal é i_modal = i - 1 e j_modal = j - 1.
    for i in 1:N
        for j in (i + 1):N
            # Condição: i_modal + j_modal deve ser ímpar
            # (i - 1) + (j - 1) ímpar <=> i + j par.
            # O critério de paridade é o mesmo para i+j ou i_modal+j_modal.
            if iseven(i + j) 
                # O valor do elemento é 2 * i_modal + 1 = 2*(i-1) + 1 = 2i - 1
                D[i, j] = 2i - 1
            end
        end
    end
    return D
end

end