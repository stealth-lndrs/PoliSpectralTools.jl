module Chebyshev  # Início do módulo

using LinearAlgebra
using FFTW

export DCT_Chb_Gauss, DCT_Chb_Lob, D_Cheby, ChebyQuadrature  


"""
    DCT_Chb_Lob(ys, direction)

Transformada de Chebyshev-Lobatto (DCT-I) para pontos de Cheby-Lobatto. 
Calcula a transformação nodal (física) para espectral e vice-versa.

# Argumentos
- `ys`: Vetor coluna (nodal ou espectral).
- `direction`: 1 (Nodal -> Espectral) ou 2 (Espectral -> Nodal).

# Retorno
- `y`: O vetor ou matriz transformado(a).
"""
function DCT_Chb_Lob(ys, direction::Int)
    ys = Float64.(ys)
    N = length(ys)
    
    if N == 1
        return ys
    end
    
    n = N - 1
    
    # --- Transformada Direta (direction == 1): (Mantida correta) ---
    if direction == 1 
        yp = reverse(ys)
        yh = FFTW.r2r(yp, FFTW.REDFT00, 1)
        
        # Normalização: Fator 1/n para interior, 1/2n para extremidades
        yh[2:n] .= yh[2:n] ./ n 
        yh[[1, N]] .= yh[[1, N]] ./ (2 * n)
        
        return yh
        
    # --- Transformada Inversa (direction == 2): ORDEM CORRIGIDA ---
    elseif direction == 2
        # 1. Aplicar a REVERSÃO da Normalização ANTES da FFTW
        
        # Criamos uma cópia do vetor de coeficientes para pré-escalar
        y_scaled = copy(ys)
        
        # Termos internos (c_1 até c_{n-1}): Multiplicar por n (Inverso de 1/n)
        y_scaled[2:n] .= y_scaled[2:n] .* 1
        
        # Termos das extremidades (c_0 e c_n): Multiplicar por 2n (Inverso de 1/2n)
        y_scaled[[1, N]] .= y_scaled[[1, N]] .* (2 * 1)
        
        # 2. Aplicar DCT-I não normalizada
        y = FFTW.r2r(y_scaled, FFTW.REDFT00, 1)
        
        # 3. Inverter a ordem (de volta para a ordem física)
        y = reverse(y)

        # 4. Ajuste final de normalização (Fator 2 do FFTW)
        y = y ./ 2.0 
        
        return y

    else
        error("Direction must be 1 (nodal to spectral) or 2 (spectral to nodal)")
    end
end




# using FFTW
# using LinearAlgebra # Para sqrt

"""
    DCT_Chb_Gauss(y, direction)

Transformada de Chebyshev (DCT-II) para pontos de Chebyshev-Gauss (raízes).
Replica a normalização padrão do Matlab (DCT-II / IDCT-II).

# Argumentos
- `y`: Vetor coluna ou matriz de entrada (nodal ou espectral).
- `direction`: 1 (Nodal -> Espectral / DCT-II) ou 2 (Espectral -> Nodal / IDCT-II).

# Retorno
- `B`: O vetor ou matriz transformado(a).
"""
function DCT_Chb_Gauss(y, direction::Int)
    y = Float64.(y) # Garante o tipo Float64
    N = size(y, 1)

    # 1. Tratamento Trivial e Constantes
    if N == 1
        return y
    end
    
    # Fatores de Escala
    scale_direct = sqrt(2 / N)
    scale_inverse = sqrt(N / 2)

    # --- Transformada Direta (Nodal -> Espectral) ---
    # --- Transformada Direta (Nodal -> Espectral) ---
    if direction == 1
        yp = reverse(y)
        B = FFTW.r2r(yp, FFTW.REDFT10, 1)

        # 1. Normalização Total
        # A normalização correta é multiplicar por 2/N. 
        # O resultado do FFTW está com um fator de N. Vamos dividir por N.
        B = B ./ N 
        
        # 2. Multiplicar os termos internos por 2 (para dar 2/N)
        # O primeiro elemento (c0) é o único que não recebe este fator 2.
        #B[2:end, :] .= B[2:end, :] .* 2.0 
        
        # 3. Ajuste do Endpoint (c0): Multiplicar por 1/2 (já que o N/2 está corrigido)
        # B[1, :] precisa ser dividido por 2 no final
        B[1, :] .= B[1, :] ./ 2.0
        
        return B

# --- Transformada Inversa (Espectral -> Nodal) ---
   elseif direction == 2
        
        # 1. Aplicar a REVERSÃO da Normalização ANTES da FFTW
        y_scaled = copy(y)
        
        # A. Multiplicar o c0 (índice 1) por 2.0
        # Isso compensa a divisão por 2.0 que o c0 recebeu na Direta.
        y_scaled[1, :] .= y_scaled[1, :] .* 2.0
        
        # 2. Aplicar IDCT-II (FFTW.REDFT01). A SAÍDA JÁ VEM MULTIPLICADA POR N.
        y = FFTW.r2r(y_scaled, FFTW.REDFT01, 1)

        # 3. Inverter a ordem
        y = reverse(y)

        # 4. Ajuste final de normalização (Fator 2 do FFTW).
        # Este fator 2.0 final corrige a simetria de 2/N da Direta.
        y = y ./ 2.0 
        
        return y

    else
        error("Direction must be 1 (DCT-II) or 2 (IDCT-II).")
    end
end



# ====================================================================
# 2. Matriz de Diferenciação Modal de CHEBYSHEV I (D_Cheby)
#    Aplica-se aos coeficientes a[0] a a[N-1]
#    Dimensão: N x N
# ====================================================================
"""
    D_Cheby(N::Int)

Cria a matriz de diferenciação modal de Chebyshev I de dimensão N x N.
Os elementos D[i, j] transformam a[j] em a'[i].
Indices (i, j) vão de 1 a N, correspondendo aos índices modais 0 a N-1.
"""
function D_Cheby(n::Int)
    N = n+1;
    D = zeros(Float64, N, N)

    # i e j são os índices baseados em 1 (1 a N).
    # O índice modal (ordem do polinômio) é i_modal = i - 1 e j_modal = j - 1.
    for i in 1:N
        i_modal = i - 1
        for j in (i + 1):N
            j_modal = j - 1
            
            # Condição de paridade: i_modal + j_modal deve ser ímpar
            if isodd(i_modal + j_modal) 
                
                # Caso i_modal = 0 (Primeira linha da matriz - i=1)
                if i_modal == 0
                    # Condição adicional: j_modal deve ser ímpar.
                    if isodd(j_modal)
                        # Elemento é j_modal
                        D[i, j] = j_modal
                    end
                # Caso 0 < i_modal < j_modal (Linhas seguintes - i > 1)
                else
                    # Elemento é 2 * j_modal
                    D[i, j] = 2 * j_modal
                end
            end
        end
    end
    return D
end


"""
    SpectralQuadrature(n::Int, tipo::Int = 1)

Calcula os nós (z) e pesos (w) para integração espectral de f(x)
no intervalo [-1, 1] (ou seja, ∫ f(x) dx) usando `n` pontos.


"""
function ChebyQuadrature(a::Real, b::Real,n::Int, tipo::Int = 1)
    
    local z, w

    if n == 0
        error("Número de pontos 'n' deve ser >= 1.")
    elseif n == 1
        w = [2.0]
        if tipo == 1
            z = [0.0]
        elseif tipo == 2
            z = [0.0]
        elseif tipo == 3
            z = [-1.0] # Nó de Radau padrão
        end
        return (z, w)
    end


    if tipo == 1
        # --- 1. Quadratura de Fejér (Tipo 1) ---
        # Gauss-Cheby
        k = 1:n
        z = @. cos((2*k - 1) * pi / (2*n))
        theta = @. (2*k - 1) * pi / (2*n)
        
        w = zeros(n)
        for k_idx in 1:n
            theta_k = theta[k_idx]
            soma_j = 0.0
            for j in 1:floor(Int, (n-1)/2)
                soma_j += cos(2*j * theta_k) / (4*j^2 - 1)
            end
            w[k_idx] = (2/n) * (1 - 2 * soma_j)
        end

    elseif tipo == 2
        # --- 2. Quadratura de Clenshaw-Curtis (Trefethen O(N^2)) ---
        # Lobatto-Cheby
        N = n - 1 # Ordem do Polinômio
        
        k = 0:N
        theta = @. pi * k / N
        z = - cos.(theta)
        
        w = zeros(n)
        ii = 2:n-1
        v = ones(n - 2)
        
        theta_ii = theta[ii]
        
        if iseven(N)
            w[1] = 1.0 / (N^2 - 1)
            w[end] = w[1]
            for k_ = 1:(N/2 - 1)
                v .-= 2 .* cos.(2 * k_ .* theta_ii) ./ (4*k_^2 - 1)
            end
            v .-= cos.(N .* theta_ii) ./ (N^2 - 1)
        else
            w[1] = 1.0 / N^2
            w[end] = w[1]
            for k_ = 1:((N-1)/2)
                v .-= 2 .* cos.(2 * k_ .* theta_ii) ./ (4*k_^2 - 1)
            end
        end
        w[ii] .= 2 .* v ./ N
        
    elseif tipo == 3
        # --- 3. Quadratura de Fejér-Radau (Nó em -1) ---
        # (Implementação O(N^2) correta, baseada no seu código MATLAB)
        
        N_poly = n - 1  # Ordem do polinômio (n pontos)
        M = 2*N_poly + 1 # = 2*n - 1
        
        # 1. Nós
        k = 0:N_poly
        theta = (2 * pi / M) .* k
        z = -cos.(theta) # z[1] é -1
        
        # 2. Pesos
        w = zeros(n)
        
        # Coeficientes para j par: 2, 4, ..., N_poly
        j_even = 2:2:N_poly
        coeffs = 2.0 ./ (1.0 .- j_even.^2)
        
        # Loop sobre cada nó `k`
        for k_idx in 1:n
            theta_k = theta[k_idx]
            
            # Soma interna (só sobre j par)
            inner_sum = sum(coeffs .* cos.(j_even .* theta_k))
            
            # Termo constante + soma
            w[k_idx] = 1.0 + inner_sum
        end
        
        # 3. Normalização final
        w .*= (4.0 / M)
        
        # 4. Ajuste para o ponto x = -1
        w[1] /= 2.0
        
    else
        error("Tipo de quadratura inválido: $tipo. Use 1, 2, ou 3.")
    end
    
    z_scal = @. ((b - a) * z + b + a) / 2
    w_scal = @. w * (b - a) / 2

    return (z_scal, w_scal)
end

end
