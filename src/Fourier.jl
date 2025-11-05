module Fourier

using FFTW
using ToeplitzMatrices


export D_Fourier, Multi_Diff_Mat, Bary_Trig_Mat, I_Fourier, fourier_quad


# ====================================================================
# 3. Matriz de Diferenciação Modal de FOURIER REAL (D_Fourier)
#    Aplica-se aos coeficientes de Fourier Reais
#    Dimensão: (2m+1) x (2m+1) ou N x N (onde N é ímpar)
# ====================================================================
"""
    D_Fourier(M::Int)

Cria a matriz de diferenciação modal de Fourier Real.
M é a maior ordem de harmônicos (sin(Mx) e cos(Mx)).
O número total de coeficientes é N = 2M + 1.
A ordenação típica dos coeficientes é: a[0], a[1], b[1], a[2], b[2], ..., a[M], b[M].
Os índices modais i e j na fórmula correspondem a i_modal = i-1, j_modal = j-1 (0 a N-1).
"""
function D_Fourier(M::Int)
    N = 2 * M + 1
    D = zeros(Float64, N, N)

    # A fórmula fornecida usa a indexação baseada nos harmônicos k=1:M
    # A matriz é decomposta em blocos 2x2 para a diferenciação de (a_k, b_k)
    # A derivada de:
    # a_0 (k=0) é 0.
    # a_k*cos(kx) + b_k*sin(kx) é -k*a_k*sin(kx) + k*b_k*cos(kx)
    #               = k*b_k*cos(kx) - k*a_k*sin(kx)
    # Coeficientes da derivada: a'_k = k*b_k; b'_k = -k*a_k

    # A ordenação de índices (i, j) na matriz:
    # i=1 (a_0)
    # i=2 (a_1), i=3 (b_1)
    # i=4 (a_2), i=5 (b_2)
    # ...
    # O bloco 2k+1, 2k corresponde à transição entre b_k e a_k
    # A indexação de bloco é: 2k (linha/coluna para a_k) e 2k+1 (linha/coluna para b_k)

    for k in 1:M
        # a_k é o coeficiente na posição 2k
        # b_k é o coeficiente na posição 2k + 1

        # D[2k, 2k+1] = k  (a'_k = k * b_k)
        # O coeficiente a_k da derivada (linha 2k) recebe contribuição de b_k (coluna 2k+1)
        D[2k, 2k + 1] = k

        # D[2k+1, 2k] = -k (b'_k = -k * a_k)
        # O coeficiente b_k da derivada (linha 2k+1) recebe contribuição de a_k (coluna 2k)
        D[2k + 1, 2k] = -k
    end
    return D
end



"""
    Multi_Diff_Mat(N, m_max)

Calcula as matrizes de diferenciação espectral de Fourier para todas as ordens
de derivada de 1 até m_max.

# Argumentos
- `N`: Tamanho da matriz de diferenciação.
- `m_max`: Ordem MÁXIMA da derivada (inteiro).

# Retornos
- `x`: Vetor de pontos da grade.
- `DM_Collection`: Um vetor (lista) de matrizes, onde DM_Collection[k] é a k-ésima derivada.
"""
function Multi_Diff_Mat(N, m_max)
    
    # 1. Coleção para armazenar as matrizes de diferenciação
    DM_Collection = Matrix{Float64}[] # Inicializa um vetor de matrizes de Float64
    
    # 2. Loop de k=1 até m_max
    local x # O vetor x é o mesmo para todas as derivadas
    
    for k in 1:m_max
        # Chamamos a função fourdif para calcular APENAS a k-ésima derivada
        (x_k, D_k) = fourdif(N, k)
        
        # Armazenamos o resultado
        push!(DM_Collection, D_k)
        
        # Apenas precisamos do 'x' da primeira vez
        if k == 1
            x = x_k
        end
    end
    
    # 3. Empilhar as matrizes (Opcional, para replicar o MATLAB 3D)
    # Se você realmente precisa da matriz 3D, empilhamos:
    # DM_3D = cat(DM_Collection..., dims=3) 
    
    # Mas o mais idiomático é retornar a lista:
    return (x, DM_Collection)
end



# Differentiation of trig series
# Precisamos destes pacotes:
# - LinearAlgebra: para a matriz Toeplitz e funções trigonométricas (cot, csc)
# - FFTW: para as funções fft e ifft (necessário para m > 2)

"""
    fourdif(N, m)

Calcula a matriz de diferenciação espectral de Fourier de ordem `m` 
em uma grade com `N` pontos equi-espaçados em [0, 2π).

Tradução do código `fourdif.m` da `dmsuite` por Weideman e Reddy.

# Argumentos
- `N`: Tamanho da matriz de diferenciação (inteiro).
- `m`: Ordem da derivada (inteiro não-negativo).

# Retornos
- `x`: Vetor (coluna) de pontos da grade.
- `DM`: Matriz de diferenciação (N x N).
"""
function fourdif(N, m)
    
    # --- Configuração da Grade ---
    # `(0:N-1)` cria um intervalo (Range)
    # `collect()` o transforma em um vetor (coluna)
    x = 2π * collect(0:N-1) / N
    h = 2π / N
    
    # A unidade imaginária em Julia é `1im`
    zi = 1im 
    
    # Declaramos as variáveis que serão preenchidas pelo bloco if/else
    local col1, row1

    # --- Lógica principal para construir a primeira coluna e linha ---
    
    if m == 0
        # m = 0: Matriz identidade
        col1 = [1; zeros(N-1)]
        row1 = col1
        
    elseif m == 1
        # m = 1: Primeira derivada
        kk = collect(1:N-1)
        n1 = floor(Int, (N-1)/2)
        n2 = ceil(Int, (N-1)/2)
        
        if iseven(N) # `iseven(N)` é o mesmo que `rem(N,2)==0`
            # O ponto `.` em `cot.` aplica a função elemento a elemento
            topc = cot.((1:n2) .* (h/2))
            # `reverse` é o equivalente Julia de `flipud` para vetores
            col1 = [0; 0.5 .* ((-1) .^ kk) .* [topc; -reverse(topc[1:n1])]]
        else
            topc = csc.((1:n2) .* (h/2))
            col1 = [0; 0.5 .* ((-1) .^ kk) .* [topc; reverse(topc[1:n1])]]
        end
        row1 = -col1
        
    elseif m == 2
        # m = 2: Segunda derivada
        kk = collect(1:N-1)
        n1 = floor(Int, (N-1)/2)
        n2 = ceil(Int, (N-1)/2)
        
        if iseven(N)
            topc = csc.((1:n2) .* (h/2)) .^ 2
            col1 = [(-π^2 / (3*h^2) - 1/6); -0.5 .* ((-1) .^ kk) .* [topc; reverse(topc[1:n1])]]
        else
            topc = csc.((1:n2) .* (h/2)) .* cot.((1:n2) .* (h/2))
            col1 = [(-π^2 / (3*h^2) + 1/12); -0.5 .* ((-1) .^ kk) .* [topc; -reverse(topc[1:n1])]]
        end
        row1 = col1
        
    else
        # m > 2: Usa FFT
        N1 = floor(Int, (N-1)/2)
        
        # Esta linha lida com a frequência de Nyquist (para N par)
        # `ones(0)` em Julia é `[]`, `ones(1)` é `[1.0]`. Funciona.
        N2 = (-N/2) * ((m+1) % 2) .* ones((N+1) % 2)
        
        # Criamos o vetor de números de onda (como coluna)
        mwave = 1im .* [0:N1; N2; -N1:-1]
        
        # `fft` e `ifft` em Julia operam em colunas por padrão
        v_in = [1; zeros(N-1)]
        col1_complex = ifft((mwave .^ m) .* fft(v_in))
        
        # `real.` aplica a função `real` em cada elemento
        col1 = real.(col1_complex)
        
        if iseven(m)
            row1 = col1
        else
            col1 = [0; col1[2:N]]
            row1 = -col1
        end
    end

    # --- Construção da Matriz de Toeplitz ---
    # `Toeplitz` (com 'T' maiúsculo) cria um tipo de matriz especial "lazy"
    # `Matrix(...)` a converte para uma matriz densa padrão, como no MATLAB
    DM = Matrix(Toeplitz(col1, row1))
    
    return (x, DM)
end


 

"""
    Bary_Trig_Mat(fk, x)

Calcula o interpolante trigonométrico e a Matriz de Interpolação P.
Usa a fórmula baricêntrica trigonométrica, tratando os termos NaN/Inf
(o "truque de Berrut") para retornar uma Matriz P limpa.

VERSÃO CORRIGIDA (Denominador correto).
"""
function Bary_Trig_Mat(fk, x)
    fk_vec = vec(fk)
    x_vec = vec(x)
    
    N = length(fk_vec) # Número de pontos da grade
    M = length(x_vec)  # Número de pontos de avaliação
    
    # 1. Pontos da grade original
    xk = (2π / N) .* collect(0:N-1)
    
    # 2. Pesos baricêntricos trigonométricos (Vetor 1xN)
    w = ((-1) .^ collect(0:N-1))' # Transposto para (1xN)
    
    # 3. Termos angulares (MxN)
    delta_x_half = (x_vec ./ 2) .- (xk ./ 2)'
    
    # 4. Cálculo dos termos D_raw (sem o truque do eps)
    local D_raw
    if iseven(N)
        D_raw = 1 ./ tan.(delta_x_half)
    else
        D_raw = 1 ./ sin.(delta_x_half)
    end
    
    # 5. Cálculo dos termos C_ik = D_ik * w_k
    # (MxN) .* (1xN) -> (MxN)
    C = D_raw .* w 

    # 6. Cálculo do Denominador: Den[i] = sum_j(C_ij)
    # sum(C, dims=2) é (Mx1)
    Den = sum(C, dims=2)

    # 7. Cálculo de P (Com NaNs/Infs)
    # P[i, k] = C[i, k] / Den[i]
    P = C ./ Den 
    
    # --- O "TRUQUE DE BERRUT": Corrigir os NaNs ---
    
    # 8. Localizar os pontos de coincidência
    coincidencia_mat = isapprox.(delta_x_half, 0.0, atol=1e-15) .| isapprox.(abs.(delta_x_half), π, atol=1e-15)

    # 9. Corrigir a Matriz P nas linhas onde há coincidência
    for i in 1:M
        if any(coincidencia_mat[i, :])
            j_coincidente = findfirst(coincidencia_mat[i, :])
            P[i, :] .= 0.0 
            P[i, j_coincidente] = 1.0
        end
    end

    # 10. Cálculo final
    t = P * fk_vec
    
    return t, P
end


"""
    I_Fourier(M::Int)

Cria a matriz de integração modal de Fourier Real (N x N, N=2M+1).
Assume que a função a ser integrada tem média zero (a[0] = 0).
A constante de integração (coeficiente b[0]) é definida como zero.

Integração:
∫ (a_k*cos(kx) + b_k*sin(kx)) dx 
   = (a_k/k)*sin(kx) - (b_k/k)*cos(kx)
"""
function I_Fourier(M::Int)
    N = 2 * M + 1
    I = zeros(Float64, N, N)

    # Coeficientes da integral: a'_k = -b_k / k; b'_k = a_k / k
    # O coeficiente a[0] (índice 1) é ignorado (assume-se 0).
    # O coeficiente b[0] (índice 1) é 0 (constante de integração).

    for k in 1:M
        idx_ak = 2k
        idx_bk = 2k + 1
        I[1,2*k+1]=1/k;

        # I[idx_ak, idx_bk] = -1.0/k  (a'_k = -b_k / k)
        # Linha a'_k (idx_ak) recebe contribuição da coluna b_k (idx_bk)
        I[idx_ak, idx_bk] = -1.0 / k

        # I[idx_bk, idx_ak] = 1.0/k (b'_k = a_k / k)
        # Linha b'_k (idx_bk) recebe contribuição da coluna a_k (idx_ak)
        I[idx_bk, idx_ak] = 1.0 / k
    end
    return I
end

"""
Fourier quadrature
Calcula a integral definida de função periódica in [0,2*pi]
"""

"""
    fourier_quad(f::Function, N::Int; interval::Tuple{Float64, Float64}=(0.0, 2.0*pi))

Calcula a integral definida de uma função periódica `f` usando a 
quadratura de Fourier (Regra Trapezoidal Periódica) com N pontos.

Este método tem precisão espectral para funções periódicas suaves.

# Argumentos
- `f`: A função a ser integrada (ex: `x -> sin(x)`).
- `N`: O número de pontos de quadratura. Para consistência com seus coeficientes
       (a0, a1, b1, ..., aM, bM), este N deve ser N = 2*M + 1.
- `interval`: Tupla opcional `(a, b)` especificando o domínio de integração.
              O padrão é `(0.0, 2.0*pi)`. A função `f` deve ser
              periódica neste intervalo, i.e., `f(a) ≈ f(b)`.

# Retorna
- `Float64`: A estimativa da integral definida de `f` em `[a, b]`.
"""
function fourier_quad(f::Function, N::Int; interval::Tuple{Real, Real}=(0.0, 2.0*pi))
    if N <= 0
        error("N (ordem da expansão) deve ser um inteiro positivo.")
    end
    
    a, b = interval
    L = b - a # Comprimento do intervalo
    
    # O passo (delta x)
    dx = L / N
    
    # Cria N nós igualmente espaçados de 'a' até (b - dx).
    # O nó 'b' não é incluído, pois f(b) = f(a) pela periodicidade.
    # Usamos range() para criar os nós de forma eficiente.
    nodes = range(a, step=dx, length=N)
    
    # Avalia a função em todos os nós usando broadcasting (o 'f.(nodes)')
    f_values = f.(nodes)
    
    # A integral é a soma dos valores * o passo
    integral = sum(f_values) * dx
    
    return integral
end

end
