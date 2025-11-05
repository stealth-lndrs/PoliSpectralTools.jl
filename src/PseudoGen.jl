module PseudoGen

using LinearAlgebra

export Bary_Interp, Generalized_Diff_Mat



"""
    Bary_Interp(xk, fk, xnew)

Calcula o interpolador polinomial dos dados (xk, fk) usando o método baricêntrico.

# Argumentos de entrada
- `xk`: Vetor de coordenadas x dos dados (assumidos distintos).
- `fk`: Vetor de coordenadas y (valores da função) dos dados.
- `xnew`: Vetor de valores x onde o polinômio interpolador deve ser avaliado.

# Retorno
- `fnew`: Vetor de valores interpolados.
- `P`: Matriz de Interpolação.
"""
function Bary_Interp(xk, fk, xnew)
    # 1. Garantir que os dados sejam vetores coluna (Julia usa colunas por padrão)
    xnew = vec(xnew)
    xk = vec(xk)
    fk = vec(fk)

    N = length(xk)
    M = length(xnew)

    # 2. Cálculo dos Pesos (w)
    # Em Julia, a função `prod` não é vetorizada da mesma forma que no Matlab,
    # e a construção da matriz de diferenças precisa de uma abordagem Julia-idiomática.
    
    # Criar a matriz de diferenças D, onde D[j, k] = xk[j] - xk[k]
    D = xk .- xk'

    # Onde j = k, D[j, j] = 0, mas para o peso w(k), o produto deve excluir este termo.
    # No Matlab, a linha D(L) = ones(N,1) resolve isso.
    # Em Julia, fazemos o seguinte para garantir o produto correto:
    for i in 1:N
        D[i, i] = 1.0  # Temporariamente, substitui a diagonal por 1 para não afetar o prod
    end

    # D vai ser uma matriz N x N
    # O produto é feito ao longo das linhas para obter o produto (xk[k] - xk[j]) para j != k
    w = 1.0 ./ prod(D, dims=2) # prod(D, dims=2) calcula o produto de cada linha
    w = vec(w) # Retorna w a um vetor

    # 3. Cálculo da Matriz de Interpolação (P)
    
    # deltaX[j, k] = xnew[j] - xk[k] (Matlab usa a transposta, aqui a ordem é mais direta)
    # O broadcasting com ponto (.) garante a operação entre o vetor coluna e o vetor linha transposto
    deltaX = xnew .- xk'

    # P[j, k] = w[k] / (xnew[j] - xk[k])
    P = w' ./ deltaX
    
    # Normalização: P[j, :] = P[j, :] / sum(P[j, :])
    P = P ./ sum(P, dims=2)
    
    # 4. Remoção de NaNs (Ocorre quando xnew[j] = xk[k])
    # Neste caso, a interpolação é trivialmente fnew = fk.
    # P[j, k] deve ser 1 e os outros P[j, l] devem ser 0, o que é garantido pela normalização.
    
    # isnan.(P) cria uma máscara booleana
    P[isnan.(P)] .= 1.0
    
    # 5. Cálculo do resultado
    fnew = P * fk
    
    return fnew, P
end


"""
    Generalized_Diff_Mat(xs)

Calcula a Matriz de Diferenciação (D) de primeira ordem para 
qualquer conjunto de pontos nodais 'xs' (assumidos distintos), 
usando a interpolação polinomial de Lagrange.

# Argumentos
- `xs`: Vetor de coordenadas x dos pontos nodais.

# Retorno
- `D`: A Matriz de Diferenciação (Matriz Quadrada de tamanho N x N).
"""
function Generalized_Diff_Mat(xs)
    # 1. Preparação dos Dados
    xs = Float64.(xs)
    N = length(xs)
    
    # 2. Matrizes de Diferença (sem modificação)
    # dx_zeros: Matriz onde dx_zeros[i, j] = xs[i] - xs[j] (com 0 na diagonal)
    dx_zeros = xs .- xs' 
    
    # 3. Cálculo dos Coeficientes de Lagrange (aj): PRODUTO SEM A DIAGONAL
    aj = ones(N) # Inicializa aj com 1.0

    for j in 1:N
        # Para cada ponto j, o produto é (x_j - x_i) para todos i != j
        # O prod é feito apenas nos termos fora da diagonal.
        
        # A forma mais limpa é: calcular o produto de (xs[j] - xs[k]) onde k != j
        # Usamos o dX com 1s na diagonal temporariamente:
        dX_temp = copy(dx_zeros) 
        dX_temp[j, j] = 1.0 # Garante que o termo diagonal é 1 (para o produto)
        
        # O produto de cada coluna j (dims=1)
        aj[j] = prod(dX_temp, dims=1)[j]
    end
    
    # 4. Pré-cálculo da Matriz de Diferenciação D
    # Inicializamos D com zeros e calculamos APENAS os termos fora da diagonal.
    D = zeros(N, N)

    for i in 1:N
        for j in 1:N
            if i != j
                # D[i, j] = (aj[i] / aj[j]) * (1 / (xs[i] - xs[j]))
                D[i, j] = (aj[i] / aj[j]) / dx_zeros[i, j]
            end
        end
    end

    # 5. Correção da Diagonal (D[i, i] = - sum(D[i, j]) para j != i)
    # sum(D, dims=2) calcula a soma dos termos fora da diagonal.
    # D = D - Diagonal(...) garante que D[i, i] = D[i, i] - sum(linha_i)
    
    # D_base é a matriz com D[i, i] = 0.0
    # O passo de subtração calcula D[i, i] = 0.0 - (soma dos off-diagonais)
    D = D - Diagonal(vec(sum(D, dims=2)))

    return D
end

end
