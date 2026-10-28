"""
sistemas_iterativos.py
======================
Métodos ITERATIVOS para [A]{x} = {b}: partem de um chute e refinam a solução.

    jacobi        -> cada x_i novo usa SÓ valores da iteração anterior
    gauss_seidel  -> cada x_i novo usa IMEDIATAMENTE os x_j já atualizados;
                     aceita relaxação lambda:
                         x_i <- lambda * x_i(GS) + (1 - lambda) * x_i(antigo)
                         lambda = 1       -> Gauss-Seidel puro
                         0 < lambda < 1   -> sub-relaxação (ajuda a convergir)
                         1 < lambda < 2   -> sobre-relaxação, SOR (acelera)

    dominancia_diagonal          -> testa |a_ii| > soma_{j!=i} |a_ij|
    reordenar_para_dominancia    -> tenta trocar a ordem das equações para que
                                    a matriz fique diagonalmente dominante

CONVERGÊNCIA: dominância diagonal estrita é condição SUFICIENTE (não
necessária). Matrizes simétricas positivas definidas também convergem com
Gauss-Seidel. É a mesma ideia da iteração de ponto fixo: x = g(x).

CRITÉRIO DE PARADA ('criterio'):
    "relativo" -> max_i |(x_i novo - x_i antigo)/x_i novo| * 100 < tol  (Chapra, %)
    "absoluto" -> max_i |x_i novo - x_i antigo| < tol
"""

import numpy as np


def dominancia_diagonal(A):
    """Retorna True se |a_ii| > soma_{j != i} |a_ij| em TODAS as linhas."""
    A = np.asarray(A, dtype=float)
    for i in range(A.shape[0]):
        fora_da_diagonal = np.sum(np.abs(A[i, :])) - abs(A[i, i])
        if abs(A[i, i]) <= fora_da_diagonal:
            return False
    return True


def reordenar_para_dominancia(A, b):
    """Procura, para cada coluna i, uma equação cujo maior coeficiente (em
    módulo) esteja nessa coluna e domine a linha. Retorna (A_nova, b_nova,
    ordem) ou levanta ValueError se não houver reordenação dominante.

    'ordem[i]' = índice (Python) da equação original que foi para a linha i.
    """
    A, b = np.asarray(A, dtype=float), np.asarray(b, dtype=float)
    n = A.shape[0]
    ordem = [None] * n
    for linha in range(n):
        coluna_dominante = int(np.argmax(np.abs(A[linha, :])))
        resto = np.sum(np.abs(A[linha, :])) - abs(A[linha, coluna_dominante])
        # a equação só serve se o maior coeficiente dominar a soma dos outros
        if abs(A[linha, coluna_dominante]) > resto and ordem[coluna_dominante] is None:
            ordem[coluna_dominante] = linha
    if None in ordem:
        raise ValueError("Não existe reordenação que torne a matriz diagonalmente dominante.")
    return A[ordem, :], b[ordem], ordem


def _erro_iterativo(x_novo, x_antigo, criterio):
    """Erro de uma componente segundo o critério escolhido."""
    if criterio == "absoluto" or x_novo == 0.0:
        return abs(x_novo - x_antigo)
    return abs((x_novo - x_antigo) / x_novo) * 100.0


def _iterar(A, b, x0, tol, max_iter, relaxacao, criterio, usar_atualizados, nome):
    """Laço comum de Jacobi (usar_atualizados=False) e Gauss-Seidel (True)."""
    A, b = np.array(A, dtype=float), np.array(b, dtype=float)
    n = len(b)
    # verifica a diagonal: x_i = (...)/a_ii exige a_ii != 0
    for i in range(n):
        if A[i, i] == 0.0:
            raise ZeroDivisionError(f"{nome}: a_ii = 0 na linha {i + 1}; reordene as equações.")
    x = np.zeros(n) if x0 is None else np.array(x0, dtype=float)
    historico = []
    for k in range(1, max_iter + 1):
        # Jacobi precisa de uma cópia congelada da iteração anterior
        x_anterior = x.copy()
        erro_max = 0.0
        for i in range(n):
            # vetor usado para os vizinhos: atualizado (GS) ou antigo (Jacobi)
            fonte = x if usar_atualizados else x_anterior
            soma = b[i]
            for j in range(n):
                if j != i:
                    soma -= A[i, j] * fonte[j]
            x_metodo = soma / A[i, i]
            # relaxação: média ponderada entre o valor novo e o antigo
            x[i] = relaxacao * x_metodo + (1.0 - relaxacao) * x_anterior[i]
            erro_max = max(erro_max, _erro_iterativo(x[i], x_anterior[i], criterio))
        linha = {"k": k}
        for i in range(n):
            linha[f"x{i + 1}"] = x[i]
        linha["erro"] = erro_max
        historico.append(linha)
        if erro_max < tol:
            return x, historico
    print(f"AVISO: {nome} NÃO convergiu em {max_iter} iterações "
          f"(dominância diagonal: {dominancia_diagonal(A)}).")
    return x, historico


def jacobi(A, b, x0=None, tol=1e-6, max_iter=1000, criterio="relativo"):
    """Iteração de Jacobi: x_i^(k+1) = (b_i - soma_{j!=i} a_ij x_j^(k)) / a_ii.
    Retorna x, historico (colunas: k, x1..xn, erro)."""
    return _iterar(A, b, x0, tol, max_iter, 1.0, criterio, False, "Jacobi")


def gauss_seidel(A, b, x0=None, tol=1e-6, max_iter=1000, relaxacao=1.0, criterio="relativo"):
    """Gauss-Seidel com relaxação opcional.

        x_i^(k+1) = (b_i - soma_{j<i} a_ij x_j^(k+1) - soma_{j>i} a_ij x_j^(k)) / a_ii
        x_i^(k+1) <- lambda x_i^(k+1) + (1 - lambda) x_i^(k)

    Retorna x, historico (colunas: k, x1..xn, erro).
    """
    return _iterar(A, b, x0, tol, max_iter, relaxacao, criterio, True, "Gauss-Seidel")
