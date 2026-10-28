"""
sistemas_lineares.py
====================
Métodos DIRETOS para sistemas lineares [A]{x} = {b}.

    Eliminação:     gauss_ingenua, gauss_pivoteamento, gauss_jordan
    Determinante:   determinante, cramer, cramer_2x2
    Decomposição:   decomposicao_lu (Doolittle), decomposicao_crout,
                    resolver_lu, inversa_lu
    Especiais:      thomas (tridiagonal), cholesky + resolver_cholesky (simétrica)
    Auxiliares:     substituicao_progressiva, substituicao_regressiva, residuo,
                    novo_contador_flops, mostrar_sistema, vetores <-> tridiagonal

CONVENÇÃO DE ÍNDICES: o Python começa em 0. Nos comentários, "linha k+1" é a
linha k da apostila (que começa em 1). Nas mensagens impressas usamos sempre a
numeração da apostila (1..n) para facilitar a comparação com o caderno.

CONTAGEM DE FLOPs (opcional): passe um dicionário criado por
novo_contador_flops() no argumento 'flops'. Cada função soma ali as
multiplicações/divisões ("mult_div") e somas/subtrações ("soma_sub") que fez.
Atribuições, comparações, trocas de linha e índices de laço NÃO contam.

PRECISÃO REDUZIDA (opcional): passe 'algarismos=4' para arredondar o resultado
de CADA operação a 4 algarismos significativos (simula uma calculadora de
4 dígitos). Com algarismos=None (padrão) nada é arredondado.

IMPORTANTE: todas as funções trabalham numa CÓPIA em float de A e b. Assim a
matriz original nunca é alterada e não ocorre o erro de usar matriz de inteiros
(np.array([[1, 2], ...]) sem ponto decimal vira int e trunca as contas!).
"""

import math

import numpy as np


# ---------------------------------------------------------------------------
# Auxiliares
# ---------------------------------------------------------------------------
def novo_contador_flops():
    """Cria o dicionário de contagem de operações de ponto flutuante."""
    return {"mult_div": 0, "soma_sub": 0}


def _contar(flops, mult_div=0, soma_sub=0):
    """Soma operações ao contador, se ele foi fornecido."""
    if flops is not None:
        flops["mult_div"] += mult_div
        flops["soma_sub"] += soma_sub


def _arredondar(valor, algarismos):
    """Arredonda 'valor' para 'algarismos' significativos (None = não mexe)."""
    if algarismos is None:
        return valor
    return float(f"{valor:.{algarismos}g}")


def limite_pivo(A):
    """Limiar abaixo do qual um pivô é considerado numericamente nulo.

    Usa um critério RELATIVO à escala da matriz (e não um número fixo):
        limite = n * eps_maquina * max|a_ij|
    em que eps_maquina ~ 2.2e-16 para float64. Assim o teste funciona igual
    se a matriz estiver em N, kN ou mN.
    """
    A = np.asarray(A, dtype=float)
    return A.shape[0] * np.finfo(float).eps * np.max(np.abs(A))


def mostrar_sistema(A, b, saida=print, casas=6):
    """Imprime a matriz aumentada [A | b], uma linha por equação."""
    for i in range(len(b)):
        coef = "  ".join(f"{A[i, j]:>{casas + 7}.{casas}g}" for j in range(A.shape[1]))
        saida(f"   [ {coef} | {b[i]:>{casas + 7}.{casas}g} ]")


def residuo(A, x, b):
    """Vetor resíduo r = A x - b (deve ser ~0 se x estiver certo)."""
    return np.asarray(A, dtype=float) @ np.asarray(x, dtype=float) - np.asarray(b, dtype=float)


# ---------------------------------------------------------------------------
# Substituições (resolvem sistemas triangulares)
# ---------------------------------------------------------------------------
def substituicao_regressiva(U, y, flops=None, algarismos=None):
    """Resolve U x = y com U triangular SUPERIOR, de baixo para cima.

        x_n = y_n / u_nn
        x_i = (y_i - soma_{j=i+1..n} u_ij x_j) / u_ii,   i = n-1, ..., 1
    """
    n = len(y)
    x = np.zeros(n)
    x[n - 1] = _arredondar(y[n - 1] / U[n - 1, n - 1], algarismos)
    _contar(flops, mult_div=1)
    # o laço anda de baixo para cima porque x_i usa os x_j já calculados (j > i)
    for i in range(n - 2, -1, -1):
        soma = y[i]
        for j in range(i + 1, n):
            soma = _arredondar(soma - _arredondar(U[i, j] * x[j], algarismos), algarismos)
            _contar(flops, mult_div=1, soma_sub=1)
        x[i] = _arredondar(soma / U[i, i], algarismos)
        _contar(flops, mult_div=1)
    return x


def substituicao_progressiva(L, b, flops=None):
    """Resolve L y = b com L triangular INFERIOR, de cima para baixo.

        y_i = (b_i - soma_{j=1..i-1} l_ij y_j) / l_ii,   i = 1, ..., n
    """
    n = len(b)
    y = np.zeros(n)
    for i in range(n):
        soma = b[i]
        for j in range(i):
            soma -= L[i, j] * y[j]
            _contar(flops, mult_div=1, soma_sub=1)
        y[i] = soma / L[i, i]
        _contar(flops, mult_div=1)
    return y


# ---------------------------------------------------------------------------
# Eliminação gaussiana
# ---------------------------------------------------------------------------
def _eliminacao_progressiva(A, b, pivotear, flops, mostrar_etapas, saida, algarismos, trocas):
    """Transforma [A|b] num sistema triangular superior equivalente.

    Índices (iguais aos da apostila):
        k -> coluna do pivô e linha pivô (fica fixo enquanto zera a coluna k)
        i -> linha abaixo do pivô onde o elemento a_ik será zerado
        j -> colunas da linha i que precisam ser atualizadas
    """
    n = len(b)
    limite = limite_pivo(A)
    for k in range(n - 1):
        # ---- pivoteamento parcial: traz para a linha k o maior |a_pk| -------
        if pivotear:
            p = k
            for i in range(k + 1, n):
                if abs(A[i, k]) > abs(A[p, k]):
                    p = i
            if p != k:
                A[[k, p], :] = A[[p, k], :]     # troca as linhas k e p de A
                b[k], b[p] = b[p], b[k]         # e os termos de b
                if trocas is not None:
                    trocas.append((k + 1, p + 1))
        # ---- pivô nulo: a eliminação não pode continuar ----------------------
        if abs(A[k, k]) <= limite:
            raise ZeroDivisionError(f"Pivô nulo (ou ~0) na etapa k = {k + 1}: "
                                    f"a_kk = {A[k, k]:.3e}. "
                                    + ("Matriz singular." if pivotear else
                                       "Use pivoteamento parcial."))
        # ---- zera a coluna k abaixo da diagonal ------------------------------
        for i in range(k + 1, n):
            fator = _arredondar(A[i, k] / A[k, k], algarismos)   # multiplicador m_ik
            _contar(flops, mult_div=1)
            for j in range(k + 1, n):
                A[i, j] = _arredondar(A[i, j] - _arredondar(fator * A[k, j], algarismos),
                                      algarismos)
                _contar(flops, mult_div=1, soma_sub=1)
            b[i] = _arredondar(b[i] - _arredondar(fator * b[k], algarismos), algarismos)
            _contar(flops, mult_div=1, soma_sub=1)
            A[i, k] = 0.0          # zero exato (não precisa calcular)
        if mostrar_etapas:
            info = f" (linha pivô escolhida p = {p + 1})" if pivotear else ""
            saida(f"  Etapa k = {k + 1}{info}, pivô a_kk = {A[k, k]:.6g}:")
            mostrar_sistema(A, b, saida)
    if abs(A[n - 1, n - 1]) <= limite:
        raise ZeroDivisionError("Último pivô a_nn ~ 0: matriz singular.")


def gauss_ingenua(A, b, flops=None, mostrar_etapas=False, saida=print, algarismos=None):
    """Eliminação gaussiana ingênua (SEM troca de linhas) + substituição regressiva.

    Para se encontrar pivô nulo ou muito pequeno (não tenta corrigir).
    Custo: ~ 2n^3/3 flops na eliminação + n^2 na substituição.
    Retorna o vetor solução x.
    """
    A, b = np.array(A, dtype=float), np.array(b, dtype=float)
    _eliminacao_progressiva(A, b, False, flops, mostrar_etapas, saida, algarismos, None)
    return substituicao_regressiva(A, b, flops, algarismos)


def gauss_pivoteamento(A, b, flops=None, mostrar_etapas=False, saida=print,
                       algarismos=None, trocas=None):
    """Eliminação gaussiana com PIVOTEAMENTO PARCIAL.

    Em cada etapa k, procura nas linhas k..n o maior |a_ik| e troca essa linha
    com a linha k. Isso evita dividir por pivôs pequenos (que amplificam erros
    de arredondamento). Não muda o sistema físico, só a ordem das equações.

    trocas : lista opcional; recebe tuplas (k, p) de cada troca feita.
    Retorna o vetor solução x.
    """
    A, b = np.array(A, dtype=float), np.array(b, dtype=float)
    _eliminacao_progressiva(A, b, True, flops, mostrar_etapas, saida, algarismos, trocas)
    return substituicao_regressiva(A, b, flops, algarismos)


def gauss_jordan(A, b):
    """Gauss-Jordan com pivoteamento parcial: zera acima E abaixo do pivô e
    normaliza a diagonal para 1. No fim, b já é a solução (sem substituição).
    Custo ~ n^3 (cerca de 50% a mais que Gauss)."""
    A, b = np.array(A, dtype=float), np.array(b, dtype=float)
    n = len(b)
    limite = limite_pivo(A)
    for k in range(n):
        p = k + int(np.argmax(np.abs(A[k:, k])))      # linha com maior |a_ik|
        if abs(A[p, k]) <= limite:
            raise ZeroDivisionError("Gauss-Jordan: matriz singular.")
        A[[k, p], :] = A[[p, k], :]
        b[k], b[p] = b[p], b[k]
        # normaliza a linha pivô: a_kk passa a valer 1
        pivo = A[k, k]
        A[k, :] = A[k, :] / pivo
        b[k] = b[k] / pivo
        # zera a coluna k em TODAS as outras linhas (acima e abaixo)
        for i in range(n):
            if i != k:
                fator = A[i, k]
                A[i, :] = A[i, :] - fator * A[k, :]
                b[i] = b[i] - fator * b[k]
    return b


# ---------------------------------------------------------------------------
# Determinante e regra de Cramer
# ---------------------------------------------------------------------------
def determinante(A):
    """det(A) = (-1)^(no de trocas) * produto dos pivôs da eliminação."""
    U = np.array(A, dtype=float)
    n = U.shape[0]
    sinal = 1.0
    for k in range(n - 1):
        p = k + int(np.argmax(np.abs(U[k:, k])))
        if U[p, k] == 0.0:
            return 0.0                    # coluna toda nula: matriz singular
        if p != k:
            U[[k, p], :] = U[[p, k], :]
            sinal = -sinal                # cada troca de linha inverte o sinal
        for i in range(k + 1, n):
            fator = U[i, k] / U[k, k]
            U[i, k:] = U[i, k:] - fator * U[k, k:]
    produto = 1.0
    for k in range(n):
        produto *= U[k, k]
    return sinal * produto


def cramer(A, b):
    """Regra de Cramer: x_i = det(A_i)/det(A), onde A_i é A com a coluna i
    trocada por b. Só vale a pena para sistemas pequenos (2x2, 3x3)."""
    A, b = np.array(A, dtype=float), np.array(b, dtype=float)
    det_A = determinante(A)
    if det_A == 0.0:
        raise ZeroDivisionError("Cramer: det(A) = 0, sistema sem solução única.")
    x = np.zeros(len(b))
    for i in range(len(b)):
        A_i = A.copy()
        A_i[:, i] = b
        x[i] = determinante(A_i) / det_A
    return x


def cramer_2x2(a, b, c, d, r1, r2):
    """Resolve  a*x1 + b*x2 = r1  e  c*x1 + d*x2 = r2  pela regra de Cramer.

        D = a d - b c,   x1 = (r1 d - b r2)/D,   x2 = (a r2 - r1 c)/D

    Funciona com números ou com arrays NumPy (usado no fractal de Bairstow).
    """
    D = a * d - b * c
    return (r1 * d - b * r2) / D, (a * r2 - r1 * c) / D


# ---------------------------------------------------------------------------
# Decomposição LU
# ---------------------------------------------------------------------------
def decomposicao_lu(A, flops=None):
    """Decomposição LU de Doolittle (L com 1 na diagonal), sem pivoteamento.

    É a própria eliminação gaussiana aplicada só em A: os multiplicadores
    f_ik = a_ik/a_kk, que a eliminação "jogaria fora", são guardados em L.
    Retorna L, U com A = L U.
    """
    U = np.array(A, dtype=float)
    n = U.shape[0]
    L = np.eye(n)
    limite = limite_pivo(U)
    for k in range(n - 1):
        if abs(U[k, k]) <= limite:
            raise ZeroDivisionError(f"LU: pivô nulo na etapa k = {k + 1}.")
        for i in range(k + 1, n):
            fator = U[i, k] / U[k, k]
            L[i, k] = fator                 # guarda o multiplicador em L
            _contar(flops, mult_div=1)
            for j in range(k + 1, n):
                U[i, j] = U[i, j] - fator * U[k, j]
                _contar(flops, mult_div=1, soma_sub=1)
            U[i, k] = 0.0
    # o último pivô não entra no laço, mas também precisa ser não nulo
    if abs(U[n - 1, n - 1]) <= limite:
        raise ZeroDivisionError(f"LU: último pivô u_nn = {U[n - 1, n - 1]:.3e} ~ 0 "
                                "(matriz singular).")
    return L, U


def decomposicao_crout(A):
    """Decomposição de Crout: L geral e U com 1 na diagonal (A = L U).

    Para cada coluna j:
        l_ij = a_ij - soma_{k<j} l_ik u_kj                 (i = j..n)
        u_ji = (a_ji - soma_{k<j} l_jk u_ki) / l_jj         (i = j+1..n)
    """
    A = np.array(A, dtype=float)
    n = A.shape[0]
    L = np.zeros((n, n))
    U = np.eye(n)
    for j in range(n):
        for i in range(j, n):
            L[i, j] = A[i, j] - sum(L[i, k] * U[k, j] for k in range(j))
        if L[j, j] == 0.0:
            raise ZeroDivisionError(f"Crout: l_jj = 0 na coluna {j + 1}.")
        for i in range(j + 1, n):
            U[j, i] = (A[j, i] - sum(L[j, k] * U[k, i] for k in range(j))) / L[j, j]
    return L, U


def resolver_lu(L, U, b, flops=None):
    """Resolve A x = b já tendo A = L U:
        1) L d = b  por substituição progressiva
        2) U x = d  por substituição regressiva
    """
    d = substituicao_progressiva(L, np.asarray(b, dtype=float), flops)
    return substituicao_regressiva(U, d, flops)


def inversa_lu(A, flops=None):
    """Inversa de A coluna por coluna, com UMA única decomposição LU.

    A coluna j de A^-1 é a solução de A c_j = e_j (coluna j da identidade).
    L e U são reaproveitadas para todas as n colunas: custo ~ n^3 no total.
    """
    L, U = decomposicao_lu(A, flops)
    n = L.shape[0]
    inversa = np.zeros((n, n))
    for coluna in range(n):
        e = np.zeros(n)
        e[coluna] = 1.0                      # vetor unitário e_j
        inversa[:, coluna] = resolver_lu(L, U, e, flops)
    return inversa


# ---------------------------------------------------------------------------
# Matrizes especiais
# ---------------------------------------------------------------------------
def tridiagonal_para_vetores(A):
    """Extrai de uma matriz tridiagonal os vetores e (sub), f (diag), g (super).

    Convenção da apostila (com índice do Python):
        e[k] = a_{k,k-1}  (e[0] não existe -> 0)
        f[k] = a_{k,k}
        g[k] = a_{k,k+1}  (g[n-1] não existe -> 0)
    """
    A = np.asarray(A, dtype=float)
    n = A.shape[0]
    e, f, g = np.zeros(n), np.zeros(n), np.zeros(n)
    for k in range(n):
        f[k] = A[k, k]
        if k > 0:
            e[k] = A[k, k - 1]
        if k < n - 1:
            g[k] = A[k, k + 1]
    return e, f, g


def vetores_para_tridiagonal(e, f, g):
    """Monta a matriz cheia a partir de e, f, g (só para conferência)."""
    n = len(f)
    A = np.zeros((n, n))
    for k in range(n):
        A[k, k] = f[k]
        if k > 0:
            A[k, k - 1] = e[k]
        if k < n - 1:
            A[k, k + 1] = g[k]
    return A


def thomas(e, f, g, r, flops=None, historico=None):
    """Algoritmo de Thomas (TDMA) para sistemas TRIDIAGONAIS.

    É a decomposição LU especializada: guarda só 3n-2 números e faz ~8n
    operações (contra ~2n^3/3 da eliminação cheia). Três laços simples:

        1) Decomposição:  e_k = e_k / f_{k-1};   f_k = f_k - e_k g_{k-1}
        2) Progressiva:   r_k = r_k - e_k r_{k-1}
        3) Regressiva:    x_n = r_n / f_n;   x_k = (r_k - g_k x_{k+1}) / f_k

    Parâmetros: e (subdiagonal), f (diagonal), g (superdiagonal), r (lado
    direito), todos com tamanho n. Os vetores originais não são alterados.
    historico: lista opcional; recebe uma linha por passo de cada laço
    (é a "tabela feita à mão" da APC 3, com índices k da apostila, 1..n).
    """
    e, f, g, r = (np.array(v, dtype=float) for v in (e, f, g, r))
    n = len(f)
    # 1) decomposição (L guarda e_k, U guarda f_k e g_k)
    for k in range(1, n):
        e[k] = e[k] / f[k - 1]
        f[k] = f[k] - e[k] * g[k - 1]
        _contar(flops, mult_div=2, soma_sub=1)
        if historico is not None:
            historico.append({"laço": "1-decomposição", "k": k + 1, "e_k": e[k], "f_k": f[k]})
    # 2) substituição progressiva
    for k in range(1, n):
        r[k] = r[k] - e[k] * r[k - 1]
        _contar(flops, mult_div=1, soma_sub=1)
        if historico is not None:
            historico.append({"laço": "2-progressiva", "k": k + 1, "r_k": r[k]})
    # 3) substituição regressiva
    x = np.zeros(n)
    x[n - 1] = r[n - 1] / f[n - 1]
    _contar(flops, mult_div=1)
    if historico is not None:
        historico.append({"laço": "3-regressiva", "k": n, "x_k": x[n - 1]})
    for k in range(n - 2, -1, -1):
        x[k] = (r[k] - g[k] * x[k + 1]) / f[k]
        _contar(flops, mult_div=2, soma_sub=1)
        if historico is not None:
            historico.append({"laço": "3-regressiva", "k": k + 1, "x_k": x[k]})
    return x


def cholesky(A, flops=None):
    """Decomposição de Cholesky A = L L^T (A simétrica e positiva definida).

    Relações de recorrência (linha k, colunas i < k):
        l_ki = (a_ki - soma_{j<i} l_ij l_kj) / l_ii
        l_kk = sqrt(a_kk - soma_{j<k} l_kj^2)
    Custo ~ n^3/3 (metade da LU), pois só calcula L.
    """
    A = np.array(A, dtype=float)
    n = A.shape[0]
    if not np.allclose(A, A.T):
        raise ValueError("Cholesky: a matriz não é simétrica.")
    L = np.zeros((n, n))
    for k in range(n):
        # elementos fora da diagonal da linha k
        for i in range(k):
            soma = 0.0
            for j in range(i):
                soma += L[i, j] * L[k, j]
                _contar(flops, mult_div=1, soma_sub=1)
            L[k, i] = (A[k, i] - soma) / L[i, i]
            _contar(flops, mult_div=1, soma_sub=1)
        # elemento da diagonal
        soma = 0.0
        for j in range(k):
            soma += L[k, j] ** 2
            _contar(flops, mult_div=1, soma_sub=1)
        valor = A[k, k] - soma
        _contar(flops, soma_sub=1)
        if valor <= 0.0:
            raise ValueError("Cholesky: a matriz não é positiva definida "
                             f"(raiz de número <= 0 na linha {k + 1}).")
        L[k, k] = math.sqrt(valor)
    return L


def resolver_cholesky(L, b, flops=None):
    """Resolve A x = b com A = L L^T:  L y = b (progressiva) e L^T x = y (regressiva)."""
    y = substituicao_progressiva(L, np.asarray(b, dtype=float), flops)
    return substituicao_regressiva(L.T, y, flops)
