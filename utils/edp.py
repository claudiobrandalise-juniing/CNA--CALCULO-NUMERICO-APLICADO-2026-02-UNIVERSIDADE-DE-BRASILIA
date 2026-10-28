"""
edp.py
======
Diferenças finitas para a equação do calor.

1) REGIME PERMANENTE 2D (elíptica): Laplace / Poisson
       d2T/dx2 + d2T/dy2 + q'''/k = 0
   Malha uniforme (dx = dy), diferenças centrais em cada nó interno:
       T_ij = (T_(i-1)j + T_(i+1)j + T_i(j-1) + T_i(j+1) + q''' dx^2 / k) / 4
   Resolvido por Gauss-Seidel aplicado direto na malha (método de Liebmann),
   sem montar a matriz. Relaxação opcional (SOR).
       -> gauss_seidel_laplace_2d, residuo_laplace_2d

2) REGIME TRANSIENTE 1D (parabólica):
       dT/dt = alfa d2T/dx2 + q'''/(rho cp)
   Esquema com parâmetro beta (PPC 3):
       T^{p+1}_m - T^p_m = Fo [ beta (T^p_{m+1} + T^p_{m-1} - 2 T^p_m)
                         + (1 - beta)(T^{p+1}_{m+1} + T^{p+1}_{m-1} - 2 T^{p+1}_m) ] + A
       Fo = alfa dt / dx^2,   A = q''' dt / (rho cp) = q''' alfa dt / k
       beta = 0   -> implícito (incondicionalmente estável)
       beta = 1/2 -> Crank-Nicolson (2a ordem no tempo, estável)
       beta = 1   -> explícito (estável só se Fo <= 1/2)
   O sistema de cada passo é TRIDIAGONAL -> algoritmo de Thomas.
       -> conducao_transiente_1d

CONVENÇÃO DE ÍNDICES 2D: T[i, j] com i na direção x (horizontal) e j na
direção y (vertical); j = 0 é a borda de BAIXO. Para desenhar com imshow use
plt.imshow(T.T, origin="lower"), que põe x na horizontal e y para cima.
"""

import numpy as np

from utils.sistemas_lineares import thomas


# ---------------------------------------------------------------------------
# 1) Regime permanente 2D
# ---------------------------------------------------------------------------
def gauss_seidel_laplace_2d(T, tol=1e-6, max_iter=100000, relaxacao=1.0, fonte=None):
    """Gauss-Seidel (Liebmann) na malha 2D.

    Parâmetros
    ----------
    T         : array 2D (float) com as temperaturas iniciais; as BORDAS já
                devem conter as condições de contorno (Dirichlet) e não mudam.
    tol       : para quando a maior variação entre duas varreduras < tol
    relaxacao : lambda (1 = Gauss-Seidel puro, 1 < lambda < 2 = SOR)
    fonte     : None ou array 2D com o termo q''' dx^2 / k de cada nó

    Retorna
    -------
    T     : campo convergido (cópia; o array de entrada não é alterado)
    erros : lista com a maior variação de cada varredura (len = iterações)
    """
    T = np.array(T, dtype=float)              # cópia em float (evita truncar!)
    nx, ny = T.shape
    S = np.zeros_like(T) if fonte is None else np.asarray(fonte, dtype=float)
    erros = []
    for _ in range(max_iter):
        erro = 0.0
        # varre os nós internos linha a linha (j = vertical, i = horizontal)
        for j in range(1, ny - 1):
            for i in range(1, nx - 1):
                antigo = T[i, j]
                # média dos 4 vizinhos (usando os valores JÁ atualizados)
                novo = (T[i - 1, j] + T[i + 1, j] + T[i, j - 1] + T[i, j + 1] + S[i, j]) / 4.0
                T[i, j] = relaxacao * novo + (1.0 - relaxacao) * antigo
                erro = max(erro, abs(T[i, j] - antigo))
        erros.append(erro)
        if erro < tol:
            return T, erros
    print(f"AVISO: Gauss-Seidel 2D NÃO convergiu em {max_iter} varreduras.")
    return T, erros


def residuo_laplace_2d(T, fonte=None):
    """Maior |r_ij| com r_ij = 4 T_ij - (soma dos 4 vizinhos) - S_ij nos nós
    internos. Verifica se as equações discretizadas foram satisfeitas."""
    T = np.asarray(T, dtype=float)
    S = np.zeros_like(T) if fonte is None else np.asarray(fonte, dtype=float)
    r = (4.0 * T[1:-1, 1:-1] - T[:-2, 1:-1] - T[2:, 1:-1]
         - T[1:-1, :-2] - T[1:-1, 2:] - S[1:-1, 1:-1])
    return float(np.max(np.abs(r)))


# ---------------------------------------------------------------------------
# 2) Regime transiente 1D
# ---------------------------------------------------------------------------
def _linha_contorno(contorno, T_no, T_vizinho, Fo, beta, dx, k, A):
    """Coeficientes (diagonal, vizinho, lado direito) de um nó de contorno.

    Tipos aceitos (tuplas):
        ("adiabatico",) ou ("simetria",) -> fluxo nulo: nó fantasma = vizinho
        ("conveccao", h, T_inf)          -> balanço de energia em meia célula
        ("fluxo", q_entrando)            -> fluxo de calor q'' [W/m2] entrando
        ("temperatura", T_fixa)          -> Dirichlet
    Nos três primeiros o nó tem meia célula (dx/2), daí os fatores 2.
    """
    implicito = 1.0 - beta                       # peso do instante p+1
    tipo = contorno[0]
    if tipo == "temperatura":
        return 1.0, 0.0, contorno[1]
    # parte comum de adiabático, convecção e fluxo
    diagonal = 1.0 + 2.0 * Fo * implicito
    vizinho = -2.0 * Fo * implicito
    lado_direito = T_no + 2.0 * Fo * beta * (T_vizinho - T_no) + A
    if tipo in ("adiabatico", "simetria"):
        return diagonal, vizinho, lado_direito
    if tipo == "conveccao":
        h, T_inf = contorno[1], contorno[2]
        Bi_no = h * dx / k                     # Biot da MALHA (h dx / k), não h L / k!
        diagonal += 2.0 * Fo * implicito * Bi_no
        lado_direito += -2.0 * Fo * beta * Bi_no * T_no + 2.0 * Fo * Bi_no * T_inf
        return diagonal, vizinho, lado_direito
    if tipo == "fluxo":
        lado_direito += 2.0 * Fo * contorno[1] * dx / k
        return diagonal, vizinho, lado_direito
    raise ValueError(f"Contorno desconhecido: {contorno}")


def matriz_tridiagonal_1d(n_nos, Fo, beta, dx, k, contorno_esq, contorno_dir):
    """Vetores (e, f, g) da matriz do esquema beta (não dependem do tempo).

    Útil para conferir a montagem "no papel": com beta = 0 e 7 nós, reproduz a
    matriz [A] do PPC 3 (Fo aqui é o Fo_dt = alfa dt/dx^2 da apostila).
    """
    implicito = 1.0 - beta
    e = np.full(n_nos, -Fo * implicito)              # subdiagonal
    f = np.full(n_nos, 1.0 + 2.0 * Fo * implicito)   # diagonal
    g = np.full(n_nos, -Fo * implicito)              # superdiagonal
    e[0] = g[-1] = 0.0                               # posições que não existem
    # a 1a e a última linha vêm do balanço nos nós de contorno
    f[0], g[0], _ = _linha_contorno(contorno_esq, 0.0, 0.0, Fo, beta, dx, k, 0.0)
    f[-1], e[-1], _ = _linha_contorno(contorno_dir, 0.0, 0.0, Fo, beta, dx, k, 0.0)
    return e, f, g


def conducao_transiente_1d(L, n_nos, alfa, k, dt, n_passos, T_inicial,
                           contorno_esq=("adiabatico",), contorno_dir=("adiabatico",),
                           geracao=0.0, beta=0.0, salvar_a_cada=1):
    """Condução transiente 1D numa parede de espessura L.

    Parâmetros
    ----------
    L, n_nos       : comprimento [m] e número de nós (dx = L/(n_nos - 1))
    alfa, k        : difusividade [m2/s] e condutividade [W/(m K)]
    dt, n_passos   : passo de tempo [s] e número de passos
    T_inicial      : número (campo uniforme) ou array com n_nos valores
    contorno_esq/dir : tuplas descritas em _linha_contorno (nó 0 e nó n-1)
    geracao        : q''' [W/m3] (uniforme)
    beta           : 0 implícito, 0.5 Crank-Nicolson, 1 explícito
    salvar_a_cada  : guarda o campo a cada tantos passos

    Retorna
    -------
    x        : posições dos nós [m]
    t_salvos : instantes guardados [s]
    T_salvos : array [instante, nó] com as temperaturas
    """
    dx = L / (n_nos - 1)
    x = np.linspace(0.0, L, n_nos)
    Fo = alfa * dt / dx ** 2                      # Fourier da malha
    A = geracao * alfa * dt / k                   # termo de geração por passo
    # estabilidade: só preocupa quando o peso explícito domina (beta > 1/2)
    if beta > 0.5 and Fo > 1.0 / (2.0 * (2.0 * beta - 1.0)):
        print(f"AVISO: Fo = {Fo:.3g} viola a estabilidade do esquema com beta = {beta}.")

    # ---- matriz tridiagonal (constante no tempo) -----------------------------
    e, f, g = matriz_tridiagonal_1d(n_nos, Fo, beta, dx, k, contorno_esq, contorno_dir)

    T = np.full(n_nos, float(T_inicial)) if np.ndim(T_inicial) == 0 \
        else np.array(T_inicial, dtype=float)
    t_salvos, T_salvos = [0.0], [T.copy()]
    for p in range(1, n_passos + 1):
        # ---- lado direito: parte explícita (instante p) + geração ------------
        r = np.empty(n_nos)
        r[1:-1] = T[1:-1] + Fo * beta * (T[2:] + T[:-2] - 2.0 * T[1:-1]) + A
        # ---- lado direito dos nós de contorno (balanço de energia) ---------------
        r[0] = _linha_contorno(contorno_esq, T[0], T[1], Fo, beta, dx, k, A)[2]
        r[-1] = _linha_contorno(contorno_dir, T[-1], T[-2], Fo, beta, dx, k, A)[2]
        # ---- resolve o sistema tridiagonal do passo -------------------------
        T = thomas(e, f, g, r)
        # guarda a cada 'salvar_a_cada' passos e SEMPRE o último passo
        if p % salvar_a_cada == 0 or p == n_passos:
            t_salvos.append(p * dt)
            T_salvos.append(T.copy())
    return x, np.array(t_salvos), np.array(T_salvos)
