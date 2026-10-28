"""
sistemas_nao_lineares.py
========================
Newton-Raphson para sistemas de equações NÃO lineares F(x) = 0.

Ideia: em cada iteração, troca o sistema não linear pela sua aproximação
linear (série de Taylor de 1a ordem) em torno do chute atual x^(k):

    F(x^(k) + dx) ~ F(x^(k)) + J(x^(k)) dx = 0
    =>  J(x^(k)) dx = -F(x^(k))           (sistema LINEAR em dx)
    =>  x^(k+1) = x^(k) + dx

em que J é a matriz Jacobiana, J_ij = dF_i / dx_j.

O sistema linear de cada iteração é resolvido por:
    - regra de Cramer, se n = 2 (como pede a atividade das tubulações);
    - eliminação gaussiana com pivoteamento parcial, se n > 2.
"""

import numpy as np

from utils.sistemas_lineares import cramer_2x2, gauss_pivoteamento


def jacobiana_numerica(F, x, h_relativo=1e-7):
    """Jacobiana por diferenças progressivas, coluna por coluna:

        J[:, j] ~ [F(x + h e_j) - F(x)] / h
    """
    x = np.asarray(x, dtype=float)
    F0 = np.asarray(F(x), dtype=float)
    n = len(x)
    J = np.zeros((len(F0), n))
    for j in range(n):
        h = h_relativo * max(1.0, abs(x[j]))
        x_perturbado = x.copy()
        x_perturbado[j] += h
        J[:, j] = (np.asarray(F(x_perturbado), dtype=float) - F0) / h
    return J


def newton_sistemas(F, J, x0, tol=1e-8, max_iter=100, criterio="passo"):
    """Newton-Raphson multivariável.

    Parâmetros
    ----------
    F        : função que recebe o vetor x e devolve o vetor F(x)
    J        : função que devolve a Jacobiana J(x), ou None para usar a
               Jacobiana numérica (diferenças finitas)
    x0       : chute inicial (lista ou array)
    criterio : "passo"   -> max |dx_i| < tol
               "residuo" -> max |F_i(x)| < tol

    Retorna
    -------
    x, historico  (colunas: k, x1..xn, F1..Fn, dx1..dxn, erro)
    """
    x = np.array(x0, dtype=float)
    n = len(x)
    historico = []
    for k in range(1, max_iter + 1):
        Fx = np.asarray(F(x), dtype=float)                 # resíduos
        Jx = jacobiana_numerica(F, x) if J is None else np.asarray(J(x), dtype=float)
        # ---- resolve J dx = -F ------------------------------------------------
        if n == 2:
            dx = np.array(cramer_2x2(Jx[0, 0], Jx[0, 1], Jx[1, 0], Jx[1, 1], -Fx[0], -Fx[1]))
        else:
            dx = gauss_pivoteamento(Jx, -Fx)
        x_atual = x                       # valor usado nesta iteração
        x = x_atual + dx                  # atualização de Newton
        erro = np.max(np.abs(dx)) if criterio == "passo" else np.max(np.abs(F(x)))
        # ---- registra a iteração ---------------------------------------------
        linha = {"k": k}
        for i in range(n):
            linha[f"x{i + 1}"] = x_atual[i]
        for i in range(n):
            linha[f"F{i + 1}"] = Fx[i]
        for i in range(n):
            linha[f"dx{i + 1}"] = dx[i]
        linha["erro"] = erro
        historico.append(linha)
        if not np.all(np.isfinite(x)):
            raise ArithmeticError("Newton para sistemas divergiu; mude o chute inicial.")
        if erro < tol:
            return x, historico
    print(f"AVISO: Newton para sistemas NÃO convergiu em {max_iter} iterações.")
    return x, historico
