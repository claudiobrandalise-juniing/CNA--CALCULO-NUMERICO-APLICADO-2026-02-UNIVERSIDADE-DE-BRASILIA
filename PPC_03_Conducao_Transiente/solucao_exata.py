"""
solucao_exata.py (PPC 03)
=========================
Solução exata (separação de variáveis) da condução transiente numa parede
plana SEM geração, com simetria em x* = 0 e convecção em x* = 1:

    theta(x*, Fo) = soma_n C_n cos(lambda_n x*) exp(-lambda_n^2 Fo)
    C_n = 4 sin(lambda_n) / (2 lambda_n + sin(2 lambda_n))
    lambda_n  raízes de  lambda tan(lambda) = Bi      (Bi = h L / k)

Os autovalores lambda_n são obtidos pelo MÉTODO DA BISSEÇÃO (utils/raizes.py).
Truque: em vez de lambda tan(lambda) - Bi, que tem assíntotas verticais, usa-se
a forma equivalente SEM singularidades

    g(lambda) = lambda sin(lambda) - Bi cos(lambda) = 0

e o fato de que a n-ésima raiz está sempre em [(n-1) pi, (n-1) pi + pi/2],
onde g troca de sinal. Assim a bisseção nunca falha.
"""

import math

import numpy as np

from utils.raizes import bissecao


def autovalores(Bi, n_termos):
    """Primeiras n_termos raízes de lambda tan(lambda) = Bi (Bi > 0)."""
    def g(lam):
        return lam * math.sin(lam) - Bi * math.cos(lam)
    raizes = []
    for n in range(1, n_termos + 1):
        inicio = (n - 1) * math.pi
        lam, _ = bissecao(g, inicio, inicio + math.pi / 2, tol=1e-13, max_iter=200, criterio="passo")
        raizes.append(lam)
    return np.array(raizes)


def theta_exato(x_adim, Fo, Bi, n_termos=200, lambdas=None):
    """theta(x*, Fo) pela série truncada em n_termos (x_adim pode ser array).

    Para Fo muito pequeno a série converge devagar; com Fo >= 0.01 e 200 termos
    o erro de truncamento é desprezível.
    """
    lam = autovalores(Bi, n_termos) if lambdas is None else lambdas
    C = 4.0 * np.sin(lam) / (2.0 * lam + np.sin(2.0 * lam))
    x = np.atleast_1d(np.asarray(x_adim, dtype=float))
    # soma dos termos: cada coluna é um termo n, cada linha uma posição x
    termos = C * np.cos(np.outer(x, lam)) * np.exp(-lam ** 2 * Fo)
    return termos.sum(axis=1)
