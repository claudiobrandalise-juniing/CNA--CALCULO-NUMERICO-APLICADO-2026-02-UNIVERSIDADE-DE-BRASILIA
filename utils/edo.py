"""
edo.py
======
Métodos de passo único para problemas de valor inicial (PVI):

    dy/dt = f(t, y),   y(t0) = y0

y pode ser um número (uma EDO) ou um vetor (sistema de EDOs). Uma EDO de
ordem 2 ou mais vira sistema de 1a ordem trocando variáveis
(ex.: y'' = g -> y1 = y, y2 = y', y1' = y2, y2' = g).

    euler -> y_{i+1} = y_i + f(t_i, y_i) h                          erro global O(h)
    heun  -> preditor (Euler) + corretor (média das inclinações)    O(h^2)
    rk4   -> Runge-Kutta clássico de 4a ordem                       O(h^4)

RK4 (APC 1):
    k1 = f(t_i, y_i)
    k2 = f(t_i + h/2, y_i + k1 h/2)
    k3 = f(t_i + h/2, y_i + k2 h/2)
    k4 = f(t_i + h, y_i + k3 h)
    y_{i+1} = y_i + (k1 + 2 k2 + 2 k3 + k4) h / 6

Todas retornam (t, y): t com n_passos + 1 instantes e y com os valores
correspondentes (1D se y0 for número; 2D [instante, variável] se for vetor).
"""

import numpy as np


def _passo_euler(f, t, y, h):
    """Euler explícito: usa só a inclinação no início do passo."""
    return y + f(t, y) * h


def _passo_heun(f, t, y, h):
    """Heun (sem iterar o corretor): média das inclinações do início e do fim."""
    inclinacao_inicio = f(t, y)
    y_preditor = y + inclinacao_inicio * h
    inclinacao_fim = f(t + h, y_preditor)
    return y + (inclinacao_inicio + inclinacao_fim) / 2.0 * h


def _passo_rk4(f, t, y, h):
    """Runge-Kutta de 4a ordem: média ponderada de 4 inclinações."""
    k1 = f(t, y)
    k2 = f(t + h / 2.0, y + k1 * h / 2.0)
    k3 = f(t + h / 2.0, y + k2 * h / 2.0)
    k4 = f(t + h, y + k3 * h)
    return y + (k1 + 2.0 * k2 + 2.0 * k3 + k4) * h / 6.0


def _integrar(passo, f, t0, y0, h, n_passos):
    """Laço de avanço no tempo comum aos três métodos."""
    escalar = np.ndim(y0) == 0                  # True se for uma única EDO
    y = np.atleast_1d(np.array(y0, dtype=float))
    # a função do usuário sempre recebe e devolve arrays do mesmo tamanho
    def f_vetorial(t, y_vetor):
        return np.atleast_1d(np.asarray(f(t, y_vetor[0] if escalar else y_vetor), dtype=float))
    t = t0 + h * np.arange(n_passos + 1)
    Y = np.zeros((n_passos + 1, y.size))
    Y[0] = y
    for i in range(n_passos):
        Y[i + 1] = passo(f_vetorial, t[i], Y[i], h)
    return t, (Y[:, 0] if escalar else Y)


def euler(f, t0, y0, h, n_passos):
    """Método de Euler explícito."""
    return _integrar(_passo_euler, f, t0, y0, h, n_passos)


def heun(f, t0, y0, h, n_passos):
    """Método de Heun (Euler melhorado)."""
    return _integrar(_passo_heun, f, t0, y0, h, n_passos)


def rk4(f, t0, y0, h, n_passos):
    """Runge-Kutta clássico de 4a ordem.

    Exemplo (sistema de 2 EDOs):
        def f(t, y):
            return np.array([y[1], -y[0]])        # oscilador harmônico
        t, Y = rk4(f, 0.0, [1.0, 0.0], 0.01, 1000)
        posicao = Y[:, 0]
    """
    return _integrar(_passo_rk4, f, t0, y0, h, n_passos)
