"""
derivadas.py
============
Derivadas numéricas por diferenças finitas (Chapra, cap. 4 e 23).

Usadas quando não se tem a derivada analítica (ex.: Newton com f caixa-preta,
Jacobiana numérica, Newton para otimização).

    progressiva:  f'(x) ~ [f(x+h) - f(x)] / h               erro O(h)
    central:      f'(x) ~ [f(x+h) - f(x-h)] / (2h)          erro O(h^2)
    2a central:   f''(x) ~ [f(x+h) - 2 f(x) + f(x-h)] / h^2  erro O(h^2)

Escolha do passo: h muito grande -> erro de truncamento; h muito pequeno ->
erro de arredondamento. Um bom padrão é h ~ 1e-5 * max(1, |x|) para a
central e h ~ 1e-7 * max(1, |x|) para a progressiva.
"""


def _passo(x, relativo):
    """Passo h proporcional à escala de x (evita h desprezível para x grande)."""
    return relativo * max(1.0, abs(x))


def derivada_progressiva(f, x, h=None):
    """f'(x) por diferença progressiva (1 avaliação extra de f)."""
    h = _passo(x, 1e-7) if h is None else h
    return (f(x + h) - f(x)) / h


def derivada_central(f, x, h=None):
    """f'(x) por diferença central (mais precisa que a progressiva)."""
    h = _passo(x, 1e-5) if h is None else h
    return (f(x + h) - f(x - h)) / (2.0 * h)


def segunda_derivada_central(f, x, h=None):
    """f''(x) por diferença central de segunda ordem."""
    h = _passo(x, 1e-4) if h is None else h
    return (f(x + h) - 2.0 * f(x) + f(x - h)) / (h * h)
