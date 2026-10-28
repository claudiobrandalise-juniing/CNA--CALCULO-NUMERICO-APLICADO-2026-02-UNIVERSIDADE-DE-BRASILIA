"""
otimizacao.py
=============
Otimização UNIDIMENSIONAL sem restrições (Chapra, cap. 13).

Todos os métodos MINIMIZAM f(x). Para MAXIMIZAR g(x), minimize f(x) = -g(x).

    razao_aurea              -> encolhe [a, b] mantendo dois pontos internos na
                                razão áurea; reaproveita uma avaliação por
                                iteração (só 1 avaliação nova de f por passo)
    interpolacao_quadratica  -> ajusta uma parábola a 3 pontos e pula para o
                                vértice; converge bem mais rápido perto do
                                mínimo, com proteções (passo áureo de reserva)
    newton_otimizacao        -> Newton aplicado a f'(x) = 0:
                                x_{k+1} = x_k - f'(x_k)/f''(x_k)

    FuncaoContada            -> "embrulha" f e conta quantas vezes ela foi
                                avaliada (medida de custo computacional)

O número de AVALIAÇÕES de f é a medida justa de custo: num problema real,
cada avaliação pode ser uma simulação inteira (minutos ou horas).
"""

import math

from utils.derivadas import derivada_central, segunda_derivada_central

# razão áurea tau = (sqrt(5) - 1)/2 ~ 0.618034
TAU = (math.sqrt(5.0) - 1.0) / 2.0


class FuncaoContada:
    """Função que conta as próprias avaliações.

    Exemplo:
        C = FuncaoContada(custo_total)
        x_otimo, historico = razao_aurea(C, 0.0, 0.12)
        print(C.n_avaliacoes)
    """

    def __init__(self, funcao):
        self.funcao = funcao
        self.n_avaliacoes = 0

    def __call__(self, x):
        self.n_avaliacoes += 1
        return self.funcao(x)


def razao_aurea(f, a, b, tol=1e-6, max_iter=200):
    """Busca pela razão áurea (golden-section search) do mínimo de f em [a, b].

    Pontos internos:  x1 = b - tau (b - a),   x2 = a + tau (b - a)
    Se f(x1) > f(x2): o mínimo NÃO está em [a, x1]  -> a = x1
    Senão:            o mínimo NÃO está em [x2, b]  -> b = x2
    Por causa da razão áurea, um dos pontos internos antigos vira um dos novos,
    então cada iteração custa uma única avaliação de f.

    Retorna x_otimo = (a + b)/2 e o historico
    (colunas: k, a, b, x1, x2, f(x1), f(x2), largura).
    Hipótese: f é unimodal em [a, b] (um único mínimo).
    """
    x1 = b - TAU * (b - a)
    x2 = a + TAU * (b - a)
    f1, f2 = f(x1), f(x2)
    historico = []
    k = 0
    while abs(b - a) > tol and k < max_iter:
        if f1 > f2:
            # descarta [a, x1]; o antigo x2 vira o novo x1 (reaproveita f2)
            a = x1
            x1, f1 = x2, f2
            x2 = a + TAU * (b - a)
            f2 = f(x2)
        else:
            # descarta [x2, b]; o antigo x1 vira o novo x2 (reaproveita f1)
            b = x2
            x2, f2 = x1, f1
            x1 = b - TAU * (b - a)
            f1 = f(x1)
        k += 1
        historico.append({"k": k, "a": a, "b": b, "x1": x1, "x2": x2,
                          "f(x1)": f1, "f(x2)": f2, "largura": b - a})
    if abs(b - a) > tol:
        print(f"AVISO: razão áurea parou em {max_iter} iterações sem atingir a tolerância.")
    return (a + b) / 2.0, historico


def interpolacao_quadratica(f, x1, x2, x3, tol=1e-6, max_iter=200):
    """Interpolação quadrática (parabólica) para o mínimo de f.

    A parábola p(x) = A x^2 + B x + D pelos 3 pontos é obtida por diferenças
    divididas (sem resolver sistema 3x3):
        d1 = (f2 - f1)/(x2 - x1)
        A  = [(f3 - f2)/(x3 - x2) - d1] / (x3 - x1)
        B  = d1 - A (x1 + x2)
        x_novo = -B / (2A)          (vértice; é mínimo se A > 0)

    Atualização da tripla (mantém x1 < x2 < x3 com f2 menor que f1 e f3):
        x_novo < x2 e f_novo <  f2 -> (x1, x_novo, x2)
        x_novo < x2 e f_novo >= f2 -> (x_novo, x2, x3)
        x_novo > x2 e f_novo <  f2 -> (x2, x_novo, x3)
        x_novo > x2 e f_novo >= f2 -> (x1, x2, x_novo)

    PROTEÇÕES: se A <= 0 (parábola sem mínimo), se o vértice cair fora de
    (x1, x3) ou praticamente em cima de um ponto já existente, faz-se um passo
    da razão áurea dentro do maior subintervalo. Parada: x3 - x1 < tol.

    Retorna x_otimo (= x2 final) e historico
    (colunas: k, x1, x2, x3, x_novo, f(x_novo), passo, largura, variacao_x2).
    """
    if not x1 < x2 < x3:
        raise ValueError("Interpolação quadrática: exige x1 < x2 < x3.")
    f1, f2, f3 = f(x1), f(x2), f(x3)
    if not (f2 < f1 and f2 < f3):
        raise ValueError("Interpolação quadrática: exige f(x2) < f(x1) e f(x2) < f(x3) "
                         "(o mínimo precisa estar cercado). Ajuste os pontos.")
    historico = []
    k = 0
    while x3 - x1 >= tol and k < max_iter:
        k += 1
        # ---- vértice da parábola por diferenças divididas --------------------
        d1 = (f2 - f1) / (x2 - x1)
        A = ((f3 - f2) / (x3 - x2) - d1) / (x3 - x1)
        B = d1 - A * (x1 + x2)
        passo = "parabólico"
        # ---- proteções ------------------------------------------------------
        if A <= 1e-14 * max(abs(f1), abs(f2), abs(f3), 1.0):
            passo = "áureo (A<=0)"
        else:
            x_novo = -B / (2.0 * A)
            if not x1 < x_novo < x3:
                passo = "áureo (fora)"
            elif min(abs(x_novo - x1), abs(x_novo - x2), abs(x_novo - x3)) < tol / 2.0:
                passo = "áureo (coincide)"
        if passo != "parabólico":
            # passo da razão áurea no maior dos dois subintervalos
            if x3 - x2 > x2 - x1:
                x_novo = x2 + (1.0 - TAU) * (x3 - x2)
            else:
                x_novo = x2 - (1.0 - TAU) * (x2 - x1)
        f_novo = f(x_novo)
        x2_antigo = x2
        # ---- atualização da tripla -------------------------------------------
        if x_novo < x2:
            if f_novo < f2:
                x3, f3, x2, f2 = x2, f2, x_novo, f_novo
            else:
                x1, f1 = x_novo, f_novo
        else:
            if f_novo < f2:
                x1, f1, x2, f2 = x2, f2, x_novo, f_novo
            else:
                x3, f3 = x_novo, f_novo
        historico.append({"k": k, "x1": x1, "x2": x2, "x3": x3, "x_novo": x_novo,
                          "f(x_novo)": f_novo, "passo": passo, "largura": x3 - x1,
                          "variacao_x2": abs(x2 - x2_antigo)})
    if x3 - x1 >= tol:
        print(f"AVISO: interpolação quadrática parou em {max_iter} iterações.")
    return x2, historico


def newton_otimizacao(f, x0, df=None, d2f=None, tol=1e-8, max_iter=100):
    """Newton para otimização: acha o zero de f'(x).

        x_{k+1} = x_k - f'(x_k) / f''(x_k)

    df e d2f são as derivadas analíticas; se forem None, usa diferenças
    centrais. f''(x*) > 0 indica mínimo; f''(x*) < 0 indica máximo.
    Retorna x_otimo, historico (colunas: k, x_k, f'(x_k), f''(x_k), erro).
    """
    primeira = df if df is not None else (lambda x: derivada_central(f, x))
    segunda = d2f if d2f is not None else (lambda x: segunda_derivada_central(f, x))
    x = x0
    historico = []
    for k in range(1, max_iter + 1):
        d1, d2 = primeira(x), segunda(x)
        if d2 == 0:
            raise ZeroDivisionError("Newton (otimização): f''(x) = 0.")
        x_novo = x - d1 / d2
        erro = abs(x_novo - x)
        historico.append({"k": k, "x_k": x, "f'(x_k)": d1, "f''(x_k)": d2, "erro": erro})
        x = x_novo
        if erro < tol:
            if segunda(x) < 0:
                print("AVISO: f''(x*) < 0, o ponto encontrado é um MÁXIMO.")
            return x, historico
    print(f"AVISO: Newton (otimização) NÃO convergiu em {max_iter} iterações.")
    return x, historico
