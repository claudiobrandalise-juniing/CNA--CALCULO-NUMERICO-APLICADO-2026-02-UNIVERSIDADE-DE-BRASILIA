"""
raizes.py
=========
Métodos para encontrar zeros de funções de UMA variável, f(x) = 0.

MÉTODOS INTERVALARES (precisam de [xl, xu] com f(xl)*f(xu) < 0; sempre convergem
se f for contínua, porque a raiz nunca sai do intervalo):
    bissecao                  -> xr = (xl + xu)/2
    falsa_posicao             -> xr = reta entre (xl, f(xl)) e (xu, f(xu))
    falsa_posicao_modificada  -> idem, mas divide por 2 o f do extremo que
                                 fica "preso" 2 vezes seguidas (acelera)

MÉTODOS ABERTOS (precisam só de chute(s) inicial(is); convergem mais rápido,
mas podem divergir se o chute for ruim):
    ponto_fixo      -> x_{k+1} = g(x_k)                    (convergência linear)
    newton_raphson  -> x_{k+1} = x_k - f(x_k)/f'(x_k)       (quadrática)
    secante         -> Newton com f' trocada pela inclinação da secante (~1.62)
    muller          -> parábola por 3 pontos; acha raízes complexas também

AUXILIAR:
    varredura_mudanca_sinal -> "método gráfico" automático: procura subintervalos
                               onde f troca de sinal (bom para achar chutes)

CRITÉRIOS DE PARADA (argumento 'criterio'):
    "funcao"   -> |f(x)| < tol                       (resíduo pequeno)
    "passo"    -> |x_novo - x_antigo| < tol          (variação absoluta)
    "relativo" -> |(x_novo - x_antigo)/x_novo|*100 < tol   (erro aproximado, %)

TODOS os métodos devolvem (raiz, historico), em que historico é uma lista de
dicionários (uma linha por iteração), pronta para Relatorio.tabela().
"""

import cmath
import math


# ---------------------------------------------------------------------------
# Funções internas compartilhadas (princípio DRY do Manual de Boas Práticas)
# ---------------------------------------------------------------------------
def _calcular_erro(criterio, x_novo, x_antigo, f_novo):
    """Calcula o erro de acordo com o critério escolhido."""
    if criterio == "funcao":
        return abs(f_novo)
    # para 'passo' e 'relativo' é preciso um valor anterior
    if x_antigo is None:
        return math.inf          # primeira iteração: ainda não há como comparar
    if criterio == "passo":
        return abs(x_novo - x_antigo)
    if criterio == "relativo":
        if x_novo == 0:          # evita divisão por zero
            return abs(x_novo - x_antigo)
        return abs((x_novo - x_antigo) / x_novo) * 100.0
    raise ValueError("criterio deve ser 'funcao', 'passo' ou 'relativo'")


def _avisar_nao_convergiu(nome_metodo, max_iter):
    """Mensagem padrão quando o limite de iterações é atingido."""
    print(f"AVISO: {nome_metodo} NÃO convergiu em {max_iter} iterações. "
          "Verifique o chute inicial, a tolerância ou o intervalo.")


# ---------------------------------------------------------------------------
# Métodos intervalares
# ---------------------------------------------------------------------------
def _metodo_intervalar(f, xl, xu, tol, max_iter, criterio, regra, modificada, nome):
    """Laço comum da bisseção e da falsa posição (normal ou modificada).

    regra = "bissecao"      -> xr = (xl + xu) / 2
    regra = "falsa_posicao" -> xr = xu - f(xu)*(xl - xu)/(f(xl) - f(xu))
    """
    # avalia os extremos uma única vez; depois os valores são reaproveitados
    fl = f(xl)
    fu = f(xu)
    # a troca de sinal garante (para f contínua) pelo menos uma raiz em [xl, xu]
    if fl * fu > 0:
        raise ValueError(f"{nome}: f(xl) = {fl:.4g} e f(xu) = {fu:.4g} têm o mesmo "
                         "sinal; o intervalo não garante uma raiz.")

    historico = []
    xr = None                    # ainda não existe aproximação anterior
    preso_l = preso_u = 0        # contadores da falsa posição modificada

    for k in range(1, max_iter + 1):
        xr_antigo = xr
        # --- nova aproximação dentro do intervalo -------------------------
        if regra == "bissecao":
            xr = (xl + xu) / 2.0
        else:
            xr = xu - fu * (xl - xu) / (fl - fu)
        fr = f(xr)
        erro = _calcular_erro(criterio, xr, xr_antigo, fr)
        historico.append({"k": k, "xl": xl, "xu": xu, "xr": xr, "f(xr)": fr, "erro": erro})

        # --- parada: raiz exata ou erro abaixo da tolerância ---------------
        if fr == 0.0 or erro < tol:
            return xr, historico

        # --- atualização do intervalo: mantém a troca de sinal -------------
        if fl * fr < 0:
            # raiz está entre xl e xr -> xr vira o novo limite superior
            xu, fu = xr, fr
            preso_u, preso_l = 0, preso_l + 1
            # (modificada) xl ficou parado 2 vezes: reduz f(xl) pela metade
            if modificada and preso_l >= 2:
                fl = fl / 2.0
        else:
            # raiz está entre xr e xu -> xr vira o novo limite inferior
            xl, fl = xr, fr
            preso_l, preso_u = 0, preso_u + 1
            if modificada and preso_u >= 2:
                fu = fu / 2.0

    _avisar_nao_convergiu(nome, max_iter)
    return xr, historico


def bissecao(f, xl, xu, tol=1e-6, max_iter=200, criterio="funcao"):
    """Método da bisseção.

    Parâmetros
    ----------
    f        : função f(x)
    xl, xu   : limites inferior e superior com f(xl)*f(xu) < 0
    tol      : tolerância do critério de parada
    max_iter : limite de iterações (evita laço infinito)
    criterio : "funcao", "passo" ou "relativo"

    Retorna
    -------
    raiz, historico  (colunas: k, xl, xu, xr, f(xr), erro)

    Obs.: o erro máximo após k iterações é (xu - xl)/2^k, então o número de
    iterações para um erro absoluto E é k = log2((xu - xl)/E).
    """
    return _metodo_intervalar(f, xl, xu, tol, max_iter, criterio,
                              "bissecao", False, "Bisseção")


def falsa_posicao(f, xl, xu, tol=1e-6, max_iter=200, criterio="funcao"):
    """Método da falsa posição (regula falsi).

    A nova aproximação é o ponto onde a reta que liga (xl, f(xl)) e (xu, f(xu))
    cruza o eixo x:  xr = xu - f(xu)*(xl - xu)/(f(xl) - f(xu)).
    Retorna raiz, historico (mesmas colunas da bisseção).
    """
    return _metodo_intervalar(f, xl, xu, tol, max_iter, criterio,
                              "falsa_posicao", False, "Falsa posição")


def falsa_posicao_modificada(f, xl, xu, tol=1e-6, max_iter=200, criterio="funcao"):
    """Falsa posição modificada (Chapra, fig. 5.15).

    Se um dos extremos fica parado em duas iterações seguidas, o valor de f
    nesse extremo é dividido por 2. Isso evita a convergência lenta da falsa
    posição em funções muito curvas (um lado do intervalo "trava").
    """
    return _metodo_intervalar(f, xl, xu, tol, max_iter, criterio,
                              "falsa_posicao", True, "Falsa posição modificada")


# ---------------------------------------------------------------------------
# Métodos abertos
# ---------------------------------------------------------------------------
def ponto_fixo(g, x0, tol=1e-6, max_iter=500, criterio="passo"):
    """Iteração de ponto fixo: reescreve f(x) = 0 como x = g(x).

    Converge se |g'(x)| < 1 perto da raiz (quanto menor, mais rápido).
    Retorna raiz, historico (colunas: k, x_k, x_k+1, erro).
    """
    historico = []
    x = x0
    for k in range(1, max_iter + 1):
        x_novo = g(x)
        # no ponto fixo o "resíduo" é x - g(x) (só calculado se for o critério)
        erro = _calcular_erro(criterio, x_novo, x,
                              x_novo - g(x_novo) if criterio == "funcao" else None)
        historico.append({"k": k, "x_k": x, "x_k+1": x_novo, "erro": erro})
        x = x_novo
        if erro < tol:
            return x, historico
    _avisar_nao_convergiu("Ponto fixo", max_iter)
    return x, historico


def newton_raphson(f, df, x0, tol=1e-6, max_iter=100, criterio="passo"):
    """Método de Newton-Raphson.

    Ideia: aproxima f por sua reta tangente em x_k e toma o zero dessa reta:
        x_{k+1} = x_k - f(x_k) / f'(x_k)

    Parâmetros
    ----------
    f, df : função e sua derivada analítica (para derivada numérica, use
            utils.derivadas.derivada_central)
    x0    : chute inicial

    Retorna raiz, historico (colunas: k, x_k, f(x_k), f'(x_k), x_k+1, erro)
    """
    historico = []
    x = x0
    for k in range(1, max_iter + 1):
        fx = f(x)
        dfx = df(x)
        # tangente horizontal: o passo de Newton não existe
        if dfx == 0:
            raise ZeroDivisionError(f"Newton-Raphson: f'(x) = 0 em x = {x}. Mude o chute.")
        x_novo = x - fx / dfx
        erro = _calcular_erro(criterio, x_novo, x, f(x_novo) if criterio == "funcao" else None)
        historico.append({"k": k, "x_k": x, "f(x_k)": fx, "f'(x_k)": dfx,
                          "x_k+1": x_novo, "erro": erro})
        x = x_novo
        if erro < tol:
            return x, historico
    _avisar_nao_convergiu("Newton-Raphson", max_iter)
    return x, historico


def secante(f, x_anterior, x_atual, tol=1e-6, max_iter=100, criterio="passo"):
    """Método da secante.

    Igual a Newton, mas a derivada é trocada pela inclinação da reta secante
    pelos dois últimos pontos (não precisa de f'):
        x_{k+1} = x_k - f(x_k) * (x_k - x_{k-1}) / (f(x_k) - f(x_{k-1}))

    Útil quando f vem de um experimento, de uma tabela ou de outro programa
    (caixa-preta). Retorna raiz, historico
    (colunas: k, x_k-1, x_k, f(x_k), x_k+1, erro).
    """
    historico = []
    f_anterior = f(x_anterior)       # cada ponto é avaliado uma única vez
    for k in range(1, max_iter + 1):
        f_atual = f(x_atual)
        # secante horizontal: não cruza o eixo x
        if f_atual == f_anterior:
            raise ZeroDivisionError("Secante: f(x_k) = f(x_k-1); a secante é horizontal.")
        x_novo = x_atual - f_atual * (x_atual - x_anterior) / (f_atual - f_anterior)
        erro = _calcular_erro(criterio, x_novo, x_atual,
                              f(x_novo) if criterio == "funcao" else None)
        historico.append({"k": k, "x_k-1": x_anterior, "x_k": x_atual, "f(x_k)": f_atual,
                          "x_k+1": x_novo, "erro": erro})
        # desloca a janela de dois pontos
        x_anterior, f_anterior = x_atual, f_atual
        x_atual = x_novo
        if erro < tol:
            return x_atual, historico
    _avisar_nao_convergiu("Secante", max_iter)
    return x_atual, historico


def muller(f, x0, x1, x2, tol=1e-8, max_iter=100, criterio="passo"):
    """Método de Müller.

    Passa uma parábola p(x) = a(x - x2)^2 + b(x - x2) + c pelos três últimos
    pontos e usa a raiz da parábola mais próxima de x2 como nova aproximação.
    Como usa raiz quadrada complexa, também encontra raízes complexas.

    Passos (iguais à apostila):
        h0 = x1 - x0,  h1 = x2 - x1
        d0 = (f1 - f0)/h0,  d1 = (f2 - f1)/h1        (inclinações das secantes)
        a  = (d1 - d0)/(h1 + h0),  b = a*h1 + d1,  c = f2
        D  = sqrt(b^2 - 4ac)
        E  = b + D se |b + D| > |b - D|, senão b - D  (evita cancelamento)
        x3 = x2 - 2c/E

    Retorna raiz (float se for real, complex se não), historico
    (colunas: k, x0, x1, x2, x3, f(x3), erro).
    """
    historico = []
    for k in range(1, max_iter + 1):
        f0, f1, f2 = f(x0), f(x1), f(x2)
        h0, h1 = x1 - x0, x2 - x1
        delta0, delta1 = (f1 - f0) / h0, (f2 - f1) / h1
        a = (delta1 - delta0) / (h1 + h0)
        b = a * h1 + delta1
        c = f2
        # raiz complexa: o discriminante pode ser negativo
        D = cmath.sqrt(b * b - 4 * a * c)
        # escolhe o denominador de maior módulo (menor erro de arredondamento)
        E = b + D if abs(b + D) > abs(b - D) else b - D
        if E == 0:
            raise ZeroDivisionError("Müller: denominador nulo; mude os pontos iniciais.")
        x3 = x2 + (-2 * c / E)
        # se a parte imaginária é desprezível, trabalha com número real
        if abs(x3.imag) <= 1e-14 * max(1.0, abs(x3.real)):
            x3 = x3.real
        f3 = f(x3)
        erro = _calcular_erro(criterio, x3, x2, f3)
        historico.append({"k": k, "x0": x0, "x1": x1, "x2": x2, "x3": x3,
                          "f(x3)": f3, "erro": erro})
        # descarta o ponto mais antigo: (x0, x1, x2) <- (x1, x2, x3)
        x0, x1, x2 = x1, x2, x3
        if erro < tol:
            return x3, historico
    _avisar_nao_convergiu("Müller", max_iter)
    return x2, historico


# ---------------------------------------------------------------------------
# Auxiliar: método gráfico automatizado
# ---------------------------------------------------------------------------
def varredura_mudanca_sinal(f, a, b, n_subintervalos=1000):
    """Divide [a, b] em partes iguais e devolve os subintervalos onde f troca de
    sinal. Cada par (xl, xu) devolvido serve de intervalo inicial para a
    bisseção/falsa posição ou de chute para Newton/secante.

    Cuidado: duas raízes muito próximas dentro do mesmo subintervalo não trocam
    o sinal e passam despercebidas; e assíntotas verticais (ex.: tangente)
    também trocam o sinal sem serem raízes. Sempre confira com um gráfico.
    """
    intervalos = []
    passo = (b - a) / n_subintervalos
    x_esq, f_esq = a, f(a)
    for i in range(1, n_subintervalos + 1):
        x_dir = a + i * passo
        f_dir = f(x_dir)
        if f_esq * f_dir <= 0:
            intervalos.append((x_esq, x_dir))
        x_esq, f_esq = x_dir, f_dir
    return intervalos
