"""
polinomios.py
=============
Ferramentas para polinômios e o método de Bairstow.

CONVENÇÃO DOS COEFICIENTES (a mesma da apostila de Bairstow): lista em ordem
CRESCENTE de potência,

    coef = [a0, a1, a2, ..., an]   <->   f(x) = a0 + a1 x + a2 x^2 + ... + an x^n

Exemplo: 2x^3 + 3x^2 - 4  ->  coef = [-4, 0, 3, 2]
(Atenção: o numpy.polyval usa a ordem contrária, decrescente.)

Funções:
    avaliar_polinomio        -> Horner: f(x) com n multiplicações
    derivar_polinomio        -> coeficientes de f'(x)
    somar_polinomios, multiplicar_polinomios, polinomio_de_raizes
    determinante_polinomial  -> det de uma matriz cujos elementos são polinômios
                                (gera o polinômio característico det(lambda^2 M +
                                lambda C + K) sem precisar expandir à mão)
    divisao_sintetica        -> Briot-Ruffini: f(x)/(x - t) = Q(x) + resto
    divisao_quadratica       -> f(x)/(x^2 - r x - s): coeficientes b_i
    bairstow                 -> todas as raízes (reais e complexas)
    mapa_convergencia_bairstow -> varre (r0, s0) e conta iterações (fractal)
"""

import cmath

import numpy as np

from utils.sistemas_lineares import cramer_2x2


# ---------------------------------------------------------------------------
# Operações básicas
# ---------------------------------------------------------------------------
def avaliar_polinomio(coef, x):
    """Método de Horner: f(x) = a0 + x(a1 + x(a2 + ... + x an)).

    Usa só n multiplicações e n somas (em vez de calcular cada x^i).
    Funciona com x real, complexo ou array.
    """
    resultado = 0.0
    for a in reversed(coef):        # começa pelo coeficiente de maior grau
        resultado = resultado * x + a
    return resultado


def derivar_polinomio(coef):
    """Coeficientes de f'(x):  d/dx (a_i x^i) = i a_i x^(i-1)."""
    return [i * coef[i] for i in range(1, len(coef))]


def somar_polinomios(p, q):
    """Soma termo a termo (completa com zeros o polinômio de menor grau)."""
    n = max(len(p), len(q))
    return [(p[i] if i < len(p) else 0.0) + (q[i] if i < len(q) else 0.0) for i in range(n)]


def multiplicar_polinomios(p, q):
    """Produto (convolução dos coeficientes): c_{i+j} += p_i q_j."""
    produto = [0.0] * (len(p) + len(q) - 1)
    for i, p_i in enumerate(p):
        for j, q_j in enumerate(q):
            produto[i + j] += p_i * q_j
    return produto


def polinomio_de_raizes(raizes, coef_lider=1.0):
    """Monta os coeficientes de coef_lider * (x - r1)(x - r2)...(x - rn).

    Raízes complexas devem vir em pares conjugados para o resultado ser real.
    """
    coef = [coef_lider]
    for r in raizes:
        coef = multiplicar_polinomios(coef, [-r, 1.0])     # multiplica por (x - r)
    # se todas as partes imaginárias se cancelaram, devolve números reais
    if all(abs(complex(c).imag) < 1e-9 * max(1.0, abs(c)) for c in coef):
        coef = [complex(c).real for c in coef]
    return coef


def determinante_polinomial(matriz):
    """Determinante de uma matriz quadrada cujos elementos são polinômios
    (cada elemento é uma lista de coeficientes crescentes).

    Usa expansão de Laplace pela primeira linha:
        det(M) = soma_j (-1)^j * m_0j * det(M sem linha 0 e coluna j)
    É exato (só soma e multiplica polinômios) e serve bem até ~5x5.
    """
    n = len(matriz)
    if n == 1:
        return list(matriz[0][0])
    resultado = [0.0]
    for j in range(n):
        # menor complementar: remove a linha 0 e a coluna j
        menor = [[matriz[i][c] for c in range(n) if c != j] for i in range(1, n)]
        termo = multiplicar_polinomios(matriz[0][j], determinante_polinomial(menor))
        if j % 2 == 1:
            termo = [-t for t in termo]                # sinal (-1)^j
        resultado = somar_polinomios(resultado, termo)
    return resultado


# ---------------------------------------------------------------------------
# Deflação (divisão sintética)
# ---------------------------------------------------------------------------
def divisao_sintetica(coef, t):
    """Briot-Ruffini: divide f(x) por (x - t).

    Recorrência (apostila, eq. 11):  b_n = a_n;  b_i = a_i + t b_{i+1}
    Retorna (quociente, resto): quociente = [b1, ..., bn] (grau n-1) e resto = b0.
    """
    n = len(coef) - 1
    b = [0.0] * (n + 1)
    b[n] = coef[n]
    for i in range(n - 1, -1, -1):
        b[i] = coef[i] + t * b[i + 1]
    return b[1:], b[0]


def divisao_quadratica(coef, r, s):
    """Divide f(x) por x^2 - r x - s (apostila, eq. 13):

        b_n = a_n
        b_{n-1} = a_{n-1} + r b_n
        b_i = a_i + r b_{i+1} + s b_{i+2},   i = n-2, ..., 0

    Quociente = [b2, ..., bn];  resto = b0 + b1 (x - r).
    r e s podem ser números ou arrays (para varrer muitos pares de uma vez).
    """
    n = len(coef) - 1
    b = [0.0] * (n + 1)
    b[n] = coef[n]
    b[n - 1] = coef[n - 1] + r * b[n]
    for i in range(n - 2, -1, -1):
        b[i] = coef[i] + r * b[i + 1] + s * b[i + 2]
    return b


def _correcao_bairstow(coef, r, s):
    """Um passo de Newton-Raphson do Bairstow: devolve (dr, ds, b).

    As derivadas parciais saem de uma 2a divisão sintética (eq. 16-17):
        c_n = b_n;  c_{n-1} = b_{n-1} + r c_n;  c_i = b_i + r c_{i+1} + s c_{i+2}
        c1 = db0/dr,  c2 = db0/ds = db1/dr,  c3 = db1/ds
    Sistema 2x2 (eq. 15), resolvido por Cramer:
        c2 dr + c3 ds = -b1
        c1 dr + c2 ds = -b0
    """
    b = divisao_quadratica(coef, r, s)
    n = len(coef) - 1
    c = [0.0] * (n + 1)
    c[n] = b[n]
    c[n - 1] = b[n - 1] + r * c[n]
    for i in range(n - 2, 0, -1):
        c[i] = b[i] + r * c[i + 1] + s * c[i + 2]
    dr, ds = cramer_2x2(c[2], c[3], c[1], c[2], -b[1], -b[0])
    return dr, ds, b


def _raizes_quadratica(r, s):
    """Raízes de x^2 - r x - s = 0:  x = (r +- sqrt(r^2 + 4s)) / 2 (eq. 19)."""
    delta = cmath.sqrt(r * r + 4.0 * s)
    return [_simplificar((r + delta) / 2.0), _simplificar((r - delta) / 2.0)]


def _simplificar(z):
    """Devolve float quando a parte imaginária é desprezível."""
    z = complex(z)
    return z.real if abs(z.imag) <= 1e-12 * max(1.0, abs(z.real)) else z


# ---------------------------------------------------------------------------
# Método de Bairstow
# ---------------------------------------------------------------------------
def bairstow(coef, r0, s0, tol=1e-10, max_iter=200):
    """Encontra TODAS as raízes de um polinômio pelo método de Bairstow.

    Algoritmo (fluxograma da apostila):
      1. chuta (r, s) e divide f(x) por x^2 - r x - s;
      2. corrige (r, s) com Newton-Raphson até o resto (b0, b1) zerar,
         isto é, até max(|dr|, |ds|) < tol;
      3. tira duas raízes de x^2 - r x - s pela fórmula quadrática;
      4. repete no quociente (deflação), reaproveitando (r, s) como chute;
      5. sobrou grau 2 -> fórmula quadrática; sobrou grau 1 -> x = -a0/a1
         (que é o x = -s/r da apostila, escrevendo o resto como r x + s).

    Parâmetros
    ----------
    coef   : [a0, a1, ..., an] (ordem CRESCENTE), an != 0
    r0, s0 : chute inicial do fator quadrático
    Retorna (raizes, historico); historico tem uma linha por fator encontrado:
    fator, iteracoes, r, s, raiz_1, raiz_2.
    """
    a = [float(c) for c in coef]
    if a[-1] == 0.0:
        raise ValueError("Bairstow: o coeficiente de maior grau não pode ser zero.")
    raizes, historico = [], []
    r, s = float(r0), float(s0)
    fator = 0
    while len(a) - 1 >= 3:                          # grau >= 3: usa Bairstow
        fator += 1
        for iteracao in range(1, max_iter + 1):
            dr, ds, b = _correcao_bairstow(a, r, s)
            r, s = r + dr, s + ds
            if not (np.isfinite(r) and np.isfinite(s)):
                raise ArithmeticError("Bairstow divergiu; tente outro (r0, s0).")
            if max(abs(dr), abs(ds)) < tol:
                break
        else:
            print(f"AVISO: Bairstow não convergiu no fator {fator} em {max_iter} iterações.")
        par = _raizes_quadratica(r, s)
        raizes.extend(par)
        historico.append({"fator": fator, "iteracoes": iteracao, "r": r, "s": s,
                          "raiz_1": par[0], "raiz_2": par[1]})
        # deflação: o quociente da divisão pelo fator convergido
        a = divisao_quadratica(a, r, s)[2:]
    # o que sobrou tem grau 2 ou 1: solução direta
    if len(a) - 1 == 2:
        raizes.extend(_raizes_quadratica(-a[1] / a[2], -a[0] / a[2]))
    elif len(a) - 1 == 1:
        raizes.append(-a[0] / a[1])
    return raizes, historico


def mapa_convergencia_bairstow(coef, R0, S0, tol=1e-8, max_iter=100, limite=1e6):
    """Varre uma grade de chutes (r0, s0) e registra quantas iterações o
    Bairstow leva para achar o PRIMEIRO fator quadrático (fractal de Bairstow).

    Usa a MESMA correção _correcao_bairstow, só que aplicada a arrays: todas as
    células da grade avançam juntas (muito mais rápido que um laço por ponto).

    Retorna
    -------
    iteracoes : array int com o número de iterações (-1 = não convergiu)
    r, s      : arrays com o fator final de cada chute (para pintar as bacias)
    """
    r = np.array(R0, dtype=float)
    s = np.array(S0, dtype=float)
    iteracoes = np.full(r.shape, -1, dtype=int)
    ativo = np.ones(r.shape, dtype=bool)            # células ainda iterando
    with np.errstate(all="ignore"):                 # divisões por 0 viram inf/nan
        for k in range(1, max_iter + 1):
            dr, ds, _ = _correcao_bairstow(coef, r, s)
            # só atualiza quem ainda está ativo
            r = np.where(ativo, r + dr, r)
            s = np.where(ativo, s + ds, s)
            passo = np.maximum(np.abs(dr), np.abs(ds))
            convergiu = ativo & (passo < tol)
            iteracoes[convergiu] = k
            divergiu = ativo & ~(np.isfinite(r) & np.isfinite(s) & np.isfinite(passo)
                                 & (np.abs(r) < limite) & (np.abs(s) < limite))
            ativo &= ~(convergiu | divergiu)
            if not ativo.any():
                break
    return iteracoes, r, s
