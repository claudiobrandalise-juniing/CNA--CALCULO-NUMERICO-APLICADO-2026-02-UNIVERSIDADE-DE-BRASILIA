"""
testar_tudo.py
==============
Confere TODOS os métodos da biblioteca utils/ contra respostas conhecidas.

Nenhum resolvedor pronto é usado: as respostas de referência são FABRICADAS
(ex.: escolhe-se x_ref e calcula-se b = A x_ref; escolhem-se as raízes e
monta-se o polinômio) ou vêm de soluções analíticas.

Como rodar (a partir da raiz do repositório):
    python testes/testar_tudo.py

Saída esperada: uma linha [OK] por teste e "TODOS OS TESTES PASSARAM".
"""

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from utils import (edo, edp, otimizacao, polinomios, raizes, sistemas_iterativos,  # noqa: E402
                   sistemas_lineares as sl, sistemas_nao_lineares)

falhas = []


def verificar(nome, condicao, detalhe=""):
    """Registra e imprime o resultado de um teste."""
    print(f"[{'OK' if condicao else 'FALHOU'}] {nome} {detalhe}")
    if not condicao:
        falhas.append(nome)


# ---------------------------------------------------------------------------
# Raízes de uma variável
# ---------------------------------------------------------------------------
def f_wallis(x):
    """x^3 - 2x - 5 (exemplo clássico de Wallis), raiz 2.0945514815423265."""
    return x ** 3 - 2 * x - 5


RAIZ_WALLIS = 2.0945514815423265
for metodo in (raizes.bissecao, raizes.falsa_posicao, raizes.falsa_posicao_modificada):
    x, h = metodo(f_wallis, 2.0, 3.0, tol=1e-12, criterio="funcao")
    verificar(metodo.__name__, abs(x - RAIZ_WALLIS) < 1e-10, f"({len(h)} iterações)")

x, h = raizes.newton_raphson(f_wallis, lambda x: 3 * x ** 2 - 2, 3.0, tol=1e-12)
verificar("newton_raphson", abs(x - RAIZ_WALLIS) < 1e-12, f"({len(h)} iterações)")
x, h = raizes.secante(f_wallis, 2.0, 3.0, tol=1e-12)
verificar("secante", abs(x - RAIZ_WALLIS) < 1e-12, f"({len(h)} iterações)")
x, h = raizes.muller(f_wallis, 1.0, 2.0, 3.0, tol=1e-12)
verificar("muller (raiz real)", abs(x - RAIZ_WALLIS) < 1e-12, f"({len(h)} iterações)")
x, h = raizes.muller(lambda z: z ** 3 - 1, 0.0, 0.5j, 1.0j, tol=1e-12)
verificar("muller (raiz complexa de z^3 = 1)", abs(x ** 3 - 1) < 1e-10, f"(z = {x:.6f})")
x, h = raizes.ponto_fixo(lambda x: math.exp(-x), 0.0, tol=1e-12)
verificar("ponto_fixo (x = e^-x)", abs(x - 0.5671432904097838) < 1e-10, f"({len(h)} iterações)")

# falsa posição comum é lenta em x^10 - 1; a modificada deve ser bem mais rápida
_, h_fp = raizes.falsa_posicao(lambda x: x ** 10 - 1, 0.0, 1.3, tol=1e-8, max_iter=5000)
_, h_fpm = raizes.falsa_posicao_modificada(lambda x: x ** 10 - 1, 0.0, 1.3, tol=1e-8)
verificar("falsa posição modificada acelera", len(h_fpm) < len(h_fp),
          f"({len(h_fpm)} x {len(h_fp)} iterações)")
intervalos = raizes.varredura_mudanca_sinal(math.sin, 0.5, 10.0, 200)
verificar("varredura_mudanca_sinal (raízes de sin em [0.5, 10])", len(intervalos) == 3)

# ---------------------------------------------------------------------------
# Polinômios e Bairstow
# ---------------------------------------------------------------------------
quociente, resto = polinomios.divisao_sintetica([-4, 0, 3, 2], -1.0)   # (2x^3+3x^2-4)/(x+1)
verificar("divisao_sintetica (Briot-Ruffini)", quociente == [-1.0, 1.0, 2.0] and resto == -3.0)
lam2 = [1.0, 0.0, 1.0]                                                  # lambda^2 + 1
detp = polinomios.determinante_polinomial([[lam2, [-1.0]], [[-1.0], lam2]])
verificar("determinante_polinomial", np.allclose(detp, [0, 0, 2, 0, 1]))

raizes_ref = [-3.0, -1.0, 2.0, 4.0, 5.0, 1 + 2j, 1 - 2j]
coef = polinomios.polinomio_de_raizes(raizes_ref)
encontradas, h = polinomios.bairstow(coef, -1.0, -1.0, tol=1e-12)
erro_max = max(min(abs(r - e) for e in encontradas) for r in raizes_ref)
verificar("bairstow (7a ordem, raízes conhecidas)", erro_max < 1e-8, f"(erro máx {erro_max:.1e})")
coef_chapra = [1.25, -3.875, 2.125, 2.75, -3.5, 1.0]                    # Chapra, ex. 7.3
encontradas, _ = polinomios.bairstow(coef_chapra, -1.0, -1.0)
ref = [2.0, -1.0, 0.5, 1 + 0.5j, 1 - 0.5j]
erro_max = max(min(abs(r - e) for e in encontradas) for r in ref)
verificar("bairstow (exemplo do Chapra)", erro_max < 1e-8)
iters, r_fim, s_fim = polinomios.mapa_convergencia_bairstow(
    coef_chapra, np.array([[-1.0, 0.5]]), np.array([[-1.0, 0.3]]))
verificar("mapa_convergencia_bairstow (vetorizado)", np.all(iters > 0))

# ---------------------------------------------------------------------------
# Sistemas lineares diretos
# ---------------------------------------------------------------------------
gerador = np.random.default_rng(2026)
n = 6
A = gerador.uniform(-1, 1, (n, n)) + n * np.eye(n)      # bem condicionada
x_ref = gerador.uniform(-10, 10, n)
b = A @ x_ref                                           # b fabricado
for nome, solver in (("gauss_ingenua", sl.gauss_ingenua),
                     ("gauss_pivoteamento", sl.gauss_pivoteamento),
                     ("gauss_jordan", sl.gauss_jordan),
                     ("cramer", sl.cramer)):
    x = solver(A, b)
    verificar(nome, np.max(np.abs(x - x_ref)) < 1e-10)

# pivoteamento resolve um sistema com pivô inicial nulo
A_pivo_nulo = np.array([[0.0, 2.0, 1.0], [1.0, 1.0, 1.0], [2.0, 1.0, 3.0]])
x = sl.gauss_pivoteamento(A_pivo_nulo, A_pivo_nulo @ np.array([1.0, 2.0, 3.0]))
verificar("gauss_pivoteamento com a11 = 0", np.allclose(x, [1, 2, 3]))
try:
    sl.gauss_ingenua(A_pivo_nulo, [1.0, 1.0, 1.0])
    verificar("gauss_ingenua detecta pivô nulo", False)
except ZeroDivisionError:
    verificar("gauss_ingenua detecta pivô nulo", True)

verificar("determinante 3x3", abs(sl.determinante([[2, -1, 0], [-1, 2, -1], [0, -1, 2]]) - 4) < 1e-12)

L, U = sl.decomposicao_lu(A)
verificar("decomposicao_lu (L U = A)", np.max(np.abs(L @ U - A)) < 1e-12)
verificar("resolver_lu", np.max(np.abs(sl.resolver_lu(L, U, b) - x_ref)) < 1e-10)
Lc, Uc = sl.decomposicao_crout(A)
verificar("decomposicao_crout (L U = A, diag(U) = 1)",
          np.max(np.abs(Lc @ Uc - A)) < 1e-12 and np.allclose(np.diag(Uc), 1))
C = sl.inversa_lu(A)
verificar("inversa_lu (A A^-1 = I)", np.max(np.abs(A @ C - np.eye(n))) < 1e-12)

e, f, g = sl.tridiagonal_para_vetores(sl.vetores_para_tridiagonal(
    [0, -1, -1, -1], [2.04, 2.04, 2.04, 2.04], [-1, -1, -1, 0]))
A_tri = sl.vetores_para_tridiagonal(e, f, g)
x_tri = sl.thomas(e, f, g, [40.8, 0.8, 0.8, 200.8])
verificar("thomas (sistema tridiagonal da APC 3)",
          np.max(np.abs(sl.residuo(A_tri, x_tri, [40.8, 0.8, 0.8, 200.8]))) < 1e-12,
          f"(x = {np.round(x_tri, 3)})")

M = gerador.uniform(-1, 1, (n, n))
A_spd = M @ M.T + n * np.eye(n)                          # simétrica positiva definida
L_ch = sl.cholesky(A_spd)
verificar("cholesky (L L^T = A)", np.max(np.abs(L_ch @ L_ch.T - A_spd)) < 1e-12)
verificar("resolver_cholesky",
          np.max(np.abs(sl.resolver_cholesky(L_ch, A_spd @ x_ref) - x_ref)) < 1e-10)

flops = sl.novo_contador_flops()
sl.gauss_ingenua(np.eye(5) * 4 + 1, np.ones(5), flops=flops)
teoria_md = sum(m * (m + 2) for m in range(1, 5)) + 5 * 6 // 2
teoria_ss = sum(m * (m + 1) for m in range(1, 5)) + 5 * 4 // 2
verificar("contagem de FLOPs = fórmula teórica (n = 5)",
          flops["mult_div"] == teoria_md and flops["soma_sub"] == teoria_ss,
          f"({flops['mult_div']} mult/div, {flops['soma_sub']} soma/sub)")

# precisão reduzida: com 4 algarismos o erro deve ser maior que com float64
x4 = sl.gauss_pivoteamento(A, b, algarismos=4)
verificar("precisão reduzida (4 algarismos) muda o resultado",
          1e-6 < np.max(np.abs(x4 - x_ref)) < 1.0)

# ---------------------------------------------------------------------------
# Sistemas lineares iterativos
# ---------------------------------------------------------------------------
for nome, kwargs in (("jacobi", {}), ("gauss_seidel", {}),
                     ("gauss_seidel (lambda = 1.2)", {"relaxacao": 1.2})):
    solver = sistemas_iterativos.jacobi if nome == "jacobi" else sistemas_iterativos.gauss_seidel
    x, h = solver(A, b, tol=1e-10, criterio="absoluto", **kwargs)
    verificar(nome, np.max(np.abs(x - x_ref)) < 1e-8, f"({len(h)} iterações)")
A_apc3 = [[2, -6, -1], [-3, -1, 7], [-8, 1, -2]]
A_re, b_re, ordem = sistemas_iterativos.reordenar_para_dominancia(A_apc3, [-38, -34, -20])
verificar("reordenar_para_dominancia (APC 3)",
          ordem == [2, 0, 1] and sistemas_iterativos.dominancia_diagonal(A_re))

# ---------------------------------------------------------------------------
# Sistemas não lineares
# ---------------------------------------------------------------------------
def F2(v):
    """x^2 + y^2 = 4 e x y = 1."""
    return np.array([v[0] ** 2 + v[1] ** 2 - 4, v[0] * v[1] - 1])


def J2(v):
    return np.array([[2 * v[0], 2 * v[1]], [v[1], v[0]]])


x, h = sistemas_nao_lineares.newton_sistemas(F2, J2, [2.0, 0.5], tol=1e-12)
verificar("newton_sistemas 2x2 (Cramer)", np.max(np.abs(F2(x))) < 1e-12, f"({len(h)} iterações)")
x, h = sistemas_nao_lineares.newton_sistemas(F2, None, [2.0, 0.5], tol=1e-10)
verificar("newton_sistemas com Jacobiana numérica", np.max(np.abs(F2(x))) < 1e-9)


def F3(v):
    """Sistema 3x3 com solução (1, 2, 3)."""
    x1, x2, x3 = v
    return np.array([x1 ** 2 + x2 - 3, x2 * x3 - 6, x1 + x2 + x3 ** 2 - 12])


x, h = sistemas_nao_lineares.newton_sistemas(F3, None, [1.5, 1.5, 2.5], tol=1e-10)
verificar("newton_sistemas 3x3 (Gauss com pivoteamento)", np.allclose(x, [1, 2, 3], atol=1e-8))

# ---------------------------------------------------------------------------
# Otimização (máximo de 2 sin x - x^2/10 em [0, 4]: x* = 1.4275517787, Chapra 13.1)
# ---------------------------------------------------------------------------
def menos_g(x):
    return -(2 * math.sin(x) - x ** 2 / 10)


X_OTIMO = 1.4275517787645942
x, _ = otimizacao.razao_aurea(menos_g, 0.0, 4.0, tol=1e-8)
verificar("razao_aurea", abs(x - X_OTIMO) < 1e-7)
x, h = otimizacao.interpolacao_quadratica(menos_g, 0.0, 1.0, 4.0, tol=1e-8)
verificar("interpolacao_quadratica", abs(x - X_OTIMO) < 1e-7, f"({len(h)} iterações)")
x, _ = otimizacao.newton_otimizacao(menos_g, 2.5)
verificar("newton_otimizacao (derivadas numéricas)", abs(x - X_OTIMO) < 1e-6)
contada = otimizacao.FuncaoContada(menos_g)
otimizacao.razao_aurea(contada, 0.0, 4.0, tol=1e-6)
verificar("FuncaoContada (razão áurea gasta 1 avaliação por iteração)", contada.n_avaliacoes < 40,
          f"({contada.n_avaliacoes} avaliações)")

# ---------------------------------------------------------------------------
# EDOs: ordem de convergência em y' = -y, y(0) = 1 -> y(1) = e^-1
# ---------------------------------------------------------------------------
for metodo, ordem_teorica in ((edo.euler, 1), (edo.heun, 2), (edo.rk4, 4)):
    erros = []
    for n_passos in (20, 40):
        t, y = metodo(lambda t, y: -y, 0.0, 1.0, 1.0 / n_passos, n_passos)
        erros.append(abs(y[-1] - math.exp(-1)))
    ordem_observada = math.log2(erros[0] / erros[1])
    verificar(f"{metodo.__name__} (ordem ~ {ordem_teorica})",
              abs(ordem_observada - ordem_teorica) < 0.15, f"(ordem observada {ordem_observada:.2f})")
t, Y = edo.rk4(lambda t, y: np.array([y[1], -y[0]]), 0.0, [1.0, 0.0], 0.01, 628)
verificar("rk4 em sistema (oscilador harmônico)", abs(Y[-1, 0] - math.cos(t[-1])) < 1e-8)

# ---------------------------------------------------------------------------
# EDPs
# ---------------------------------------------------------------------------
# Poisson 2D com solução fabricada T = 1 - x^2 - y^2 (lap T = -4, exata na malha)
n_lado = 11
dx = 1.0 / (n_lado - 1)
xs = np.linspace(0, 1, n_lado)
X, Y2 = np.meshgrid(xs, xs, indexing="ij")
T_exato = 1 - X ** 2 - Y2 ** 2
T0 = T_exato.copy()
T0[1:-1, 1:-1] = 0.0                                    # chute ruim no interior
fonte = np.full_like(T0, 4.0 * dx ** 2)                 # q''' dx^2 / k com q'''/k = 4
T_num, erros = edp.gauss_seidel_laplace_2d(T0, tol=1e-12, fonte=fonte, relaxacao=1.5)
verificar("gauss_seidel_laplace_2d (Poisson, solução fabricada)",
          np.max(np.abs(T_num - T_exato)) < 1e-9, f"({len(erros)} varreduras)")
verificar("residuo_laplace_2d", edp.residuo_laplace_2d(T_num, fonte) < 1e-9)

# Regime permanente 1D com geração: T(x) = T_inf + q L/h + q (L^2 - x^2)/(2k)
L, k_mat, alfa, h_conv, T_inf, q = 0.01, 20.0, 1e-5, 500.0, 300.0, 1e7
for beta in (0.0, 0.5):
    x_nos, t_s, T_s = edp.conducao_transiente_1d(
        L, 21, alfa, k_mat, 1.0, 3000, T_inf, ("adiabatico",),
        ("conveccao", h_conv, T_inf), geracao=q, beta=beta, salvar_a_cada=3000)
    T_analitico = T_inf + q * L / h_conv + q * (L ** 2 - x_nos ** 2) / (2 * k_mat)
    verificar(f"conducao_transiente_1d -> permanente (beta = {beta})",
              np.max(np.abs(T_s[-1] - T_analitico)) < 1e-6)

# ---------------------------------------------------------------------------
print()
if falhas:
    print(f"{len(falhas)} TESTE(S) FALHARAM: {falhas}")
    sys.exit(1)
print("TODOS OS TESTES PASSARAM")
