"""
Aula Prática 02 - Métodos de Newton-Raphson e da Secante
=========================================================
Temperatura de equilíbrio T de uma superfície aquecida que perde calor por
convecção e radiação:

    P = h A (T - T_inf) + eps sigma A (T^4 - T_inf^4)
    f(T) = h A (T - T_inf) + eps sigma A (T^4 - T_inf^4) - P = 0
    f'(T) = h A + 4 eps sigma A T^3

Executar (a partir da raiz do repositório):
    python Aulas_Praticas/AP02_Newton_Secante/main.py
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# ---- acesso à biblioteca utils/ (raiz do repositório) ------------------------
PASTA = Path(__file__).resolve().parent
SAIDA = PASTA / "outputs"
sys.path.insert(0, str(PASTA.parents[1]))
from utils.padrao_grafico import configurar_estilo, salvar_figura  # noqa: E402
from utils.raizes import newton_raphson, secante  # noqa: E402
from utils.saida import Relatorio  # noqa: E402

# ---- DADOS DE ENTRADA (enunciado) --------------------------------------------
AREA = 0.01          # área da superfície A [m2]
H_CONV = 10.0        # coeficiente de convecção h [W/(m2 K)]
EMISSIVIDADE = 0.8   # emissividade eps [-]
T_INF = 300.0        # temperatura ambiente [K]
POTENCIA = 20.0      # potência recebida P [W]
SIGMA = 5.67e-8      # constante de Stefan-Boltzmann [W/(m2 K4)]
T0_NEWTON = 350.0    # chute inicial de Newton [K]
T0_SEC, T1_SEC = 400.0, 425.0   # chutes da secante (lados opostos da raiz) [K]
TOL = 1e-6           # critério |T_{k+1} - T_k| < TOL [K]
MAX_ITER = 50


def residuo_termico(T):
    """f(T): potência perdida - potência recebida [W] (zero no equilíbrio)."""
    return (H_CONV * AREA * (T - T_INF)
            + EMISSIVIDADE * SIGMA * AREA * (T ** 4 - T_INF ** 4) - POTENCIA)


def derivada_residuo(T):
    """f'(T) analítica [W/K]: condutância convectiva + radiativa linearizada."""
    return H_CONV * AREA + 4.0 * EMISSIVIDADE * SIGMA * AREA * T ** 3


def main():
    rel = Relatorio()
    configurar_estilo()
    rel.titulo("AP02 - Temperatura de equilíbrio (convecção + radiação)")

    # ---- conferência da Tabela 1 do enunciado --------------------------------
    rel.subtitulo("Conferência da Tabela 1 (f(T) avaliada pela rotina)")
    rel.tabela([{"T [K]": T, "f(T) [W]": residuo_termico(T)}
                for T in (350, 375, 400, 425, 450, 475, 500)], casas=9)

    # ---- Parte 1: Newton-Raphson ---------------------------------------------
    rel.subtitulo("Parte 1 - Newton-Raphson (T0 = 350 K)")
    T_newton, hist_newton = newton_raphson(residuo_termico, derivada_residuo, T0_NEWTON,
                                           tol=TOL, max_iter=MAX_ITER, criterio="passo")
    rel.tabela(hist_newton, colunas=["k", "x_k", "f(x_k)", "x_k+1", "erro"], casas=12)
    rel.escrever(f"\nT de equilíbrio (Newton) = {T_newton:.6f} K em {len(hist_newton)} iterações")

    # ---- Parte 2: secante (f tratada como caixa-preta, sem f') -------------
    rel.subtitulo("Parte 2 - Secante (T0 = 400 K, T1 = 425 K)")
    T_sec, hist_sec = secante(residuo_termico, T0_SEC, T1_SEC,
                              tol=TOL, max_iter=MAX_ITER, criterio="passo")
    rel.tabela(hist_sec, colunas=["k", "x_k", "f(x_k)", "x_k+1", "erro"], casas=12)
    rel.escrever(f"\nT de equilíbrio (secante) = {T_sec:.6f} K em {len(hist_sec)} iterações")

    # ---- comparação --------------------------------------------------------
    rel.subtitulo("Comparação")
    rel.escrever(f"|T_newton - T_secante| = {abs(T_newton - T_sec):.2e} K")
    rel.escrever(f"f(T_newton) = {residuo_termico(T_newton):.2e} W   "
                 f"f(T_secante) = {residuo_termico(T_sec):.2e} W")
    perda_conv = H_CONV * AREA * (T_newton - T_INF)
    perda_rad = EMISSIVIDADE * SIGMA * AREA * (T_newton ** 4 - T_INF ** 4)
    rel.escrever(f"Balanço: convecção = {perda_conv:.4f} W + radiação = {perda_rad:.4f} W "
                 f"= {perda_conv + perda_rad:.4f} W (P = {POTENCIA} W)")

    # ---- gráficos -------------------------------------------------------------
    T_graf = np.linspace(300, 500, 300)
    fig, ax = plt.subplots()
    ax.plot(T_graf, residuo_termico(T_graf), label="f(T)")
    ax.axhline(0, color="k", lw=0.8)
    ax.plot(T_newton, 0, "o", label=f"T = {T_newton:.3f} K")
    ax.set(xlabel="T [K]", ylabel="f(T) [W]", title="Resíduo do balanço térmico")
    ax.legend()
    salvar_figura(fig, SAIDA / "funcao_residuo.png")

    fig, ax = plt.subplots()
    ax.semilogy([h["k"] for h in hist_newton], [h["erro"] for h in hist_newton], "o-",
                label="Newton-Raphson")
    ax.semilogy([h["k"] for h in hist_sec], [h["erro"] for h in hist_sec], "s-", label="Secante")
    ax.axhline(TOL, color="k", ls="--", lw=0.8, label="tolerância")
    ax.set(xlabel="iteração k", ylabel="|T(k+1) - T(k)| [K]", title="Convergência dos métodos abertos")
    ax.legend()
    salvar_figura(fig, SAIDA / "convergencia.png")

    rel.salvar(SAIDA / "resultados.txt")


if __name__ == "__main__":
    main()
