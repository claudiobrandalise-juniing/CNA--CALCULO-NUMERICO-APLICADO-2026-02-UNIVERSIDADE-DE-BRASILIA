"""
E07 - Velocidade ótima da ventoinha  [OTIMIZAÇÃO 1D + RAIZ DENTRO DA FUNÇÃO OBJETIVO]
====================================================================================
ENUNCIADO (estilo prova)
    Girar a ventoinha mais rápido gasta mais energia (P_fan = c v^3), mas esfria
    o chip, o que reduz a potência de fuga (que cresce exponencialmente com Tj).
    Para cada velocidade do ar v [m/s]:

        hA(v)    = 1.2 v^0.8 + 0.5                       [W/K]
        R_tot(v) = R_js + 1/hA(v)                        [K/W]
        Tj(v)    : raiz de  g(Tj) = Tj - T_inf - P(Tj) R_tot(v) = 0
        P(Tj)    = P_din + P_fuga0 exp[gamma (Tj - T_ref)]

    Potência elétrica total a minimizar:  E(v) = c v^3 + P(Tj(v)),  1.5 <= v <= 6 m/s.

    a) Encontre o ótimo por razão áurea e por interpolação quadrática.
    b) Se, por durabilidade, Tj não pode passar de 85 °C, qual a velocidade
       mínima permitida (bisseção em Tj(v) - 85 = 0)? Qual é o ótimo restrito?

MÉTODOS: utils/otimizacao.py (razão áurea, interp. quadrática, FuncaoContada)
         utils/raizes.py (Newton para Tj; bisseção para a restrição)
COMO ADAPTAR: a função objetivo pode ser qualquer coisa calculável, inclusive
    algo que exija resolver uma equação por dentro (como aqui).

Executar (a partir da raiz do repositório):
    python Exemplos_Processador/E07_velocidade_otima_ventoinha.py
"""

import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

PASTA = Path(__file__).resolve().parent
sys.path.insert(0, str(PASTA.parent))
from utils.otimizacao import FuncaoContada, interpolacao_quadratica, razao_aurea  # noqa: E402
from utils.padrao_grafico import configurar_estilo, salvar_figura  # noqa: E402
from utils.raizes import bissecao, newton_raphson  # noqa: E402
from utils.saida import Relatorio  # noqa: E402

# ---- DADOS -----------------------------------------------------------------------------
C_FAN = 0.15                 # P_fan = C_FAN v^3 [W/(m/s)^3]
P_DIN, P_FUGA0 = 70.0, 10.0  # [W]
T_REF, GAMMA = 50.0, 0.02    # [°C], [1/K]
R_JS = 0.30                  # [K/W]
T_INF = 35.0                 # [°C]
V_MIN, V_MAX = 1.5, 6.0      # faixa de velocidades [m/s]
TJ_LIMITE = 85.0             # restrição do item b [°C]
TOL = 1e-6


def conveccao(v):
    """hA(v): convecção forçada (~v^0.8) + parcela natural [W/K]."""
    return 1.2 * v ** 0.8 + 0.5


def potencia_chip(Tj):
    """P(Tj) com fuga exponencial [W]."""
    return P_DIN + P_FUGA0 * math.exp(GAMMA * (Tj - T_REF))


def temperatura_juncao(v):
    """Resolve g(Tj) = Tj - T_inf - P(Tj) R_tot = 0 por Newton (problema interno)."""
    R_tot = R_JS + 1.0 / conveccao(v)
    g = lambda Tj: Tj - T_INF - potencia_chip(Tj) * R_tot                        # noqa: E731
    dg = lambda Tj: 1.0 - GAMMA * P_FUGA0 * math.exp(GAMMA * (Tj - T_REF)) * R_tot  # noqa: E731
    Tj, _ = newton_raphson(g, dg, T_INF + P_DIN * R_tot, tol=1e-10, max_iter=50)
    return Tj


def energia_total(v):
    """Função objetivo E(v) = P_fan + P_chip [W]."""
    return C_FAN * v ** 3 + potencia_chip(temperatura_juncao(v))


def main():
    rel = Relatorio()
    configurar_estilo()
    rel.titulo("E07 - Velocidade ótima da ventoinha (otimização 1D)")

    rel.subtitulo("Exploração da função objetivo")
    rel.tabela([{"v [m/s]": v, "Tj [°C]": temperatura_juncao(v), "P_fan [W]": C_FAN * v ** 3,
                 "P_chip [W]": potencia_chip(temperatura_juncao(v)), "E [W]": energia_total(v)}
                for v in (1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0)], casas=7)

    # ---- a) ótimo sem restrição ----------------------------------------------------------
    rel.subtitulo("a) Ótimo sem restrição")
    linhas = []
    for nome in ("razão áurea", "interpolação quadrática"):
        E = FuncaoContada(energia_total)
        if nome == "razão áurea":
            v_otimo, hist = razao_aurea(E, V_MIN, V_MAX, tol=TOL)
        else:
            v_otimo, hist = interpolacao_quadratica(E, 2.0, 3.0, 4.0, tol=TOL)
        linhas.append({"método": nome, "v ótima [m/s]": v_otimo, "E mín [W]": energia_total(v_otimo),
                       "Tj [°C]": temperatura_juncao(v_otimo), "iterações": len(hist),
                       "avaliações de E": E.n_avaliacoes})
    rel.tabela(linhas, casas=8)
    k_x2 = next(h["k"] for h in hist if abs(h["x2"] - v_otimo) < TOL)
    protecoes = sum(1 for h in hist if h["passo"] != "parabólico")
    rel.escrever(f"Interpolação quadrática: o ponto central x2 já estava a menos de {TOL:g} m/s do"
                 f" ótimo na iteração {k_x2}, mas o critério x3 - x1 < tol só foi atingido depois de"
                 f" {protecoes} passos áureos de proteção (os extremos encolhem devagar). Perto do"
                 " mínimo, E(v) varia menos que o ruído do Newton interno e a parábola perde precisão.")

    # ---- b) restrição Tj <= 85 °C -------------------------------------------------------
    rel.subtitulo(f"b) Restrição Tj <= {TJ_LIMITE} °C")
    v_restricao, hist = bissecao(lambda v: temperatura_juncao(v) - TJ_LIMITE, V_MIN, V_MAX,
                                 tol=1e-8, criterio="funcao")
    v_irrestrito = linhas[0]["v ótima [m/s]"]
    v_final = max(v_irrestrito, v_restricao)
    rel.escrever(f"Velocidade mínima para Tj = {TJ_LIMITE} °C: v = {v_restricao:.5f} m/s "
                 f"({len(hist)} iterações de bisseção)")
    rel.escrever(f"O ótimo irrestrito (v = {v_irrestrito:.4f} m/s) daria Tj = "
                 f"{temperatura_juncao(v_irrestrito):.2f} °C > {TJ_LIMITE} °C -> restrição ATIVA")
    rel.escrever(f"Ótimo restrito: v = {v_final:.4f} m/s, E = {energia_total(v_final):.4f} W "
                 f"(+{energia_total(v_final) - energia_total(v_irrestrito):.3f} W em relação ao irrestrito)")
    rel.escrever(f"Abaixo de v ~ {V_MIN} m/s não há regime permanente (fuga térmica, ver E06).")

    v_graf = np.linspace(V_MIN, V_MAX, 200)
    fig, ax = plt.subplots()
    ax.plot(v_graf, [energia_total(v) for v in v_graf], "k", label="E total")
    ax.plot(v_graf, C_FAN * v_graf ** 3, label="P ventoinha")
    ax.plot(v_graf, [potencia_chip(temperatura_juncao(v)) for v in v_graf], label="P chip (com fuga)")
    ax.axvline(v_restricao, color="C3", ls=":", label=f"Tj = {TJ_LIMITE} °C")
    ax.plot(v_irrestrito, energia_total(v_irrestrito), "o", color="k")
    ax.set(xlabel="velocidade do ar v [m/s]", ylabel="potência [W]", ylim=(0, 130),
           title="Compromisso ventoinha x fuga")
    ax.legend()
    salvar_figura(fig, PASTA / "outputs" / "E07_ventoinha.png")

    rel.salvar(PASTA / "outputs" / "E07_resultados.txt")


if __name__ == "__main__":
    main()
