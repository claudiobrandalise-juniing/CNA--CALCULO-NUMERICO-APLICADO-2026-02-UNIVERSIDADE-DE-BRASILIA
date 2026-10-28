"""
E06 - Fuga térmica (thermal runaway) por corrente de fuga  [NEWTON PARA SISTEMAS]
=================================================================================
ENUNCIADO (estilo prova)
    A potência do chip tem uma parte dinâmica fixa e uma parte de FUGA que cresce
    exponencialmente com a temperatura da junção:
        P(Tj) = P_din + P_fuga0 exp[gamma (Tj - T_ref)]
    O calor vai da junção ao dissipador (R_js) e do dissipador ao ar por
    convecção (hA, depende da ventoinha) e radiação. Em regime permanente:

        F1 = P(Tj) - (Tj - Ts)/R_js                                    = 0
        F2 = (Tj - Ts)/R_js - hA (Ts - T_inf) - eps sigma A_r (Ts^4 - T_inf^4) = 0

    a) Com hA = 3 W/K, encontre (Tj, Ts) por Newton (Jacobiana analítica).
    b) Partindo de um chute alto, mostre que existe uma SEGUNDA solução e
       classifique as duas (estável/instável) pelo sinal de det(J).
    c) Reduzindo a ventoinha (hA menor), as duas soluções se aproximam e somem:
       é a fuga térmica. Encontre o hA crítico resolvendo o sistema 3x3
           F1 = 0,  F2 = 0,  det J(Tj, Ts) = 0     (incógnitas Tj, Ts, hA)
       por Newton com Jacobiana NUMÉRICA (Gauss com pivoteamento por dentro).

MÉTODO: utils/sistemas_nao_lineares.py -> newton_sistemas
COMO ADAPTAR: escreva F(x) (lista de resíduos) e, se der tempo, J(x);
    senão passe J = None e o código usa a Jacobiana numérica.

Executar (a partir da raiz do repositório):
    python Exemplos_Processador/E06_fuga_termica_newton_sistemas.py
"""

import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

PASTA = Path(__file__).resolve().parent
sys.path.insert(0, str(PASTA.parent))
from utils.padrao_grafico import configurar_estilo, salvar_figura  # noqa: E402
from utils.saida import Relatorio  # noqa: E402
from utils.sistemas_nao_lineares import newton_sistemas  # noqa: E402

# ---- DADOS (temperaturas em kelvin) --------------------------------------------------------
P_DIN = 70.0               # potência dinâmica [W]
P_FUGA0 = 10.0             # fuga na temperatura de referência [W]
T_REF = 50.0 + 273.15      # [K]
GAMMA = 0.02               # sensibilidade da fuga [1/K] (dobra a cada ~35 K)
R_JS = 0.30                # junção -> dissipador [K/W]
EPS_AREA = 0.85 * 0.08     # emissividade x área radiante [m2]
SIGMA = 5.67e-8
T_INF = 35.0 + 273.15      # ar do gabinete [K]
HA_NOMINAL = 3.0           # condutância convectiva com a ventoinha nominal [W/K]
KELVIN = 273.15


def potencia(Tj):
    """P(Tj) = P_din + P_fuga0 exp(gamma (Tj - T_ref)) [W]."""
    return P_DIN + P_FUGA0 * math.exp(GAMMA * (Tj - T_REF))


def residuos(x, hA):
    """[F1, F2] para x = [Tj, Ts]."""
    Tj, Ts = x[0], x[1]
    fluxo_js = (Tj - Ts) / R_JS
    return np.array([potencia(Tj) - fluxo_js,
                     fluxo_js - hA * (Ts - T_INF) - EPS_AREA * SIGMA * (Ts ** 4 - T_INF ** 4)])


def jacobiana(x, hA):
    """J = [[dF1/dTj, dF1/dTs], [dF2/dTj, dF2/dTs]] (analítica)."""
    Tj, Ts = x[0], x[1]
    return np.array([
        [GAMMA * P_FUGA0 * math.exp(GAMMA * (Tj - T_REF)) - 1 / R_JS, 1 / R_JS],
        [1 / R_JS, -1 / R_JS - hA - 4 * EPS_AREA * SIGMA * Ts ** 3]])


def det_jacobiana(x, hA):
    """det J (2x2, calculado à mão). det > 0: equilíbrio estável; det < 0: instável."""
    J = jacobiana(x, hA)
    return J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0]


def resolver(hA, chute):
    """Newton 2x2 (Cramer) para um valor de hA."""
    return newton_sistemas(lambda x: residuos(x, hA), lambda x: jacobiana(x, hA), chute,
                           tol=1e-9, max_iter=60)


def main():
    rel = Relatorio()
    configurar_estilo()
    rel.titulo("E06 - Fuga térmica: Newton-Raphson para sistemas")

    # ---- a) e b) duas soluções para hA nominal ---------------------------------------------
    rel.subtitulo(f"a) e b) Equilíbrios com hA = {HA_NOMINAL} W/K")
    solucoes = []
    for chute in ([350.0, 330.0], [500.0, 450.0]):
        x, hist = resolver(HA_NOMINAL, chute)
        rel.escrever(f"Chute (Tj, Ts) = {chute} K:")
        nomes = {"x1": "Tj [K]", "x2": "Ts [K]", "dx1": "dTj", "dx2": "dTs"}
        rel.tabela([{nomes.get(c, c): v for c, v in linha.items()} for linha in hist], casas=10)
        det = det_jacobiana(x, HA_NOMINAL)
        solucoes.append(x)
        rel.escrever(f"  -> Tj = {x[0] - KELVIN:.3f} °C, Ts = {x[1] - KELVIN:.3f} °C, "
                     f"P = {potencia(x[0]):.2f} W, det J = {det:.3f} -> "
                     f"{'ESTÁVEL' if det > 0 else 'INSTÁVEL'}\n")
    rel.escrever("Interpretação: na solução instável, um pequeno aumento de Tj aumenta a fuga MAIS")
    rel.escrever("do que o sistema consegue dissipar -> a temperatura dispara.")

    # ---- c) varredura de hA por continuação + ponto crítico -------------------------------
    rel.subtitulo("c) Redução da ventoinha (continuação) e hA crítico")
    ramos = {"estável": [], "instável": []}
    for nome, x in zip(ramos, solucoes):
        for hA in np.arange(HA_NOMINAL, 1.0, -0.01):
            try:
                x_novo, hist = resolver(hA, x)            # chute = solução do hA anterior
            except ArithmeticError:
                break
            if hist[-1]["erro"] > 1e-6:                   # Newton não convergiu: acabou o ramo
                break
            x = x_novo
            ramos[nome].append((hA, x[0] - KELVIN))
    for nome, pontos in ramos.items():
        rel.escrever(f"Ramo {nome}: existe até hA = {pontos[-1][0]:.2f} W/K "
                     f"(Tj = {pontos[-1][1]:.1f} °C)")

    def sistema_critico(w):
        """Incógnitas w = [Tj, Ts, hA]: equilíbrio E det J = 0 (ponto de dobra)."""
        return np.array([*residuos(w[:2], w[2]), det_jacobiana(w[:2], w[2])])

    chute = [(ramos["estável"][-1][1] + ramos["instável"][-1][1]) / 2 + KELVIN, 370.0,
             ramos["estável"][-1][0]]
    w, hist = newton_sistemas(sistema_critico, None, chute, tol=1e-10)
    rel.escrever(f"\nPonto crítico (Newton 3x3, Jacobiana numérica, {len(hist)} iterações):")
    rel.escrever(f"  hA_crítico = {w[2]:.4f} W/K, Tj = {w[0] - KELVIN:.2f} °C, Ts = {w[1] - KELVIN:.2f} °C")
    rel.escrever(f"  Abaixo de hA = {w[2]:.3f} W/K não existe regime permanente: fuga térmica.")

    fig, ax = plt.subplots()
    for nome, estilo in (("estável", "-"), ("instável", "--")):
        ax.plot([p[0] for p in ramos[nome]], [p[1] for p in ramos[nome]], estilo, label=f"ramo {nome}")
    ax.plot(w[2], w[0] - KELVIN, "ko", label=f"ponto crítico (hA = {w[2]:.3f} W/K)")
    ax.set(xlabel="hA da ventoinha [W/K]", ylabel="Tj [°C]",
           title="Equilíbrios do chip com fuga exponencial")
    ax.invert_xaxis()
    ax.legend()
    salvar_figura(fig, PASTA / "outputs" / "E06_fuga_termica.png")

    rel.salvar(PASTA / "outputs" / "E06_resultados.txt")


if __name__ == "__main__":
    main()
