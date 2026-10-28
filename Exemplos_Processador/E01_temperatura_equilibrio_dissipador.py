"""
E01 - Temperatura de equilíbrio de um processador com dissipador  [RAÍZES 1D]
==============================================================================
ENUNCIADO (estilo prova)
    Um processador dissipa P = 95 W. O calor passa da junção (chip) para a base
    do dissipador por R_jc (junção-carcaça) e R_tim (pasta térmica) e sai do
    dissipador para o ar do gabinete por convecção forçada e radiação:

        P = h A (Ts - T_inf) + eps sigma A (Ts^4 - T_inf^4)      (dissipador)
        Tj = Ts + P (R_jc + R_tim)                                (junção)

    a) Encontre Ts por bisseção, falsa posição, Newton, secante e Müller e
       compare o número de iterações.
    b) Calcule a temperatura da junção Tj.
    c) Qual a potência MÁXIMA (TDP) para que Tj não passe de Tj_max = 95 °C?

MÉTODO: zero de função de uma variável (utils/raizes.py).
COMO ADAPTAR: troque a função residuo_dissipador pela equação da prova;
    os métodos só precisam de f(x) (e de f'(x) no Newton).

Executar (a partir da raiz do repositório):
    python Exemplos_Processador/E01_temperatura_equilibrio_dissipador.py
"""

import sys
from pathlib import Path

PASTA = Path(__file__).resolve().parent
sys.path.insert(0, str(PASTA.parent))
from utils.raizes import (bissecao, falsa_posicao, muller, newton_raphson,  # noqa: E402
                          secante)
from utils.saida import Relatorio  # noqa: E402

# ---- DADOS ------------------------------------------------------------------------
POTENCIA = 95.0         # potência do processador [W]
R_JC = 0.20             # resistência junção-carcaça [K/W]
R_TIM = 0.10            # resistência da pasta térmica [K/W]
H_CONV = 40.0           # convecção forçada no dissipador [W/(m2 K)]
AREA = 0.08             # área total das aletas [m2]
EMISSIVIDADE = 0.85     # alumínio anodizado [-]
SIGMA = 5.67e-8         # Stefan-Boltzmann [W/(m2 K4)]
T_INF = 35.0 + 273.15   # ar dentro do gabinete [K]
TJ_MAX = 95.0 + 273.15  # temperatura máxima de junção [K]
TOL = 1e-8


def calor_dissipador(Ts):
    """Calor que o dissipador entrega ao ar [W] quando a base está a Ts [K]."""
    return H_CONV * AREA * (Ts - T_INF) + EMISSIVIDADE * SIGMA * AREA * (Ts ** 4 - T_INF ** 4)


def residuo_dissipador(Ts, potencia=POTENCIA):
    """f(Ts) = calor dissipado - potência gerada (zero no equilíbrio)."""
    return calor_dissipador(Ts) - potencia


def derivada_residuo(Ts):
    """f'(Ts) = h A + 4 eps sigma A Ts^3 [W/K]."""
    return H_CONV * AREA + 4.0 * EMISSIVIDADE * SIGMA * AREA * Ts ** 3


def residuo_tdp(P):
    """Item c: dado P, a junção está em TJ_MAX -> a base está em
    Ts = TJ_MAX - P (R_jc + R_tim). O resíduo é o balanço do dissipador."""
    Ts = TJ_MAX - P * (R_JC + R_TIM)
    return residuo_dissipador(Ts, potencia=P)


def main():
    rel = Relatorio()
    rel.titulo("E01 - Temperatura de equilíbrio do processador (raízes 1D)")

    # ---- a) Ts por cinco métodos ---------------------------------------------------
    # intervalo: em Ts = T_inf o dissipador não troca calor (f < 0); em Ts = 400 K
    # ele dissipa muito mais que 95 W (f > 0) -> existe raiz entre eles
    Ts_l, Ts_u = T_INF, 400.0
    rel.escrever(f"f({Ts_l:.2f}) = {residuo_dissipador(Ts_l):.3f} W,  "
                 f"f({Ts_u:.2f}) = {residuo_dissipador(Ts_u):.3f} W")
    resultados = {
        "bisseção": bissecao(residuo_dissipador, Ts_l, Ts_u, tol=TOL, criterio="funcao"),
        "falsa posição": falsa_posicao(residuo_dissipador, Ts_l, Ts_u, tol=TOL, criterio="funcao"),
        "Newton-Raphson": newton_raphson(residuo_dissipador, derivada_residuo, 350.0, tol=TOL),
        "secante": secante(residuo_dissipador, 340.0, 360.0, tol=TOL),
        "Müller": muller(residuo_dissipador, 320.0, 340.0, 360.0, tol=TOL),
    }
    rel.subtitulo("a) Temperatura da base do dissipador")
    rel.tabela([{"método": nome, "Ts [K]": Ts, "Ts [°C]": Ts - 273.15, "iterações": len(hist)}
                for nome, (Ts, hist) in resultados.items()], casas=10)
    rel.subtitulo("Iterações de Newton (convergência quadrática)")
    rel.tabela(resultados["Newton-Raphson"][1], casas=12)

    # ---- b) junção ---------------------------------------------------------------
    Ts = resultados["Newton-Raphson"][0]
    Tj = Ts + POTENCIA * (R_JC + R_TIM)
    rel.subtitulo("b) Temperatura da junção")
    rel.escrever(f"Tj = Ts + P (R_jc + R_tim) = {Ts - 273.15:.3f} + {POTENCIA * (R_JC + R_TIM):.3f}"
                 f" = {Tj - 273.15:.3f} °C")
    conv = H_CONV * AREA * (Ts - T_INF)
    rel.escrever(f"Divisão do calor: convecção = {conv:.2f} W ({100 * conv / POTENCIA:.1f} %), "
                 f"radiação = {POTENCIA - conv:.2f} W ({100 * (1 - conv / POTENCIA):.1f} %)")

    # ---- c) potência máxima ---------------------------------------------------------
    rel.subtitulo("c) Potência máxima com Tj <= 95 °C")
    P_max, hist = bissecao(residuo_tdp, 1.0, 200.0, tol=1e-10, criterio="funcao")
    Ts_max = TJ_MAX - P_max * (R_JC + R_TIM)
    rel.escrever(f"P_max = {P_max:.4f} W  (Ts = {Ts_max - 273.15:.3f} °C) em {len(hist)} iterações")
    rel.escrever(f"Conferência: f(P_max) = {residuo_tdp(P_max):.2e} W")

    rel.salvar(PASTA / "outputs" / "E01_resultados.txt")


if __name__ == "__main__":
    main()
