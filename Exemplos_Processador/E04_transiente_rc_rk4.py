"""
E04 - Aquecimento transiente do processador: modelo RC de 2 nós  [EDO / RK4]
=============================================================================
ENUNCIADO (estilo prova)
    Modelo de capacitâncias concentradas (cada bloco com T uniforme):
        nó 1 = chip + IHS (C1),  nó 2 = dissipador (C2)

        C1 dT1/dt = P(t) - (T1 - T2)/R12
        C2 dT2/dt = (T1 - T2)/R12 - (T2 - T_inf)/R2a

    a) Valide o RK4 com potência constante (P = 95 W) contra a solução analítica
       (autovalores do sistema 2x2 pela fórmula de Bhaskara).
    b) Simule uma carga real: 15 W (ocioso) -> 95 W em t = 20 s -> 15 W em t = 400 s,
       sem e com "thermal throttling" (se T1 >= 95 °C a potência cai para 65 W).
    c) Mostre o efeito do passo h: ordem 4 e limite de estabilidade do RK4.

MÉTODO: utils/edo.py -> rk4 (sistema de EDOs)
COMO ADAPTAR: escreva a função derivadas(t, T) com um balanço por nó
    (C_i dT_i/dt = calor que entra - calor que sai) e chame rk4.

Executar (a partir da raiz do repositório):
    python Exemplos_Processador/E04_transiente_rc_rk4.py
"""

import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

PASTA = Path(__file__).resolve().parent
sys.path.insert(0, str(PASTA.parent))
from utils.edo import rk4  # noqa: E402
from utils.padrao_grafico import configurar_estilo, salvar_figura  # noqa: E402
from utils.saida import Relatorio  # noqa: E402
from utils.sistemas_lineares import cramer_2x2  # noqa: E402

# ---- DADOS ------------------------------------------------------------------------------
C1, C2 = 17.0, 360.0        # capacitâncias térmicas [J/K] (IHS de cobre; dissipador de alumínio)
R12, R2A = 0.30, 0.35       # resistências [K/W] (chip->dissipador; dissipador->ar)
T_INF = 40.0                # ar do gabinete [°C]
P_OCIOSO, P_CARGA, P_LIMITADA = 15.0, 95.0, 65.0   # potências [W]
T_LIMITE = 95.0             # temperatura que dispara o throttling [°C]
H = 0.5                     # passo de tempo [s]
T_FINAL = 900.0             # duração da simulação [s]


def potencia_programada(t, T1, throttling):
    """Perfil de carga: ocioso, carga plena entre 20 s e 400 s, ocioso de novo."""
    potencia = P_CARGA if 20.0 <= t < 400.0 else P_OCIOSO
    if throttling and potencia == P_CARGA and T1 >= T_LIMITE:
        potencia = P_LIMITADA               # o processador reduz o clock
    return potencia


def derivadas(t, T, potencia):
    """Balanço de energia de cada nó: dT/dt = (entra - sai)/C."""
    T1, T2 = T
    P = potencia(t, T1)
    fluxo_12 = (T1 - T2) / R12              # do chip para o dissipador [W]
    fluxo_2a = (T2 - T_INF) / R2A           # do dissipador para o ar [W]
    return np.array([(P - fluxo_12) / C1, (fluxo_12 - fluxo_2a) / C2])


def solucao_analitica(t, P, T0):
    """T(t) = T_perm + c1 v1 e^(l1 t) + c2 v2 e^(l2 t) para P constante.

    dT/dt = A T + B, com A = [[a11, a12], [a21, a22]]. Os autovalores são as
    raízes de l^2 - tr(A) l + det(A) = 0 (Bhaskara) e v = [a12, l - a11].
    """
    a11, a12 = -1 / (C1 * R12), 1 / (C1 * R12)
    a21, a22 = 1 / (C2 * R12), -1 / (C2 * R12) - 1 / (C2 * R2A)
    B = np.array([P / C1, T_INF / (C2 * R2A)])
    T_perm = np.array(cramer_2x2(a11, a12, a21, a22, -B[0], -B[1]))   # A T = -B
    traco, determinante = a11 + a22, a11 * a22 - a12 * a21
    discriminante = math.sqrt(traco ** 2 - 4 * determinante)
    l1, l2 = (traco + discriminante) / 2, (traco - discriminante) / 2
    v1, v2 = np.array([a12, l1 - a11]), np.array([a12, l2 - a11])
    c1, c2 = cramer_2x2(v1[0], v2[0], v1[1], v2[1], *(np.asarray(T0) - T_perm))
    T = T_perm[None, :] + c1 * np.outer(np.exp(l1 * t), v1) + c2 * np.outer(np.exp(l2 * t), v2)
    return T, (l1, l2), T_perm


def main():
    rel = Relatorio()
    configurar_estilo()
    rel.titulo("E04 - Transiente do processador: modelo RC de 2 nós (RK4)")

    # ---- a) validação com potência constante -------------------------------------------
    rel.subtitulo("a) Validação: P = 95 W constante, partindo de T = T_inf")
    n_passos = int(T_FINAL / H)
    t, T_num = rk4(lambda t, T: derivadas(t, T, lambda t_, T1: P_CARGA), 0.0, [T_INF, T_INF], H, n_passos)
    T_exata, (l1, l2), T_perm = solucao_analitica(t, P_CARGA, [T_INF, T_INF])
    rel.escrever(f"Autovalores: l1 = {l1:.6f} 1/s (tau = {-1 / l1:.1f} s), "
                 f"l2 = {l2:.6f} 1/s (tau = {-1 / l2:.2f} s)")
    rel.escrever(f"Regime permanente: T1 = {T_perm[0]:.3f} °C, T2 = {T_perm[1]:.3f} °C")
    rel.escrever(f"Erro máximo RK4 x analítica (h = {H} s): {np.max(np.abs(T_num - T_exata)):.2e} °C")

    # ---- c) passo de tempo ---------------------------------------------------------------
    rel.subtitulo("c) Efeito do passo h (até t = 300 s)")
    linhas = []
    for h in (8.0, 4.0, 2.0, 1.0, 13.0, 14.0):
        t_h, T_h = rk4(lambda t, T: derivadas(t, T, lambda t_, T1: P_CARGA), 0.0, [T_INF, T_INF],
                       h, int(round(300.0 / h)))
        erro = np.max(np.abs(T_h - solucao_analitica(t_h, P_CARGA, [T_INF, T_INF])[0]))
        linhas.append({"h [s]": h, "h |l_rápido|": h * abs(l2), "erro máx [°C]": erro,
                       "situação": "estável" if erro < 100 else "INSTÁVEL"})
    rel.tabela(linhas, casas=4)
    rel.escrever("Limite de estabilidade do RK4 (autovalor real): h |l| < 2.785 -> "
                 f"h < {2.785 / abs(l2):.1f} s. Quem limita o passo é o nó RÁPIDO (chip).")

    # ---- b) carga real, sem e com throttling --------------------------------------------
    rel.subtitulo("b) Perfil de carga 15 W -> 95 W -> 15 W")
    curvas = {}
    for throttling in (False, True):
        potencia = lambda t, T1, th=throttling: potencia_programada(t, T1, th)  # noqa: E731
        T0 = [T_INF + P_OCIOSO * (R12 + R2A), T_INF + P_OCIOSO * R2A]     # equilíbrio ocioso
        t_b, T_b = rk4(lambda t, T: derivadas(t, T, potencia), 0.0, T0, H, n_passos)
        curvas[throttling] = (t_b, T_b)
        acima = t_b[T_b[:, 0] >= T_LIMITE - 0.01]
        quando = f"atinge {T_LIMITE} °C em t = {acima[0]:.1f} s" if acima.size else \
            f"nunca atinge {T_LIMITE} °C"
        rel.escrever(f"throttling {'LIGADO   ' if throttling else 'DESLIGADO'}: T1 máx = "
                     f"{T_b[:, 0].max():.2f} °C; {quando}")

    fig, ax = plt.subplots()
    for throttling, estilo in ((False, "-"), (True, "--")):
        t_b, T_b = curvas[throttling]
        rotulo = "com throttling" if throttling else "sem throttling"
        ax.plot(t_b, T_b[:, 0], estilo, color="C3", label=f"chip ({rotulo})")
        ax.plot(t_b, T_b[:, 1], estilo, color="C0", label=f"dissipador ({rotulo})")
    ax.axhline(T_LIMITE, color="k", lw=0.8, ls=":", label="limite 95 °C")
    ax.axvspan(20, 400, alpha=0.08, color="C1")
    ax.set(xlabel="t [s]", ylabel="T [°C]", title="Carga plena entre 20 s e 400 s (faixa laranja)")
    ax.legend(fontsize=8)
    salvar_figura(fig, PASTA / "outputs" / "E04_transiente.png")

    rel.salvar(PASTA / "outputs" / "E04_resultados.txt")


if __name__ == "__main__":
    main()
