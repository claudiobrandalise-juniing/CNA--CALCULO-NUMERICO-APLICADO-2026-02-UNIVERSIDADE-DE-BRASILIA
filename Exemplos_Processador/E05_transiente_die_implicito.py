"""
E05 - Pico de temperatura no die durante um "turbo"  [EDP 1D TRANSIENTE: IMPLÍCITO + THOMAS]
============================================================================================
ENUNCIADO (estilo prova)
    O die de silício (espessura L = 0,8 mm, área 12 mm x 12 mm) recebe o calor
    gerado na face ativa (x = 0) como um fluxo q'' = P/A. A outra face (x = L)
    está ligada pela pasta térmica à base do dissipador, que fica a 60 °C em
    escalas de tempo curtas (o dissipador é pesado e quase não aquece em < 1 s);
    a pasta é modelada por um coeficiente equivalente h_eq.

        dT/dt = alfa d2T/dx2
        x = 0: -k dT/dx = q''(t)        x = L: -k dT/dx = h_eq (T - T_base)

    Partindo do regime permanente com 30 W, o processador entra em "turbo"
    (95 W) por 0,1 s e volta a 30 W.
    a) Valide o código: com 95 W constantes, o regime permanente tem que ser
       T(x) = T_base + q''/h + q''(L - x)/k.
    b) Simule o turbo: temperatura máxima na face ativa e tempo de resposta.
    c) Compare o custo: explícito (precisa de dt ~ 1e-6 s) x implícito (dt = 2.5e-4 s).

MÉTODO: utils/edp.py -> conducao_transiente_1d (contorno "fluxo" e "conveccao")
COMO ADAPTAR: troque material, espessura, condições de contorno e potência;
    para potência variável, simule por TRECHOS usando o campo final de um
    trecho como condição inicial do próximo.

Executar (a partir da raiz do repositório):
    python Exemplos_Processador/E05_transiente_die_implicito.py
"""

import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

PASTA = Path(__file__).resolve().parent
sys.path.insert(0, str(PASTA.parent))
from utils.edp import conducao_transiente_1d  # noqa: E402
from utils.padrao_grafico import configurar_estilo, salvar_figura  # noqa: E402
from utils.saida import Relatorio  # noqa: E402

# ---- DADOS (silício a ~80 °C) ----------------------------------------------------------------
L_DIE = 0.8e-3              # espessura [m]
AREA_DIE = 0.012 * 0.012    # área [m2]
K_SI, RHO_SI, CP_SI = 130.0, 2330.0, 705.0   # [W/(m K)], [kg/m3], [J/(kg K)]
ALFA = K_SI / (RHO_SI * CP_SI)               # difusividade [m2/s]
H_EQ = 50000.0              # pasta térmica equivalente [W/(m2 K)]
T_BASE = 60.0               # base do dissipador [°C]
P_NORMAL, P_TURBO = 30.0, 95.0               # [W]
N_NOS = 41
DT = 2.5e-4                 # passo do implícito [s]


def fluxo(potencia):
    """q'' = P / A [W/m2]."""
    return potencia / AREA_DIE


def permanente(x, potencia):
    """Regime permanente analítico (perfil linear): T = T_base + q''/h + q''(L - x)/k."""
    q = fluxo(potencia)
    return T_BASE + q / H_EQ + q * (L_DIE - x) / K_SI


def simular(potencia, T_inicial, duracao, dt=DT, beta=0.0):
    """Um trecho de potência constante (implícito por padrão)."""
    return conducao_transiente_1d(L_DIE, N_NOS, ALFA, K_SI, dt, int(round(duracao / dt)), T_inicial,
                                  ("fluxo", fluxo(potencia)), ("conveccao", H_EQ, T_BASE), beta=beta)


def main():
    rel = Relatorio()
    configurar_estilo()
    rel.titulo("E05 - Pico de temperatura no die (condução transiente 1D implícita)")
    dx = L_DIE / (N_NOS - 1)
    tau = RHO_SI * CP_SI * L_DIE / H_EQ
    rel.escrever(f"alfa = {ALFA:.3e} m2/s, Bi = h L/k = {H_EQ * L_DIE / K_SI:.3f}, "
                 f"tau ~ rho cp L/h = {tau * 1000:.1f} ms, L^2/alfa = {L_DIE ** 2 / ALFA * 1000:.2f} ms")
    rel.escrever(f"Fo da malha com dt = {DT:g} s: {ALFA * DT / dx ** 2:.1f}  (o implícito aceita Fo >> 1/2)")

    # ---- a) validação pelo regime permanente --------------------------------------------
    rel.subtitulo("a) Validação: 95 W constantes por 0,4 s")
    x, t, T = simular(P_TURBO, T_BASE, 0.4)
    T_analitico = permanente(x, P_TURBO)
    rel.escrever(f"Analítico: T(0) = {T_analitico[0]:.4f} °C, T(L) = {T_analitico[-1]:.4f} °C")
    rel.escrever(f"Numérico : T(0) = {T[-1, 0]:.4f} °C, T(L) = {T[-1, -1]:.4f} °C;  "
                 f"max|erro| = {np.max(np.abs(T[-1] - T_analitico)):.2e} °C")

    # ---- b) turbo de 0,1 s (simulação por trechos) -----------------------------------------
    rel.subtitulo("b) Turbo: 30 W -> 95 W (0,1 s) -> 30 W (0,15 s)")
    T0 = permanente(x, P_NORMAL)                       # começa no permanente de 30 W
    _, t1, T1 = simular(P_TURBO, T0, 0.10)
    _, t2, T2 = simular(P_NORMAL, T1[-1], 0.15)        # campo final vira condição inicial
    t_total = np.concatenate([t1, t1[-1] + t2[1:]])
    T_total = np.vstack([T1, T2[1:]])
    face = T_total[:, 0]
    rel.escrever(f"Face ativa: {T0[0]:.2f} °C (30 W) -> máx {face.max():.2f} °C no fim do turbo "
                 f"(permanente de 95 W seria {permanente(0.0, P_TURBO):.2f} °C)")
    alvo = T0[0] + 0.632 * (face.max() - T0[0])
    rel.escrever(f"Tempo para 63,2 % da subida: {t_total[np.argmax(face >= alvo)] * 1000:.1f} ms "
                 f"(estimativa tau = {tau * 1000:.1f} ms)")
    rel.escrever(f"Diferença face ativa - face da pasta no pico: {T_total[np.argmax(face), 0] - T_total[np.argmax(face), -1]:.2f} K"
                 " (condução no silício) x "
                 f"{T_total[np.argmax(face), -1] - T_BASE:.2f} K (pasta)")

    # ---- c) custo: explícito x implícito -------------------------------------------------
    rel.subtitulo("c) Explícito x implícito (primeiros 10 ms do turbo)")
    dt_limite = dx ** 2 / (2 * ALFA * (1 + H_EQ * dx / K_SI))   # limite com o nó de convecção
    linhas = []
    for nome, dt, beta in (("explícito", 0.9 * dt_limite, 1.0), ("explícito", 1.2 * dt_limite, 1.0),
                           ("implícito", DT, 0.0)):
        inicio = time.perf_counter()
        with np.errstate(all="ignore"):              # o caso instável estoura (overflow)
            _, _, T_c = simular(P_TURBO, T0, 0.01, dt=dt, beta=beta)
        estavel = np.all(np.isfinite(T_c[-1]))
        linhas.append({"esquema": nome, "dt [s]": dt, "passos": int(round(0.01 / dt)),
                       "T face em 10 ms [°C]": T_c[-1, 0] if estavel else "INSTÁVEL",
                       "tempo de CPU [s]": time.perf_counter() - inicio})
    rel.tabela(linhas, casas=6)
    rel.escrever(f"Limite do explícito: dt <= dx^2/(2 alfa (1 + h dx/k)) = {dt_limite:.2e} s")

    # ---- gráfico ------------------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    ax1.plot(t_total * 1000, face, label="face ativa (x = 0)")
    ax1.plot(t_total * 1000, T_total[:, -1], label="face da pasta (x = L)")
    ax1.axvspan(0, 100, alpha=0.1, color="C3", label="turbo 95 W")
    ax1.set(xlabel="t [ms]", ylabel="T [°C]", title="Resposta do die ao turbo")
    ax1.legend()
    for instante in (0.0, 0.005, 0.02, 0.1):
        i = int(np.argmin(np.abs(t_total - instante)))
        ax2.plot(x * 1000, T_total[i], label=f"t = {instante * 1000:.0f} ms")
    ax2.set(xlabel="x [mm] (0 = face ativa)", ylabel="T [°C]", title="Perfis no die")
    ax2.legend()
    salvar_figura(fig, PASTA / "outputs" / "E05_turbo_die.png")

    rel.salvar(PASTA / "outputs" / "E05_resultados.txt")


if __name__ == "__main__":
    main()
