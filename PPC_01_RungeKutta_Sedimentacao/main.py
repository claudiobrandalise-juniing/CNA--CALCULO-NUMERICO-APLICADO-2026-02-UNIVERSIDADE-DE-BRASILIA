"""
PPC 01 - Runge-Kutta de 4a ordem: sedimentação de uma esfera em baixo Reynolds
==============================================================================
Equação adimensional do movimento (APC 1):

    St dv/dt = 1 - v - (3/8) Re v^2,        v(0) = 0

    Re -> 0 (Stokes):  v(t) = 1 - exp(-t/St)                         (exata)
    Re != 0 (Oseen):   v(t) = v* + [Q/P + (-1/v* - Q/P) exp(-P t)]^-1  (Riccati, exata)
        v* = (-1 + sqrt(1 + 3Re/2)) / (3Re/4)    (velocidade terminal)
        Q  = 3 Re / (8 St)
        P  = -(3/4)(Re/St) v* - 1/St

Tarefas do enunciado:
    1. Re -> 0: numérico x analítico para vários St
    2. refinamento do passo h (ordem de convergência)
    3. Re != 0 resolvido numericamente
    4. validação com a solução de Riccati
    5. comportamento para vários Re e desvio em relação a Re -> 0

Executar (a partir da raiz do repositório):
    python PPC_01_RungeKutta_Sedimentacao/main.py
"""

import json
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# ---- acesso à biblioteca utils/ (raiz do repositório) ------------------------
PASTA = Path(__file__).resolve().parent
SAIDA = PASTA / "outputs"
sys.path.insert(0, str(PASTA.parent))
from utils.edo import rk4  # noqa: E402
from utils.metricas_erro import rmse  # noqa: E402
from utils.padrao_grafico import configurar_estilo, salvar_figura  # noqa: E402
from utils.saida import Relatorio  # noqa: E402

# ---- leitura dos parâmetros (inputs/parametros.json) ----------------------------
with open(PASTA / "inputs" / "parametros.json", encoding="utf-8") as arquivo:
    PARAM = json.load(arquivo)


def aceleracao(St, Re):
    """Devolve f(t, v) = dv/dt = (1 - v - 3/8 Re v^2)/St para o RK4."""
    return lambda t, v: (1.0 - v - 0.375 * Re * v ** 2) / St


def v_stokes(t, St):
    """Solução exata para Re -> 0 (fator integrante): v = 1 - e^(-t/St)."""
    return 1.0 - np.exp(-t / St)


def velocidade_terminal(Re):
    """v* (raiz positiva de 3/8 Re v^2 + v - 1 = 0)."""
    return (-1.0 + math.sqrt(1.0 + 1.5 * Re)) / (0.75 * Re)


def v_oseen(t, St, Re):
    """Solução exata da equação de Riccati (arrasto de Stokes + Oseen), v(0) = 0."""
    v_term = velocidade_terminal(Re)
    Q = 3.0 * Re / (8.0 * St)
    P = -0.75 * (Re / St) * v_term - 1.0 / St
    return v_term + 1.0 / (Q / P + (-1.0 / v_term - Q / P) * np.exp(-P * t))


def integrar(St, Re, h, t_final):
    """Resolve o PVI por RK4 do instante 0 até t_final com passo h."""
    n_passos = int(round(t_final / h))
    return rk4(aceleracao(St, Re), 0.0, PARAM["v0"], h, n_passos)


def main():
    rel = Relatorio()
    configurar_estilo()
    rel.titulo("PPC 01 - RK4: sedimentação de uma esfera (Stokes e Oseen)")
    h, t_final = PARAM["h"], PARAM["t_final"]

    # ---- Tarefa 1: Re -> 0, vários St ---------------------------------------------
    rel.subtitulo(f"Tarefa 1 - Re -> 0: RK4 x solução exata (h = {h})")
    fig, ax = plt.subplots()
    linhas = []
    for St in PARAM["St_lista"]:
        t, v = integrar(St, 0.0, h, t_final)
        exata = v_stokes(t, St)
        linhas.append({"St": St, "RMSE": rmse(v, exata), "erro máximo": np.max(np.abs(v - exata)),
                       "v(t_final)": v[-1]})
        ax.plot(t, v, label=f"St = {St}")
        ax.plot(t[::150], exata[::150], "k.", ms=3)
    rel.tabela(linhas, casas=6)
    ax.set(xlabel="t*", ylabel="v*", xlim=(0, 15),
           title="Re → 0: RK4 (linhas) x exata (pontos)")
    ax.legend()
    salvar_figura(fig, SAIDA / "tarefa1_stokes_varios_St.png")

    # ---- Tarefa 2: refinamento do passo de tempo --------------------------------------
    St = PARAM["St_estudo_passo"]
    rel.subtitulo(f"Tarefa 2 - Refinamento de h (St = {St}, Re -> 0, até t = 10)")
    passos, erros = [], []
    for h_teste in PARAM["h_lista"]:
        t, v = integrar(St, 0.0, h_teste, 10.0)
        passos.append(h_teste)
        erros.append(np.max(np.abs(v - v_stokes(t, St))))
    linhas = [{"h": passos[i], "erro máximo": erros[i],
               "ordem observada": (math.log(erros[i - 1] / erros[i]) / math.log(passos[i - 1] / passos[i])
                                   if i > 0 else "-")} for i in range(len(passos))]
    rel.tabela(linhas, casas=6)
    fig, ax = plt.subplots()
    ax.loglog(passos, erros, "o-", label="RK4")
    referencia = [erros[0] * (p / passos[0]) ** 4 for p in passos]
    ax.loglog(passos, referencia, "k--", label="inclinação 4 (O(h⁴))")
    ax.set(xlabel="passo h", ylabel="erro máximo", title="Convergência do RK4 com o refinamento de h")
    ax.legend()
    ax.minorticks_off()                       # evita rótulos sobrepostos no eixo log
    salvar_figura(fig, SAIDA / "tarefa2_convergencia_h.png")

    # ---- Tarefas 3 e 4: Re != 0, validação com Riccati ------------------------------------
    St = PARAM["St_inercial"]
    rel.subtitulo(f"Tarefas 3 e 4 - Re != 0: RK4 x solução de Riccati (St = {St}, h = {h})")
    linhas, curvas = [], {}
    for Re in PARAM["Re_lista"]:
        t, v = integrar(St, Re, h, t_final)
        exata = v_oseen(t, St, Re)
        curvas[Re] = (t, v)
        linhas.append({"Re": Re, "v* terminal (exata)": velocidade_terminal(Re),
                       "v(t_final) RK4": v[-1], "RMSE": rmse(v, exata),
                       "erro máximo": np.max(np.abs(v - exata))})
    rel.tabela(linhas, casas=8)

    # ---- Tarefa 5: desvio em relação a Re -> 0 ------------------------------------------------
    rel.subtitulo("Tarefa 5 - Desvio em relação ao regime de Stokes")
    t0, v0 = integrar(St, 0.0, h, t_final)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    ax1.plot(t0, v0, "k", label="Re → 0 (Stokes)")
    linhas = []
    for Re, (t, v) in curvas.items():
        ax1.plot(t, v, label=f"Re = {Re}")
        ax2.plot(t, v0 - v, label=f"Re = {Re}")
        linhas.append({"Re": Re, "v* terminal": velocidade_terminal(Re),
                       "redução da v terminal [%]": 100 * (1 - velocidade_terminal(Re))})
    rel.tabela(linhas, casas=6)
    ax1.set(xlabel="t*", ylabel="v*", xlim=(0, 10), title="Velocidade da esfera")
    ax2.set(xlabel="t*", ylabel="v_Stokes - v", xlim=(0, 10), title="Desvio em relação a Re → 0")
    ax1.legend()
    ax2.legend()
    salvar_figura(fig, SAIDA / "tarefa5_efeito_Re.png")

    rel.salvar(SAIDA / "resultados.txt")


if __name__ == "__main__":
    main()
