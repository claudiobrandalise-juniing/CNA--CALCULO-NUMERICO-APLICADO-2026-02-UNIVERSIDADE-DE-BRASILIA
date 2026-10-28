"""
PPC 03 - Condução transiente 1D com geração interna (combustível nuclear)
==========================================================================
    dT/dt = alfa d2T/dx2 + q'''/(rho cp),   0 < x < L   (L = meia-espessura)
    x = 0: simetria (adiabático)       x = L: -k dT/dx = h (T - T_inf)

Diferenças finitas + esquema implícito (beta = 0) + algoritmo de Thomas.

    Tarefa 1: validação SEM geração contra a solução exata em série
    Tarefa 2: geração interna; limite q''' -> 0 e regime permanente analítico
    Tarefa 3: visualizações (perfis, mapa de cores, animação)
    Tarefa 4: valores realistas (placa de UO2 num reator PWR)

Executar (a partir da raiz do repositório):
    python PPC_03_Conducao_Transiente/main.py
"""

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

# ---- acesso à biblioteca utils/ (raiz do repositório) ------------------------
PASTA = Path(__file__).resolve().parent
SAIDA = PASTA / "outputs"
sys.path.insert(0, str(PASTA.parent))
sys.path.insert(0, str(PASTA))
from solucao_exata import autovalores, theta_exato  # noqa: E402
from utils.edp import conducao_transiente_1d, matriz_tridiagonal_1d  # noqa: E402
from utils.padrao_grafico import configurar_estilo, salvar_figura  # noqa: E402
from utils.saida import Relatorio  # noqa: E402
from utils.sistemas_lineares import vetores_para_tridiagonal  # noqa: E402

with open(PASTA / "inputs" / "parametros.json", encoding="utf-8") as arquivo:
    PARAM = json.load(arquivo)
VAL = PARAM["validacao_adimensional"]
COMB = PARAM["combustivel"]


def resolver_adimensional(n_nos, dFo, Fo_final, beta=0.0, salvar_a_cada=1):
    """Problema adimensional: L = alfa = k = 1, h = Bi, T_inf = 0, theta_inicial = 1."""
    n_passos = int(round(Fo_final / dFo))
    return conducao_transiente_1d(1.0, n_nos, 1.0, 1.0, dFo, n_passos, 1.0,
                                  ("simetria",), ("conveccao", VAL["Bi"], 0.0),
                                  beta=beta, salvar_a_cada=salvar_a_cada)


def erro_contra_exata(n_nos, dFo, Fo_alvo, lambdas, beta=0.0):
    """Maior |theta_num - theta_exato| no instante Fo_alvo."""
    x, Fo, theta = resolver_adimensional(n_nos, dFo, Fo_alvo, beta, salvar_a_cada=10 ** 9)
    exata = theta_exato(x, Fo_alvo, VAL["Bi"], lambdas=lambdas)
    return np.max(np.abs(theta[-1] - exata))


def tarefa1(rel, lambdas):
    """Validação sem geração + estudo de malha e de passo de tempo."""
    rel.subtitulo("Tarefa 1a - Matriz do esquema implícito para 7 nós (Fo_dt = 0.5, Bi = 1)")
    dx = 1.0 / 6.0
    e, f, g = matriz_tridiagonal_1d(7, 0.5, 0.0, dx, 1.0, ("simetria",), ("conveccao", VAL["Bi"], 0.0))
    rel.matriz(vetores_para_tridiagonal(e, f, g), "[A]", casas=5)
    rel.escrever(f"(1a linha: 1+2Fo e -2Fo; última: -2Fo e 1+2Fo+2Fo Bi_malha, Bi_malha = h dx/k = "
                 f"{VAL['Bi'] * dx:.5f})")

    rel.subtitulo(f"Tarefa 1b - Numérico x exato (Bi = {VAL['Bi']}, {VAL['n_nos_escolhido']} nós, "
                  f"dFo = {VAL['dFo_escolhido']})")
    rel.escrever("Primeiros autovalores lambda_n (bisseção em g = lambda sin - Bi cos): "
                 + ", ".join(f"{v:.6f}" for v in lambdas[:5]))
    x, Fo, theta = resolver_adimensional(VAL["n_nos_escolhido"], VAL["dFo_escolhido"], VAL["Fo_final"])
    fig, ax = plt.subplots()
    linhas = []
    for Fo_alvo in VAL["Fo_comparacao"]:
        indice = int(np.argmin(np.abs(Fo - Fo_alvo)))
        exata = theta_exato(x, Fo_alvo, VAL["Bi"], lambdas=lambdas)
        linhas.append({"Fo": Fo_alvo, "theta(0) num.": theta[indice, 0], "theta(0) exato": exata[0],
                       "erro máx": np.max(np.abs(theta[indice] - exata))})
        ax.plot(x, exata, "k-", lw=1)
        ax.plot(x[::4], theta[indice, ::4], "o", ms=4, label=f"Fo = {Fo_alvo}")
    rel.tabela(linhas, casas=8)
    ax.set(xlabel="x* = x/L", ylabel="θ = (T - T∞)/(Ti - T∞)",
           title=f"Validação sem geração (Bi = {VAL['Bi']}): linhas = exata, pontos = numérico")
    ax.legend()
    salvar_figura(fig, SAIDA / "tarefa1_validacao_perfis.png")

    # ---- estudo de malha (espacial) e de passo (temporal) em Fo = 0.2 ------------------
    rel.subtitulo("Tarefa 1c - Convergência em Fo = 0.2")
    erros_x = [erro_contra_exata(n, 1e-4, 0.2, lambdas, beta=0.5) for n in VAL["n_nos_estudo"]]
    rel.escrever("Malha (Crank-Nicolson com dFo = 1e-4 para isolar o erro espacial):")
    rel.tabela([{"n_nos": n, "dx": 1 / (n - 1), "erro máx": e_} for n, e_ in zip(VAL["n_nos_estudo"], erros_x)])
    erros_t = {beta: [erro_contra_exata(161, d, 0.2, lambdas, beta) for d in VAL["dFo_estudo"]]
               for beta in (0.0, 0.5)}
    rel.escrever("Passo de tempo (161 nós para isolar o erro temporal):")
    rel.tabela([{"dFo": d, "erro implícito": erros_t[0.0][i], "erro Crank-Nicolson": erros_t[0.5][i]}
                for i, d in enumerate(VAL["dFo_estudo"])])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    dxs = [1 / (n - 1) for n in VAL["n_nos_estudo"]]
    ax1.loglog(dxs, erros_x, "o-", label="numérico")
    ax1.loglog(dxs, [erros_x[0] * (d / dxs[0]) ** 2 for d in dxs], "k--", label="O(Δx²)")
    ax1.set(xlabel="Δx*", ylabel="erro máximo", title="Convergência espacial")
    ax1.legend()
    for beta, rotulo in ((0.0, "implícito"), (0.5, "Crank-Nicolson")):
        ax2.loglog(VAL["dFo_estudo"], erros_t[beta], "o-", label=rotulo)
    # retas de referência ancoradas no menor passo (regime assintótico)
    d_fino = VAL["dFo_estudo"][-1]
    ax2.loglog(VAL["dFo_estudo"], [erros_t[0.0][-1] * d / d_fino for d in VAL["dFo_estudo"]], "k--",
               label="O(Δt)")
    ax2.loglog(VAL["dFo_estudo"], [erros_t[0.5][-1] * (d / d_fino) ** 2 for d in VAL["dFo_estudo"]], "k:",
               label="O(Δt²)")
    ax2.set(xlabel="ΔFo", ylabel="erro máximo", title="Convergência temporal")
    ax2.legend()
    for eixo in (ax1, ax2):
        eixo.minorticks_off()                 # evita rótulos sobrepostos no eixo log
    salvar_figura(fig, SAIDA / "tarefa1_convergencia.png")

    # ---- estabilidade: explícito x implícito ----------------------------------------------
    rel.subtitulo("Tarefa 1d - Estabilidade (21 nós, até Fo = 0.2)")
    dx = 1.0 / 20
    linhas = []
    for beta, Fo_dt in ((1.0, 0.4), (1.0, 0.6), (0.0, 0.6), (0.0, 5.0)):
        _, _, theta = resolver_adimensional(21, Fo_dt * dx ** 2, 0.2, beta, salvar_a_cada=10 ** 9)
        linhas.append({"esquema": "explícito" if beta == 1.0 else "implícito", "Fo_dt": Fo_dt,
                       "max|theta| final": np.max(np.abs(theta[-1])),
                       "situação": "estável" if np.max(np.abs(theta[-1])) <= 1.0 else "INSTÁVEL"})
    rel.tabela(linhas, casas=6)


def tarefa2_3(rel):
    """Problema físico com geração: limite q -> 0, permanente e visualizações."""
    L, k, h = COMB["meia_espessura_m"], COMB["k_W_mK"], COMB["h_W_m2K"]
    alfa = k / (COMB["rho_kg_m3"] * COMB["cp_J_kgK"])
    T_inf, q = COMB["T_fluido_K"], COMB["q_geracao_W_m3"]
    dt = VAL["dFo_escolhido"] * L ** 2 / alfa           # mesmo dFo validado na tarefa 1
    n_passos = int(round(COMB["tempo_final_s"] / dt))
    n_nos = VAL["n_nos_escolhido"]
    rel.subtitulo("Tarefa 4 - Parâmetros físicos (placa de UO2, reator PWR)")
    rel.escrever(f"L = {L * 1000:.1f} mm, k = {k} W/(m K), rho = {COMB['rho_kg_m3']} kg/m3, "
                 f"cp = {COMB['cp_J_kgK']} J/(kg K)")
    rel.escrever(f"alfa = {alfa:.4e} m2/s,  Bi = hL/k = {h * L / k:.1f},  L^2/alfa = {L ** 2 / alfa:.2f} s")
    rel.escrever(f"q''' = {q:.2e} W/m3, h = {h} W/(m2 K), T_fluido = {T_inf} K;  "
                 f"dt = {dt:.4e} s, {n_passos} passos, {n_nos} nós")

    def simular(geracao):
        return conducao_transiente_1d(L, n_nos, alfa, k, dt, n_passos, COMB["T_inicial_K"],
                                      ("simetria",), ("conveccao", h, T_inf), geracao=geracao,
                                      salvar_a_cada=max(1, n_passos // 600))

    # ---- limite q -> 0: deve coincidir com o caso sem geração ---------------------------------
    rel.subtitulo("Tarefa 2a - Limite assintótico q''' -> 0 (partindo de T_i = 1000 K)")
    T_teste = 1000.0
    casos = {}
    for geracao in (0.0, COMB["q_limite_teste_W_m3"]):
        casos[geracao] = conducao_transiente_1d(L, n_nos, alfa, k, dt, int(round(20.0 / dt)), T_teste,
                                                ("simetria",), ("conveccao", h, T_inf),
                                                geracao=geracao, salvar_a_cada=10 ** 9)
    diferenca = np.max(np.abs(casos[0.0][2][-1] - casos[COMB["q_limite_teste_W_m3"]][2][-1]))
    rel.escrever(f"max|T(q = {COMB['q_limite_teste_W_m3']:g}) - T(q = 0)| em t = 20 s: {diferenca:.2e} K")

    # ---- simulação com geração e conferência com o regime permanente ----------------------
    rel.subtitulo("Tarefa 2b - Com geração: aproximação do regime permanente analítico")
    x, t, T = simular(q)
    T_permanente = T_inf + q * L / h + q * (L ** 2 - x ** 2) / (2.0 * k)
    rel.escrever(f"Permanente analítico: T(0) = {T_permanente[0]:.2f} K, T(L) = {T_permanente[-1]:.2f} K")
    rel.escrever(f"Numérico em t = {t[-1]:.1f} s: T(0) = {T[-1, 0]:.2f} K, T(L) = {T[-1, -1]:.2f} K, "
                 f"max|T - T_perm| = {np.max(np.abs(T[-1] - T_permanente)):.2e} K")
    instantes = [1, 2, 5, 10, 20, 40, 120]
    rel.tabela([{"t [s]": ts, "T centro [K]": T[int(np.argmin(np.abs(t - ts))), 0],
                 "T superfície [K]": T[int(np.argmin(np.abs(t - ts))), -1]} for ts in instantes], casas=7)

    # ---- tarefa 3: visualizações ------------------------------------------------------------
    fig, ax = plt.subplots()
    for ts in instantes:
        ax.plot(x * 1000, T[int(np.argmin(np.abs(t - ts)))], label=f"t = {ts} s")
    ax.plot(x * 1000, T_permanente, "k--", label="permanente (analítico)")
    ax.set(xlabel="x [mm] (0 = centro)", ylabel="T [K]", title="Perfis de temperatura no combustível")
    ax.legend(fontsize=8)
    salvar_figura(fig, SAIDA / "tarefa3_perfis.png")

    fig, ax = plt.subplots()
    malha = ax.pcolormesh(x * 1000, t, T, shading="auto", cmap="inferno")
    fig.colorbar(malha, label="T [K]")
    ax.set(xlabel="x [mm]", ylabel="t [s]", title="Campo T(x, t): aquecimento após ligar a geração")
    ax.set_ylim(0, 40)
    ax.grid(False)
    salvar_figura(fig, SAIDA / "tarefa3_mapa_x_t.png")

    fig, ax = plt.subplots()
    ax.plot(t, T[:, 0], label="centro (x = 0)")
    ax.plot(t, T[:, -1], label="superfície (x = L)")
    ax.axhline(T_permanente[0], color="k", ls="--", lw=0.8)
    ax.set(xlabel="t [s]", ylabel="T [K]", title="Evolução temporal no centro e na superfície")
    ax.legend()
    salvar_figura(fig, SAIDA / "tarefa3_historico.png")

    # animação GIF (apenas os primeiros 40 s, onde está o transiente)
    quadros = [i for i in range(len(t)) if t[i] <= 40.0][::4]
    fig, ax = plt.subplots(figsize=(6, 4))
    linha, = ax.plot(x * 1000, T[0])
    ax.plot(x * 1000, T_permanente, "k--", lw=0.8)
    ax.set(xlim=(0, L * 1000), ylim=(T_inf - 20, T_permanente[0] + 50), xlabel="x [mm]", ylabel="T [K]")
    titulo = ax.set_title("")

    def atualizar(i):
        linha.set_ydata(T[i])
        titulo.set_text(f"t = {t[i]:5.1f} s")
        return linha, titulo

    FuncAnimation(fig, atualizar, frames=quadros).save(SAIDA / "tarefa3_animacao.gif",
                                                      writer=PillowWriter(fps=12))
    plt.close(fig)
    print(f"[animação salva em {SAIDA / 'tarefa3_animacao.gif'}]")


def main():
    rel = Relatorio()
    configurar_estilo()
    rel.titulo("PPC 03 - Condução transiente 1D com geração (diferenças finitas + Thomas)")
    lambdas = autovalores(VAL["Bi"], VAL["n_termos_serie"])
    tarefa1(rel, lambdas)
    tarefa2_3(rel)
    rel.salvar(SAIDA / "resultados.txt")


if __name__ == "__main__":
    main()
