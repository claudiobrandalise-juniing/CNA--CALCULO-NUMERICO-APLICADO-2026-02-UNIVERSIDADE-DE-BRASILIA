"""
E08 - Constantes de tempo térmicas do processador  [RAÍZES DE POLINÔMIOS: BAIRSTOW / MÜLLER]
===========================================================================================
ENUNCIADO (estilo prova)
    Rede RC de 3 nós: die (1), IHS (2) e dissipador (3), com capacitâncias C_i e
    condutâncias G entre eles e para o ar. Sem potência, a temperatura relativa
    theta = T - T_inf obedece a

        C d(theta)/dt + G theta = 0      (mesma estrutura do PPC 2, sem a massa)

    Propondo theta = phi e^(lambda t):  (lambda C + G) phi = 0  ->  det(lambda C + G) = 0,
    um polinômio de grau 3 em lambda. As constantes de tempo são tau_i = -1/lambda_i.

    a) Monte o polinômio característico e ache as raízes por Bairstow.
    b) Confira com Müller (três pontos iniciais diferentes) e pelas relações de Vieta.
    c) Valide com uma simulação RK4 do resfriamento: no fim, theta decai como
       e^(lambda_lento t).
    d) Comente a rigidez (stiffness) do problema.

MÉTODOS: utils/polinomios.py (determinante_polinomial, bairstow), utils/raizes.py
         (muller), utils/edo.py (rk4)
COMO ADAPTAR: mude C e G; para n nós o polinômio tem grau n.

Executar (a partir da raiz do repositório):
    python Exemplos_Processador/E08_constantes_tempo_bairstow.py
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
from utils.polinomios import avaliar_polinomio, bairstow, determinante_polinomial  # noqa: E402
from utils.raizes import muller  # noqa: E402
from utils.saida import Relatorio  # noqa: E402

# ---- DADOS ------------------------------------------------------------------------------
CAPACITANCIAS = [0.5, 17.0, 360.0]           # die, IHS, dissipador [J/K]
G_DIE_IHS, G_IHS_DISS, G_DISS_AR = 20.0, 10.0, 1 / 0.35   # condutâncias [W/K]
THETA_INICIAL = [50.0, 40.0, 30.0]           # T - T_inf no início do resfriamento [K]
H, T_FINAL = 0.02, 900.0                     # passo e duração do RK4 [s]


def matrizes():
    """C (diagonal) e G (simétrica) da rede de 3 nós."""
    C = np.diag(CAPACITANCIAS)
    G = np.array([[G_DIE_IHS, -G_DIE_IHS, 0.0],
                  [-G_DIE_IHS, G_DIE_IHS + G_IHS_DISS, -G_IHS_DISS],
                  [0.0, -G_IHS_DISS, G_IHS_DISS + G_DISS_AR]])
    return C, G


def main():
    rel = Relatorio()
    configurar_estilo()
    rel.titulo("E08 - Constantes de tempo da rede RC (polinômio característico)")
    C, G = matrizes()
    n = len(CAPACITANCIAS)

    # ---- a) polinômio e Bairstow ----------------------------------------------------------
    # cada elemento de (lambda C + G) é o polinômio G_ij + C_ij lambda -> [G_ij, C_ij]
    coef = determinante_polinomial([[[G[i, j], C[i, j]] for j in range(n)] for i in range(n)])
    rel.subtitulo("a) Polinômio característico det(lambda C + G) e Bairstow")
    rel.escrever("P(lambda) = " + " + ".join(f"{coef[i]:.8g} l^{i}" for i in range(n, -1, -1)))
    lambdas, hist = bairstow(coef, -1.0, -1.0, tol=1e-12)
    rel.tabela(hist, casas=10)
    lambdas = sorted(float(np.real(l)) for l in lambdas)
    rel.tabela([{"modo": i + 1, "lambda [1/s]": l, "tau = -1/lambda [s]": -1 / l,
                 "resíduo relativo": abs(avaliar_polinomio(coef, l))
                 / sum(abs(c) * abs(l) ** i for i, c in enumerate(coef))}
                for i, l in enumerate(lambdas)], casas=8)
    rel.escrever("resíduo relativo = |P(lambda)| / soma|a_i lambda^i| (independe da escala de lambda)")
    rel.escrever("Todas reais e negativas: rede RC (G simétrica, C > 0) não oscila, só decai.")

    # ---- b) Müller e Vieta -----------------------------------------------------------------
    rel.subtitulo("b) Conferências")
    f = lambda lam: avaliar_polinomio(coef, lam)     # noqa: E731
    for pontos in ((-0.001, -0.005, -0.01), (-0.5, -1.0, -1.5), (-30.0, -50.0, -60.0)):
        raiz, h_m = muller(f, *pontos, tol=1e-12)
        rel.escrever(f"Müller partindo de {pontos}: lambda = {raiz:.10f} em {len(h_m)} iterações")
    traco = sum(G[i, i] / C[i, i] for i in range(n))
    rel.escrever(f"Vieta: soma das raízes = {sum(lambdas):.8f};  -a2/a3 = {-coef[2] / coef[3]:.8f};"
                 f"  -traço(C^-1 G) = {-traco:.8f}")
    rel.escrever(f"Vieta: produto das raízes = {math.prod(lambdas):.8f};  -a0/a3 = {-coef[0] / coef[3]:.8f}")

    # ---- c) validação por RK4 ------------------------------------------------------------
    rel.subtitulo("c) Resfriamento simulado por RK4")
    C_inv = np.array([1.0 / c for c in CAPACITANCIAS])
    t, theta = rk4(lambda t, th: -C_inv * (G @ th), 0.0, THETA_INICIAL, H, int(T_FINAL / H))
    # inclinação de ln(theta_3) no fim da simulação (só o modo lento sobrevive)
    i0, i1 = int(600 / H), int(T_FINAL / H)
    inclinacao = (math.log(theta[i1, 2]) - math.log(theta[i0, 2])) / (t[i1] - t[i0])
    rel.escrever(f"Inclinação de ln(theta_dissipador) entre 600 s e 900 s = {inclinacao:.8f} 1/s")
    rel.escrever(f"lambda_lento (Bairstow)                              = {lambdas[-1]:.8f} 1/s")

    # ---- d) rigidez ------------------------------------------------------------------------
    rel.subtitulo("d) Rigidez (stiffness)")
    razao = lambdas[0] / lambdas[-1]
    rel.escrever(f"tau_lento / tau_rápido = {razao:.0f}. O RK4 (explícito) precisa de h < 2.785/|lambda_rápido|"
                 f" = {2.785 / abs(lambdas[0]):.3f} s por estabilidade, embora a dinâmica de interesse")
    rel.escrever(f"dure ~{-5 / lambdas[-1]:.0f} s: são ~{int(-5 / lambdas[-1] / (2.785 / abs(lambdas[0])))} passos."
                 " Problemas assim pedem métodos implícitos (como no PPC 3).")

    fig, ax = plt.subplots()
    for i, nome in enumerate(("die", "IHS", "dissipador")):
        ax.semilogy(t, theta[:, i], label=nome)
    ax.semilogy(t, theta[i1, 2] * np.exp(lambdas[-1] * (t - t[i1])), "k:", label="e^(λ_lento t)")
    ax.set(xlabel="t [s]", ylabel="T - T∞ [K] (escala log)", ylim=(1e-1, 60),
           title="Resfriamento: os modos rápidos somem e sobra o lento")
    ax.legend()
    salvar_figura(fig, PASTA / "outputs" / "E08_resfriamento.png")

    rel.salvar(PASTA / "outputs" / "E08_resultados.txt")


if __name__ == "__main__":
    main()
