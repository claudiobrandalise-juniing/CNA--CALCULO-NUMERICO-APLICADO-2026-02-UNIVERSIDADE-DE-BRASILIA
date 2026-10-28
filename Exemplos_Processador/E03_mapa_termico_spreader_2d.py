"""
E03 - Mapa térmico 2D do spreader de cobre com 4 núcleos  [GAUSS-SEIDEL / SOR + FONTE]
======================================================================================
ENUNCIADO (estilo prova)
    O spreader (IHS) de cobre é uma placa quadrada de 40 mm x 40 mm, espessura
    t = 3 mm, k = 390 W/(m K). As bordas estão presas à base do dissipador a 60 °C.
    Quatro núcleos de 6 mm x 6 mm, no centro, injetam P = [30, 25, 20, 15] W.
    Regime permanente, condução 2D com geração (Poisson):

        d2T/dx2 + d2T/dy2 + q'''/k = 0   ->   T_ij = (soma dos 4 vizinhos + S_ij)/4
        S_ij = q''' dx^2 / k = P_no / (k t)   (P_no = potência atribuída ao nó)

    a) Calcule o campo de temperatura por Gauss-Seidel e encontre o ponto mais quente.
    b) Compare o número de iterações para lambda = 1.0, 1.5 e 1.85 (SOR).
    c) Verifique a conservação de energia (calor que sai pelas bordas = 90 W).

MÉTODO: utils/edp.py -> gauss_seidel_laplace_2d (com fonte e relaxação)
COMO ADAPTAR: mude dimensões, temperatura das bordas e a função fonte_nucleos.

Executar (a partir da raiz do repositório):
    python Exemplos_Processador/E03_mapa_termico_spreader_2d.py      (~10 s)
"""

import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

PASTA = Path(__file__).resolve().parent
sys.path.insert(0, str(PASTA.parent))
from utils.edp import gauss_seidel_laplace_2d, residuo_laplace_2d  # noqa: E402
from utils.padrao_grafico import configurar_estilo, salvar_figura  # noqa: E402
from utils.saida import Relatorio  # noqa: E402

# ---- DADOS -----------------------------------------------------------------------------
LADO = 0.040            # lado da placa [m]
ESPESSURA = 0.003       # t [m]
K_COBRE = 390.0         # condutividade [W/(m K)]
T_BORDA = 60.0          # bordas na base do dissipador [°C]
DX = 0.001              # espaçamento da malha [m] (dx = dy)
# núcleos: (x_min, x_max, y_min, y_max) em metros e potência [W]
NUCLEOS = [((0.013, 0.019, 0.021, 0.027), 30.0),
           ((0.021, 0.027, 0.021, 0.027), 25.0),
           ((0.013, 0.019, 0.013, 0.019), 20.0),
           ((0.021, 0.027, 0.013, 0.019), 15.0)]
TOL = 1e-6              # maior variação entre varreduras [°C]
RELAXACOES = [1.0, 1.5, 1.85]


def fonte_nucleos(x, y):
    """Termo S_ij = P_no/(k t) de cada nó; a potência de cada núcleo é dividida
    igualmente entre os nós que caem dentro dele (assim a energia é exata)."""
    S = np.zeros((len(x), len(y)))
    for (x0, x1, y0, y1), potencia in NUCLEOS:
        dentro = ((x[:, None] >= x0 - 1e-12) & (x[:, None] <= x1 + 1e-12)
                  & (y[None, :] >= y0 - 1e-12) & (y[None, :] <= y1 + 1e-12))
        S[dentro] += potencia / dentro.sum() / (K_COBRE * ESPESSURA)
    return S


def calor_pelas_bordas(T):
    """Soma k t (T_vizinho_interno - T_borda) em todas as ligações nó interno-borda.
    Pela forma discreta do teorema da divergência, isso é exatamente o calor
    gerado (a menos da tolerância do Gauss-Seidel)."""
    diferencas = ((T[1, 1:-1] - T[0, 1:-1]).sum() + (T[-2, 1:-1] - T[-1, 1:-1]).sum()
                  + (T[1:-1, 1] - T[1:-1, 0]).sum() + (T[1:-1, -2] - T[1:-1, -1]).sum())
    return K_COBRE * ESPESSURA * diferencas


def main():
    rel = Relatorio()
    configurar_estilo()
    rel.titulo("E03 - Mapa térmico 2D do spreader (Gauss-Seidel/SOR com geração)")
    n_nos = int(round(LADO / DX)) + 1
    x = np.linspace(0.0, LADO, n_nos)
    y = np.linspace(0.0, LADO, n_nos)
    S = fonte_nucleos(x, y)
    T_inicial = np.full((n_nos, n_nos), T_BORDA)          # bordas e chute = 60 °C
    rel.escrever(f"Malha {n_nos} x {n_nos} nós (dx = {DX * 1000:.1f} mm), "
                 f"{(n_nos - 2) ** 2} incógnitas; potência total = {sum(p for _, p in NUCLEOS)} W")
    lambda_otimo = 2.0 / (1.0 + math.sin(math.pi / (n_nos - 1)))
    rel.escrever(f"lambda ótimo teórico do SOR (Laplace, malha quadrada) = 2/(1+sin(pi/N)) = "
                 f"{lambda_otimo:.3f}")

    # ---- b) efeito da relaxação -------------------------------------------------------
    rel.subtitulo("b) Iterações x relaxação")
    linhas, historicos = [], {}
    for relaxacao in RELAXACOES:
        T, erros = gauss_seidel_laplace_2d(T_inicial, tol=TOL, relaxacao=relaxacao, fonte=S)
        historicos[relaxacao] = erros
        linhas.append({"lambda": relaxacao, "iterações": len(erros), "T máx [°C]": T.max(),
                       "max|r_ij|": residuo_laplace_2d(T, S)})
    rel.tabela(linhas, casas=8)

    # ---- a) ponto mais quente --------------------------------------------------------
    i_max, j_max = np.unravel_index(np.argmax(T), T.shape)
    rel.subtitulo("a) Ponto mais quente")
    rel.escrever(f"T_max = {T.max():.3f} °C em x = {x[i_max] * 1000:.0f} mm, y = {y[j_max] * 1000:.0f} mm"
                 f" (dentro do núcleo de 30 W)")
    for (x0, x1, y0, y1), potencia in NUCLEOS:
        regiao = T[(x >= x0 - 1e-12) & (x <= x1 + 1e-12)][:, (y >= y0 - 1e-12) & (y <= y1 + 1e-12)]
        rel.escrever(f"   núcleo de {potencia:4.0f} W: T média = {regiao.mean():.2f} °C, "
                     f"T máx = {regiao.max():.2f} °C")

    # ---- c) conservação de energia -----------------------------------------------------
    rel.subtitulo("c) Conservação de energia")
    rel.escrever(f"Calor que sai pelas bordas = {calor_pelas_bordas(T):.6f} W  "
                 f"(gerado = {sum(p for _, p in NUCLEOS):.6f} W)")

    # ---- gráficos -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(6, 5))
    imagem = ax.imshow(T.T, origin="lower", extent=(0, LADO * 1000, 0, LADO * 1000), cmap="inferno")
    ax.contour(x * 1000, y * 1000, T.T, levels=10, colors="white", linewidths=0.5)
    fig.colorbar(imagem, label="T [°C]")
    ax.set(xlabel="x [mm]", ylabel="y [mm]", title="Spreader de cobre: temperatura em regime permanente")
    ax.grid(False)
    salvar_figura(fig, PASTA / "outputs" / "E03_mapa_termico.png")

    fig, ax = plt.subplots()
    for relaxacao, erros in historicos.items():
        ax.semilogy(range(1, len(erros) + 1), erros, label=f"λ = {relaxacao}")
    ax.axhline(TOL, color="k", ls="--", lw=0.8)
    ax.set(xlabel="varredura", ylabel="max |ΔT| [°C]", title="Convergência: Gauss-Seidel x SOR")
    ax.legend()
    salvar_figura(fig, PASTA / "outputs" / "E03_convergencia_sor.png")

    rel.salvar(PASTA / "outputs" / "E03_resultados.txt")


if __name__ == "__main__":
    main()
