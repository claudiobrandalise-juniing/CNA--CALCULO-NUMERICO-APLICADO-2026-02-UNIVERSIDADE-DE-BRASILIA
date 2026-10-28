"""
Aula Prática 01 - Métodos da Bisseção e da Falsa Posição
=========================================================
Profundidade normal H de um canal retangular aberto (equação de Manning).

    Q = (S^(1/2)/n) (B H)^(5/3) / (B + 2H)^(2/3)
    f(H) = (S^(1/2)/n) (B H)^(5/3) / (B + 2H)^(2/3) - Q = 0

Executar (a partir da raiz do repositório):
    python Aulas_Praticas/AP01_Bissecao_FalsaPosicao/main.py
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# ---- acesso à biblioteca utils/ (raiz do repositório) ------------------------
PASTA = Path(__file__).resolve().parent
SAIDA = PASTA / "outputs"            # pasta dos resultados
sys.path.insert(0, str(PASTA.parents[1]))
from utils.padrao_grafico import configurar_estilo, salvar_figura  # noqa: E402
from utils.raizes import bissecao, falsa_posicao, falsa_posicao_modificada  # noqa: E402
from utils.saida import Relatorio  # noqa: E402

# ---- DADOS DE ENTRADA (enunciado) --------------------------------------------
S = 0.0002          # inclinação do canal [-]
B = 20.0            # largura do canal [m]
N_MANNING = 0.03    # coeficiente de rugosidade de Manning [s/m^(1/3)]
Q = 5.0             # vazão volumétrica [m3/s]
H_L, H_U = 0.5, 1.0  # intervalo inicial [m]
TOL = 1e-6          # critério: |f(Hr)| < TOL
MAX_ITER = 100      # proteção contra laço infinito


def f_manning(H):
    """Resíduo da equação de Manning: vazão calculada - vazão prescrita [m3/s].

    O coeficiente S^(1/2)/n = 0.471405 é o mesmo do enunciado.
    """
    return (S ** 0.5 / N_MANNING) * (B * H) ** (5.0 / 3.0) / (B + 2.0 * H) ** (2.0 / 3.0) - Q


def main():
    rel = Relatorio()
    configurar_estilo()

    # ---- 1. intervalo inicial: a troca de sinal garante uma raiz -------------
    rel.titulo("AP01 - Profundidade normal em canal retangular (Manning)")
    rel.escrever(f"Coeficiente S^(1/2)/n = {S ** 0.5 / N_MANNING:.6f}")
    f_l, f_u = f_manning(H_L), f_manning(H_U)
    rel.escrever(f"f(Hl = {H_L}) = {f_l:.6f}   f(Hu = {H_U}) = {f_u:.6f}")
    rel.escrever(f"f(Hl)*f(Hu) = {f_l * f_u:.6f} < 0 -> existe raiz no intervalo")

    # ---- Parte 1: bisseção ----------------------------------------------------
    rel.subtitulo("Parte 1 - Método da bisseção")
    H_bis, hist_bis = bissecao(f_manning, H_L, H_U, tol=TOL, max_iter=MAX_ITER, criterio="funcao")
    rel.tabela(hist_bis, colunas=["k", "xl", "xu", "xr", "f(xr)"])
    rel.escrever(f"\nProfundidade (bisseção): H = {H_bis:.6f} m em {len(hist_bis)} iterações")

    # ---- Parte 2: falsa posição ----------------------------------------------
    rel.subtitulo("Parte 2 - Método da falsa posição")
    H_fp, hist_fp = falsa_posicao(f_manning, H_L, H_U, tol=TOL, max_iter=MAX_ITER, criterio="funcao")
    rel.tabela(hist_fp, colunas=["k", "xl", "xu", "xr", "f(xr)"])
    rel.escrever(f"\nProfundidade (falsa posição): H = {H_fp:.6f} m em {len(hist_fp)} iterações")

    # ---- extra: falsa posição modificada (Chapra 5.4) -----------------------
    H_fpm, hist_fpm = falsa_posicao_modificada(f_manning, H_L, H_U, tol=TOL, max_iter=MAX_ITER)

    # ---- 9. velocidade média U = Q/(B H) -----------------------------------
    rel.subtitulo("9 - Velocidade média do escoamento")
    U_bis, U_fp = Q / (B * H_bis), Q / (B * H_fp)
    rel.escrever(f"U (bisseção)     = {U_bis:.6f} m/s")
    rel.escrever(f"U (falsa posição) = {U_fp:.6f} m/s")

    # ---- 10. comparação ----------------------------------------------------
    rel.subtitulo("10 - Comparação dos métodos")
    rel.tabela([
        {"método": "bisseção", "H [m]": H_bis, "iterações": len(hist_bis), "U [m/s]": U_bis},
        {"método": "falsa posição", "H [m]": H_fp, "iterações": len(hist_fp), "U [m/s]": U_fp},
        {"método": "falsa pos. modificada", "H [m]": H_fpm, "iterações": len(hist_fpm),
         "U [m/s]": Q / (B * H_fpm)},
    ], casas=8)
    rel.escrever(f"Diferença |H_bis - H_fp| = {abs(H_bis - H_fp):.2e} m")

    # ---- gráficos -------------------------------------------------------------
    H_graf = np.linspace(0.05, 1.5, 300)
    fig, ax = plt.subplots()
    ax.plot(H_graf, f_manning(H_graf), label="f(H)")
    ax.axhline(0.0, color="k", lw=0.8)
    ax.axvspan(H_L, H_U, alpha=0.12, label="intervalo inicial")
    ax.plot(H_fp, 0.0, "o", label=f"raiz H = {H_fp:.4f} m")
    ax.set(xlabel="H [m]", ylabel="f(H) [m³/s]", title="Método gráfico: resíduo da equação de Manning")
    ax.legend()
    salvar_figura(fig, SAIDA / "funcao_manning.png")

    fig, ax = plt.subplots()
    for nome, hist in (("bisseção", hist_bis), ("falsa posição", hist_fp),
                       ("falsa pos. modificada", hist_fpm)):
        ax.semilogy([h["k"] for h in hist], [abs(h["f(xr)"]) for h in hist], "o-", label=nome)
    ax.axhline(TOL, color="k", ls="--", lw=0.8, label="tolerância")
    ax.set(xlabel="iteração k", ylabel="|f(Hr)| [m³/s]", title="Convergência dos métodos intervalares")
    ax.legend()
    salvar_figura(fig, SAIDA / "convergencia.png")

    rel.salvar(SAIDA / "resultados.txt")


if __name__ == "__main__":
    main()
