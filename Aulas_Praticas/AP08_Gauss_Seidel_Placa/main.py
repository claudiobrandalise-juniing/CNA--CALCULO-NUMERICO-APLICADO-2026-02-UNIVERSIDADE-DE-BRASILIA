"""
Aula Prática 08 - Método de Gauss-Seidel: condução de calor numa placa
========================================================================
Placa quadrada (L = 1 m), regime permanente, sem geração:

    d2T/dx2 + d2T/dy2 = 0

Malha uniforme (dx = dy) e diferenças centrais -> em cada nó interno a
temperatura é a média dos 4 vizinhos:

    T(i,j) = [T(i-1,j) + T(i+1,j) + T(i,j-1) + T(i,j+1)] / 4

Bordas: inferior a 100 °C; superior, esquerda e direita a 20 °C.
O sistema é resolvido por Gauss-Seidel direto na malha (sem montar a matriz).

Esta versão corrige o "PPC SURPRESA" original: lá a matriz era criada com
np.full(..., 20) (INTEIRA), e cada média era truncada para inteiro. Por isso
o código parava em 14 iterações com o centro a 34 °C. Com float, o centro dá
40 °C, que é o valor exato por simetria (ver validação no README).

Convenção: T[i, j] com i = direção x (horizontal) e j = direção y (vertical);
j = 0 é a borda de baixo, igual ao pseudocódigo do enunciado.

Executar (a partir da raiz do repositório):
    python Aulas_Praticas/AP08_Gauss_Seidel_Placa/main.py
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# ---- acesso à biblioteca utils/ (raiz do repositório) ------------------------
PASTA = Path(__file__).resolve().parent
SAIDA = PASTA / "outputs"
sys.path.insert(0, str(PASTA.parents[1]))
from utils.edp import gauss_seidel_laplace_2d, residuo_laplace_2d  # noqa: E402
from utils.padrao_grafico import configurar_estilo, salvar_figura  # noqa: E402
from utils.saida import Relatorio  # noqa: E402

# ---- DADOS DE ENTRADA (enunciado) --------------------------------------------
T_INFERIOR = 100.0          # borda de baixo [°C]
T_OUTRAS = 20.0             # bordas superior, esquerda e direita [°C]
TOL = 1e-6                  # maior variação entre varreduras [°C]
MAX_ITER = 100000
MALHAS = [5, 10, 20]        # número n de nós internos em cada direção


def malha_inicial(n):
    """Matriz (n+2) x (n+2) com as bordas já impostas e o interior a 20 °C.

    O dtype=float é ESSENCIAL: com inteiros as médias seriam truncadas.
    """
    T = np.full((n + 2, n + 2), T_OUTRAS, dtype=float)
    T[1:n + 1, 0] = T_INFERIOR      # borda inferior (j = 0), sem os cantos
    return T


def main():
    rel = Relatorio()
    configurar_estilo()
    rel.titulo("AP08 - Gauss-Seidel: condução 2D numa placa quadrada")

    resultados, campos, historicos = [], {}, {}
    for n in MALHAS:
        T, erros = gauss_seidel_laplace_2d(malha_inicial(n), tol=TOL, max_iter=MAX_ITER)
        convergiu = erros[-1] < TOL
        centro = (n + 1) // 2           # nó central (exato quando n é ímpar)
        resultados.append({"n": n, "iterações": len(erros), "erro final [°C]": erros[-1],
                           "convergiu": convergiu,
                           "T no nó ((n+1)/2, (n+1)/2) [°C]": T[centro, centro],
                           "max|r_ij| [°C]": residuo_laplace_2d(T)})
        campos[n], historicos[n] = T, erros
        # matriz completa (com bordas) em arquivo, na orientação física:
        # 1a linha do arquivo = borda de CIMA, última linha = borda de BAIXO
        SAIDA.mkdir(parents=True, exist_ok=True)
        np.savetxt(SAIDA / f"temperaturas_n{n}.txt", T.T[::-1], fmt="%9.4f",
                   header=f"Temperaturas [°C], n = {n}. Linhas: de cima (y = L) para baixo "
                          "(y = 0); colunas: da esquerda (x = 0) para a direita (x = L).")

    rel.tabela(resultados, casas=8)
    rel.escrever("\nNó (4,4) do enunciado para n = 5 (índice [3, 3] no Python): "
                 f"T = {campos[5][3, 3]:.6f} °C")
    rel.escrever("Referência por simetria/superposição: (100 + 20 + 20 + 20)/4 = 40 °C")
    rel.subtitulo("Matriz de temperaturas para n = 5 (cima -> baixo)")
    rel.matriz(campos[5].T[::-1], casas=6)

    # ---- gráficos -------------------------------------------------------------
    fig, eixos = plt.subplots(1, 3, figsize=(12, 3.8))
    for ax, n in zip(eixos, MALHAS):
        imagem = ax.imshow(campos[n].T, origin="lower", extent=(0, 1, 0, 1),
                           cmap="inferno", vmin=T_OUTRAS, vmax=T_INFERIOR)
        ax.set(title=f"n = {n} ({len(historicos[n])} iterações)", xlabel="x [m]", ylabel="y [m]")
        ax.grid(False)
    fig.colorbar(imagem, ax=eixos, label="T [°C]", shrink=0.85)
    fig.suptitle("Distribuição de temperatura na placa (borda inferior a 100 °C)")
    fig.savefig(SAIDA / "campos_temperatura.png", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots()
    for n in MALHAS:
        ax.semilogy(range(1, len(historicos[n]) + 1), historicos[n], label=f"n = {n}")
    ax.axhline(TOL, color="k", ls="--", lw=0.8, label="tolerância")
    ax.set(xlabel="varredura", ylabel="max |T_novo - T_antigo| [°C]",
           title="Convergência de Gauss-Seidel x refinamento da malha")
    ax.legend()
    salvar_figura(fig, SAIDA / "convergencia.png")

    rel.salvar(SAIDA / "resultados.txt")


if __name__ == "__main__":
    main()

# =============================================================================
# QUESTÕES RÁPIDAS (respondidas no código, como pede o enunciado)
# =============================================================================
# a) Por que as temperaturas das bordas não são atualizadas?
#    Porque são CONDIÇÕES DE CONTORNO de Dirichlet: temperaturas impostas
#    fisicamente (mantidas constantes por algo externo). Elas não são
#    incógnitas do sistema, entram só como dados nas equações dos nós vizinhos.
#    Por isso o laço percorre apenas os nós internos (1..n).
#
# b) O que diferencia Gauss-Seidel de Jacobi?
#    Jacobi calcula TODOS os novos valores usando só a varredura anterior (precisa
#    de uma cópia da matriz). Gauss-Seidel usa IMEDIATAMENTE cada valor recém-
#    atualizado no cálculo do próximo nó, na mesma matriz T. Isso propaga a
#    informação mais rápido, converge em cerca de metade das iterações de
#    Jacobi e economiza memória.
#
# c) Como o número de iterações muda quando a malha é refinada?
#    Cresce bastante: n = 5 -> 56, n = 10 -> 178, n = 20 -> 598 varreduras
#    (aproximadamente proporcional a n^2). Cada varredura só leva informação de
#    um nó para o vizinho, então a influência da borda quente precisa de cerca
#    de n passos para atravessar a placa, e o erro decai mais devagar a cada
#    varredura porque o raio espectral da iteração se aproxima de 1. Além
#    disso, cada varredura custa n^2 atualizações, então o custo total cresce
#    ~n^4. A sobre-relaxação (SOR, lambda entre 1 e 2) reduz muito esse número.
