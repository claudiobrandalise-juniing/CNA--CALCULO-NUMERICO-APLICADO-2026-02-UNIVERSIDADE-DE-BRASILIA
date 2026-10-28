"""
E02 - Rede de resistências térmicas de um processador de 4 núcleos  [SISTEMAS LINEARES]
=======================================================================================
ENUNCIADO (estilo prova)
    Quatro núcleos em linha (nós 1-4) dissipam P = [30, 25, 20, 15] W. Cada núcleo
    troca calor com os vizinhos (G_lat) e com o spreader/IHS (nó 5, G_nucleo).
    O IHS liga-se ao dissipador (nó 6) pela pasta (G_tim) e o dissipador ao ar
    (T_amb conhecida) por G_ar. Em regime permanente, em cada nó i:

        soma_j G_ij (T_j - T_i) + P_i = 0     (o que entra = o que sai)

    Isso dá G T = b, com G simétrica e positiva definida (matriz de condutâncias).
    a) Resolva por Gauss com pivoteamento, LU, Cholesky e Gauss-Seidel.
    b) Verifique a conservação de energia (calor que vai ao ar = soma das P).
    c) Monte Z = G^-1 (matriz de "resistências de transferência", em K/W) e
       diga quanto sobe cada nó se o núcleo 2 dissipar +10 W (turbo).

MÉTODOS: utils/sistemas_lineares.py e utils/sistemas_iterativos.py
COMO ADAPTAR: escreva o balanço de cada nó, monte G e b com a função
    montar_sistema (só a lista de ligações e as potências mudam).

Executar (a partir da raiz do repositório):
    python Exemplos_Processador/E02_rede_resistencias_multinucleo.py
"""

import sys
from pathlib import Path

import numpy as np

PASTA = Path(__file__).resolve().parent
sys.path.insert(0, str(PASTA.parent))
from utils.saida import Relatorio  # noqa: E402
from utils.sistemas_iterativos import dominancia_diagonal, gauss_seidel  # noqa: E402
from utils.sistemas_lineares import (cholesky, decomposicao_lu, gauss_pivoteamento,  # noqa: E402
                                     inversa_lu, resolver_cholesky, resolver_lu, residuo)

# ---- DADOS --------------------------------------------------------------------------
POTENCIAS = [30.0, 25.0, 20.0, 15.0, 0.0, 0.0]   # P_i [W] (nós 5 e 6 não geram calor)
T_AMB = 35.0                                     # ar do gabinete [°C]
G_LAT, G_NUCLEO, G_TIM, G_AR = 1.5, 4.0, 12.0, 2.5   # condutâncias [W/K]
# ligações internas (nó i, nó j, condutância) com numeração 1..6
LIGACOES = [(1, 2, G_LAT), (2, 3, G_LAT), (3, 4, G_LAT),
            (1, 5, G_NUCLEO), (2, 5, G_NUCLEO), (3, 5, G_NUCLEO), (4, 5, G_NUCLEO),
            (5, 6, G_TIM)]
LIGACOES_AMBIENTE = [(6, G_AR)]                  # (nó, condutância para T_AMB)
NOMES = ["núcleo 1", "núcleo 2", "núcleo 3", "núcleo 4", "IHS", "dissipador"]


def montar_sistema(n, ligacoes, ligacoes_ambiente, potencias, T_amb):
    """Monta G T = b pelo balanço de energia de cada nó.

    Uma ligação G entre i e j soma G na diagonal de i e de j e -G em (i,j), (j,i).
    Uma ligação com o ambiente soma G na diagonal e G*T_amb em b (lado direito).
    """
    G = np.zeros((n, n))
    b = np.array(potencias, dtype=float)
    for i, j, g in ligacoes:
        i, j = i - 1, j - 1                     # numeração 1..n -> índice Python
        G[i, i] += g
        G[j, j] += g
        G[i, j] -= g
        G[j, i] -= g
    for i, g in ligacoes_ambiente:
        G[i - 1, i - 1] += g
        b[i - 1] += g * T_amb
    return G, b


def main():
    rel = Relatorio()
    rel.titulo("E02 - Rede de resistências térmicas (sistemas lineares)")
    n = len(POTENCIAS)
    G, b = montar_sistema(n, LIGACOES, LIGACOES_AMBIENTE, POTENCIAS, T_AMB)
    rel.matriz(G, "G [W/K]", casas=5)
    rel.escrever(f"b = {b} [W]")
    rel.escrever(f"G simétrica? {np.allclose(G, G.T)}   diagonalmente dominante? "
                 f"{dominancia_diagonal(G)}")

    # ---- a) quatro métodos ----------------------------------------------------------
    T_gauss = gauss_pivoteamento(G, b)
    L, U = decomposicao_lu(G)
    T_lu = resolver_lu(L, U, b)
    T_chol = resolver_cholesky(cholesky(G), b)
    T_gs, hist_gs = gauss_seidel(G, b, x0=np.full(n, T_AMB), tol=1e-10, criterio="absoluto")
    rel.subtitulo("a) Temperaturas [°C]")
    rel.tabela([{"nó": NOMES[i], "Gauss+pivô": T_gauss[i], "LU": T_lu[i], "Cholesky": T_chol[i],
                 "Gauss-Seidel": T_gs[i]} for i in range(n)], casas=9)
    rel.escrever(f"Gauss-Seidel: {len(hist_gs)} iterações; maior diferença entre métodos = "
                 f"{max(np.max(np.abs(T - T_gauss)) for T in (T_lu, T_chol, T_gs)):.1e} °C")
    rel.escrever(f"||G T - b||_inf = {np.max(np.abs(residuo(G, T_gauss, b))):.1e} W")
    rel.escrever("Obs.: G não é ESTRITAMENTE dominante (linha do IHS: 28 = 4*4 + 12), mas é")
    rel.escrever("simétrica positiva definida, e isso também garante a convergência de Gauss-Seidel.")
    for relaxacao in (1.5, 1.8):
        _, hist_sor = gauss_seidel(G, b, x0=np.full(n, T_AMB), tol=1e-10, criterio="absoluto",
                                   relaxacao=relaxacao)
        rel.escrever(f"Com sobre-relaxação lambda = {relaxacao}: {len(hist_sor)} iterações")

    # ---- b) conservação de energia ------------------------------------------------------
    rel.subtitulo("b) Conservação de energia")
    calor_ar = G_AR * (T_gauss[5] - T_AMB)
    rel.escrever(f"Calor do dissipador para o ar = {calor_ar:.6f} W;  soma das potências = "
                 f"{sum(POTENCIAS):.6f} W")
    rel.escrever(f"Núcleo mais quente: {NOMES[int(np.argmax(T_gauss[:4]))]} a {np.max(T_gauss[:4]):.2f} °C")

    # ---- c) matriz de resistências de transferência ----------------------------------------
    rel.subtitulo("c) Z = G^-1 [K/W]: Z_ij = aumento de T_i por 1 W no nó j")
    Z = inversa_lu(G)
    rel.matriz(Z, "Z", casas=5)
    aumento = Z[:, 1] * 10.0                     # +10 W no núcleo 2 (coluna 2 de Z)
    rel.tabela([{"nó": NOMES[i], "aumento [K]": aumento[i], "nova T [°C]": T_gauss[i] + aumento[i]}
                for i in range(n)], casas=6)

    rel.salvar(PASTA / "outputs" / "E02_resultados.txt")


if __name__ == "__main__":
    main()
