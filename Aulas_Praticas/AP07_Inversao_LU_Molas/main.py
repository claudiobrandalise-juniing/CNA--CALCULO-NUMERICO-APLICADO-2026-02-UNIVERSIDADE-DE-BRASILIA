"""
Aula Prática 07 - Inversão matricial por decomposição LU
=========================================================
Matriz de flexibilidade C = K^-1 de uma cadeia de molas (F = K u  ->  u = C F).

    K = [k0+k1   -k1      0     ...  ]     (tridiagonal, simétrica)
        [ -k1   k1+k2   -k2     ...  ]
        [  ...                       ]

Estratégia: decompor K = L U UMA vez e, para cada coluna j da identidade,
resolver L y = e_j (progressiva) e U c_j = y (regressiva). c_j é a coluna j
de C = deslocamentos causados por uma força unitária no nó j.

Executar (a partir da raiz do repositório):
    python Aulas_Praticas/AP07_Inversao_LU_Molas/main.py
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# ---- acesso à biblioteca utils/ (raiz do repositório) ------------------------
PASTA = Path(__file__).resolve().parent
SAIDA = PASTA / "outputs"
sys.path.insert(0, str(PASTA.parents[1]))
from utils.padrao_grafico import configurar_estilo, salvar_figura  # noqa: E402
from utils.saida import Relatorio  # noqa: E402
from utils.sistemas_lineares import decomposicao_lu, inversa_lu, resolver_lu  # noqa: E402

# ---- DADOS DE ENTRADA (enunciado) --------------------------------------------
RIGIDEZ_2GDL = [1000.0, 500.0, 1000.0]                          # k0, k1, k2 [N/m]
RIGIDEZ_5GDL = [1200.0, 800.0, 1500.0, 1000.0, 1800.0, 900.0]  # k0..k5 [N/m]
CARREGAMENTO = np.array([0.0, 100.0, 0.0, -50.0, 25.0])        # F [N]


def matriz_rigidez(k):
    """Monta K (n x n) de uma cadeia parede-mola-nó-...-nó-mola-parede.

    k = [k0, k1, ..., kn]: k0 liga a parede esquerda ao nó 1, k_i liga o nó i
    ao nó i+1 e k_n liga o último nó à parede direita (n + 1 molas, n nós).
    Para o caso de 2 GDL do enunciado, k = [k0, k1, k2] dá exatamente
    K = [[k0+k1, -k1], [-k1, k1+k2]].
    """
    n = len(k) - 1
    K = np.zeros((n, n))
    for i in range(n):
        K[i, i] = k[i] + k[i + 1]          # molas à esquerda e à direita do nó
        if i < n - 1:
            K[i, i + 1] = K[i + 1, i] = -k[i + 1]   # acoplamento com o vizinho
    return K


def verificar_inversa(rel, K, C):
    """Imprime max|K C - I| e max|C - C^T|."""
    n = K.shape[0]
    erro_identidade = np.max(np.abs(K @ C - np.eye(n)))
    erro_simetria = np.max(np.abs(C - C.T))
    rel.escrever(f"erro_max = max|K C - I|   = {erro_identidade:.3e}")
    rel.escrever(f"erro_sim = max|C_ij - C_ji| = {erro_simetria:.3e} m/N")
    return erro_identidade, erro_simetria


def main():
    rel = Relatorio()
    configurar_estilo()
    rel.titulo("AP07 - Matriz de flexibilidade por inversão LU")

    # ---- Parte 1: 2 graus de liberdade -----------------------------------------
    rel.subtitulo("Parte 1 - Sistema com 2 GDL")
    K2 = matriz_rigidez(RIGIDEZ_2GDL)
    rel.matriz(K2, "K [N/m]")
    det_K2 = K2[0, 0] * K2[1, 1] - K2[0, 1] * K2[1, 0]
    C2_manual = np.array([[K2[1, 1], -K2[0, 1]], [-K2[1, 0], K2[0, 0]]]) / det_K2
    rel.escrever(f"Inversa 2x2 'no papel': C = (1/det) [[k22, -k12], [-k21, k11]], det = {det_K2:g}")
    rel.matriz(C2_manual, "C manual [m/N]", casas=8)
    C2 = inversa_lu(K2)
    rel.matriz(C2, "C por LU [m/N]", casas=8)
    verificar_inversa(rel, K2, C2)
    for j in range(2):
        e = np.zeros(2)
        e[j] = 1.0
        rel.escrever(f"F = e_{j + 1} (1 N no nó {j + 1}) -> u = C F = {C2 @ e} m = coluna {j + 1} de C")

    # ---- Partes 2 e 3: 5 graus de liberdade -----------------------------------------
    rel.subtitulo("Partes 2 e 3 - Sistema com 5 GDL")
    K = matriz_rigidez(RIGIDEZ_5GDL)
    rel.matriz(K, "K [N/m]")
    L, U = decomposicao_lu(K)                      # feita UMA única vez
    rel.matriz(L, "L", casas=6)
    rel.matriz(U, "U [N/m]", casas=6)
    rel.escrever(f"max|L U - K| = {np.max(np.abs(L @ U - K)):.3e} N/m")
    C = inversa_lu(K)
    rel.matriz(C * 1e3, "C [mm/N]", casas=6)

    # ---- Parte 4: verificações ----------------------------------------------------
    rel.subtitulo("Parte 4 - Verificação e uso da matriz de flexibilidade")
    verificar_inversa(rel, K, C)
    u = C @ CARREGAMENTO                           # só uma multiplicação matriz-vetor
    u_direto = resolver_lu(L, U, CARREGAMENTO)     # conferência: resolver K u = F
    rel.escrever(f"F = {CARREGAMENTO} N")
    rel.escrever("u = C F = [" + ", ".join(f"{v * 1e3:.6f}" for v in u) + "] mm")
    rel.escrever(f"diferença para a solução direta K u = F: {np.max(np.abs(u - u_direto)):.2e} m")
    SAIDA.mkdir(parents=True, exist_ok=True)
    np.savetxt(SAIDA / "matriz_flexibilidade.txt", C, fmt="%.10e",
               header="Matriz de flexibilidade C = K^-1 [m/N] (linha i, coluna j)")

    # ---- questão f: mola central quase sem rigidez -------------------------------
    rel.subtitulo("Questão f - Mola k2 (entre os nós 2 e 3) quase sem rigidez (1500 -> 1 N/m)")
    k_mole = list(RIGIDEZ_5GDL)
    k_mole[2] = 1.0
    C_mole = inversa_lu(matriz_rigidez(k_mole))
    rel.matriz(C_mole * 1e3, "C [mm/N]", casas=6)

    # ---- questão g: sem apoios nas extremidades -------------------------------------
    rel.subtitulo("Questão g - Sem as paredes (k0 = k5 = 0)")
    k_livre = list(RIGIDEZ_5GDL)
    k_livre[0] = k_livre[-1] = 0.0
    try:
        inversa_lu(matriz_rigidez(k_livre))
    except ZeroDivisionError as erro:
        rel.escrever(f"Inversão falhou: {erro}")

    # ---- gráfico: mapa da matriz de flexibilidade ------------------------------------
    fig, ax = plt.subplots(figsize=(5.2, 4.3))
    imagem = ax.imshow(C * 1e3, cmap="viridis")
    fig.colorbar(imagem, label="C_ij [mm/N]")
    ax.set(xticks=range(5), yticks=range(5), xticklabels=range(1, 6), yticklabels=range(1, 6),
           xlabel="nó j onde a força unitária é aplicada", ylabel="nó i (deslocamento)",
           title="Matriz de flexibilidade C = K⁻¹")
    ax.grid(False)
    salvar_figura(fig, SAIDA / "matriz_flexibilidade.png")

    rel.salvar(SAIDA / "resultados.txt")


if __name__ == "__main__":
    main()
