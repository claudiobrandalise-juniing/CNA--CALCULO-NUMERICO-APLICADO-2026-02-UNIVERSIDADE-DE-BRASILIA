"""
verificar_apc3.py
=================
Confere, passo a passo, as contas que a APC 3 pede para fazer À MÃO.
Use para checar o caderno/LaTeX antes de entregar (e como cola de estudo).

    1. Inversa de A por LU e verificação A A^-1 = I
    2. Sistema tridiagonal: LU "cheia" x algoritmo de Thomas (tabela + FLOPs)
    3. Cholesky de uma matriz simétrica 3x3 inventada
    4. Gauss-Seidel sem e com relaxação (lambda = 1.2), erro relativo de 5 %

Executar (a partir da raiz do repositório):
    python APC_03/verificar_apc3.py
"""

import sys
from pathlib import Path

import numpy as np

PASTA = Path(__file__).resolve().parent
sys.path.insert(0, str(PASTA.parent))
from utils.saida import Relatorio  # noqa: E402
from utils.sistemas_iterativos import (dominancia_diagonal, gauss_seidel,  # noqa: E402
                                       reordenar_para_dominancia)
from utils.sistemas_lineares import (cholesky, decomposicao_lu, inversa_lu,  # noqa: E402
                                     novo_contador_flops, resolver_cholesky, resolver_lu,
                                     substituicao_progressiva, substituicao_regressiva,
                                     thomas, tridiagonal_para_vetores)


def item1(rel):
    """Inversa por LU, mostrando L, U e cada coluna (y e x)."""
    rel.subtitulo("Item 1 - Inversa de A por decomposição LU")
    A = np.array([[3.0, -0.1, -0.2], [0.1, 7.0, -0.3], [0.3, -0.2, 10.0]])
    L, U = decomposicao_lu(A)
    rel.matriz(L, "L", casas=8)
    rel.matriz(U, "U", casas=8)
    for j in range(3):
        e = np.zeros(3)
        e[j] = 1.0
        y = substituicao_progressiva(L, e)
        x = substituicao_regressiva(U, y)
        rel.escrever(f"coluna {j + 1}: L y = e{j + 1} -> y = {np.round(y, 8)};  U x = y -> x = {np.round(x, 8)}")
    A_inv = inversa_lu(A)
    rel.matriz(A_inv, "A^-1", casas=8)
    rel.matriz(A @ A_inv, "A A^-1 (deve ser I)", casas=4)
    rel.escrever(f"max|A A^-1 - I| = {np.max(np.abs(A @ A_inv - np.eye(3))):.2e}")


def item2(rel):
    """Tridiagonal 4x4: LU completa x Thomas."""
    rel.subtitulo("Item 2 - Sistema tridiagonal 4x4: LU x Thomas")
    A = np.array([[2.04, -1, 0, 0], [-1, 2.04, -1, 0], [0, -1, 2.04, -1], [0, 0, -1, 2.04]])
    b = np.array([40.8, 0.8, 0.8, 200.8])
    flops_lu = novo_contador_flops()
    L, U = decomposicao_lu(A, flops_lu)
    x_lu = resolver_lu(L, U, b, flops_lu)
    rel.matriz(L, "L", casas=6)
    rel.matriz(U, "U", casas=6)
    rel.escrever(f"Solução por LU: x = {np.round(x_lu, 4)}")
    tabela = []
    flops_thomas = novo_contador_flops()
    e, f, g = tridiagonal_para_vetores(A)
    x_th = thomas(e, f, g, b, flops_thomas, historico=tabela)
    rel.escrever("Tabela dos laços do algoritmo de Thomas:")
    rel.tabela(tabela, casas=8)
    rel.escrever(f"Solução por Thomas: x = {np.round(x_th, 4)}  (diferença para LU: "
                 f"{np.max(np.abs(x_th - x_lu)):.1e})")
    rel.tabela([
        {"método": "LU completa (decomp. + 2 substituições)", "mult/div": flops_lu["mult_div"],
         "soma/sub": flops_lu["soma_sub"], "total": flops_lu["mult_div"] + flops_lu["soma_sub"]},
        {"método": "Thomas", "mult/div": flops_thomas["mult_div"], "soma/sub": flops_thomas["soma_sub"],
         "total": flops_thomas["mult_div"] + flops_thomas["soma_sub"]},
    ])
    rel.escrever("Fórmulas gerais: Thomas = 5(n-1) + 1 mult/div e 3(n-1) soma/sub -> 8n - 7 flops; "
                 "LU completa ~ 2n^3/3 + 2n^2.")


def item3(rel):
    """Cholesky de uma matriz simétrica 3x3 (troque pelos seus números!)."""
    rel.subtitulo("Item 3 - Cholesky (matriz simétrica inventada)")
    A = np.array([[6.0, 15.0, 55.0], [15.0, 55.0, 225.0], [55.0, 225.0, 979.0]])
    b = A @ np.array([1.0, 2.0, 3.0])          # b escolhido para a solução ser (1, 2, 3)
    rel.matriz(A, "A (simétrica)", casas=6)
    rel.escrever(f"b = {b}")
    L = cholesky(A)
    rel.matriz(L, "L (A = L L^T)", casas=8)
    rel.escrever(f"max|L L^T - A| = {np.max(np.abs(L @ L.T - A)):.2e}")
    y = substituicao_progressiva(L, b)
    rel.escrever(f"L y = b      -> y = {np.round(y, 8)}")
    rel.escrever(f"L^T x = y    -> x = {np.round(resolver_cholesky(L, b), 8)}")


def item4(rel):
    """Gauss-Seidel com reordenação, sem e com relaxação."""
    rel.subtitulo("Item 4 - Gauss-Seidel (erro relativo de 5 %)")
    A = [[2, -6, -1], [-3, -1, 7], [-8, 1, -2]]
    b = [-38, -34, -20]
    rel.escrever(f"Ordem original diagonalmente dominante? {dominancia_diagonal(A)}")
    A_re, b_re, ordem = reordenar_para_dominancia(A, b)
    rel.escrever(f"Nova ordem das equações (originais): {[i + 1 for i in ordem]}")
    rel.matriz(A_re, "A reordenada", casas=4)
    rel.escrever(f"b reordenado = {b_re};  dominante agora? {dominancia_diagonal(A_re)}")
    for relaxacao in (1.0, 1.2):
        x, hist = gauss_seidel(A_re, b_re, x0=[0, 0, 0], tol=5.0, relaxacao=relaxacao,
                               criterio="relativo")
        rel.escrever(f"\nlambda = {relaxacao}:  (erro = máximo erro relativo aproximado, %)")
        rel.tabela(hist, casas=8)
        rel.escrever(f"x = {np.round(x, 5)} em {len(hist)} iterações (exata: [4, 8, -2])")


def main():
    rel = Relatorio()
    rel.titulo("APC 03 - Conferência das contas feitas à mão")
    item1(rel)
    item2(rel)
    item3(rel)
    item4(rel)
    rel.salvar(PASTA / "outputs" / "verificacao.txt")


if __name__ == "__main__":
    main()
