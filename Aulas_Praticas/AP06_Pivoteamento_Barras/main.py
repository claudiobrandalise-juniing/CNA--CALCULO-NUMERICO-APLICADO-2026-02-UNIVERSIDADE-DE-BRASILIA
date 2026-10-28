"""
Aula Prática 06 - Condicionamento e pivoteamento parcial em estrutura de barras
================================================================================
Equilíbrio de nós articulados (barras só com esforço axial):

  2 barras (nó A):  [sin t1  sin t2] [F1]     [Px]
                    [cos t1  cos t2] [F2] = - [Py]      det = sin(t1 - t2)

  5 barras (nós A, B, C):  A F = b  com a 1a linha começando por sin(t1),
  que é um pivô minúsculo quando as barras 1 e 2 ficam quase verticais.

Objetivos: (i) separar a SENSIBILIDADE DO PROBLEMA (condicionamento) do ERRO
DO ALGORITMO (pivô pequeno) e (ii) implementar o pivoteamento parcial.

Executar (a partir da raiz do repositório):
    python Aulas_Praticas/AP06_Pivoteamento_Barras/main.py
"""

import math
import sys
from pathlib import Path

import numpy as np

# ---- acesso à biblioteca utils/ (raiz do repositório) ------------------------
PASTA = Path(__file__).resolve().parent
SAIDA = PASTA / "outputs"
sys.path.insert(0, str(PASTA.parents[1]))
from utils.metricas_erro import norma_inf, numero_condicao_inf  # noqa: E402
from utils.saida import Relatorio  # noqa: E402
from utils.sistemas_lineares import (gauss_ingenua, gauss_pivoteamento,  # noqa: E402
                                     novo_contador_flops, residuo)

# ---- DADOS DE ENTRADA (enunciado) --------------------------------------------
CONFIGURACOES = {"A (afastadas)": (15.0, 65.0),       # (theta1, theta2) em graus
                 "B (próximas)": (10.0, 10.1),
                 "C (quase paralelas)": (0.01, 0.02)}
PHI = 45.0                                          # ângulo da barra 4 [graus]
F_REF_2 = np.array([100.0, -80.0])                  # esforços de referência, 2 barras [N]
F_REF_5 = np.array([100.0, -80.0, 20.0, 35.0, 15.0])  # idem, 5 barras [N]
PERTURBACAO = 0.01                                  # acréscimo em p_y / b_2 [N]


def matriz_2_barras(theta1, theta2):
    """Matriz A2 do equilíbrio do nó A (ângulos em graus)."""
    t1, t2 = math.radians(theta1), math.radians(theta2)
    return np.array([[math.sin(t1), math.sin(t2)],
                     [math.cos(t1), math.cos(t2)]])


def matriz_5_barras(theta1, theta2, phi=PHI):
    """Matriz A (5x5) do equilíbrio dos nós A, B e C (ângulos em graus)."""
    t1, t2, p = math.radians(theta1), math.radians(theta2), math.radians(phi)
    return np.array([
        [math.sin(t1), math.sin(t2), 1.0, 0.0, 0.0],      # nó A, direção x
        [math.cos(t1), math.cos(t2), 0.0, 0.0, 0.0],      # nó A, direção y
        [0.0, 0.0, -1.0, math.cos(p), 0.0],               # nó B, direção x
        [0.0, 0.0, 0.0, -math.sin(p), -1.0],              # nó B, direção y
        [0.0, 0.0, 0.0, 0.0, 1.0],                        # nó C, direção y
    ])


def inversa_por_colunas(A):
    """A^-1 resolvendo A x = e_j com pivoteamento (usada só para cond(A))."""
    n = A.shape[0]
    return np.column_stack([gauss_pivoteamento(A, np.eye(n)[:, j]) for j in range(n)])


def parte1(rel):
    """Determinante e sensibilidade do sistema 2x2."""
    rel.subtitulo("Parte 1 - Sistema 2x2: determinante e sensibilidade às cargas")
    linhas = []
    for nome, (t1, t2) in CONFIGURACOES.items():
        A2 = matriz_2_barras(t1, t2)
        p = -A2 @ F_REF_2                       # cargas que geram F_ref
        F = gauss_pivoteamento(A2, -p)
        linha = {"config.": nome, "|det| = |sin(t1-t2)|": abs(math.sin(math.radians(t1 - t2)))}
        # perturba uma componente de carga por vez: p_y (pedido) e p_x (complemento)
        for indice, rotulo in ((1, "py"), (0, "px")):
            p_pert = p.copy()
            p_pert[indice] += PERTURBACAO
            dF = gauss_pivoteamento(A2, -p_pert) - F
            linha[f"dF1 ({rotulo})"] = dF[0]
            linha[f"dF2 ({rotulo})"] = dF[1]
        linha["cond_inf"] = numero_condicao_inf(A2, inversa_por_colunas(A2))
        linhas.append(linha)
    rel.tabela(linhas, casas=6)
    rel.escrever(f"(perturbação aplicada: +{PERTURBACAO} N em p_y e, separadamente, em p_x)")

    # o determinante depende da escala; o número de condição não
    rel.subtitulo("Determinante x escala: mesma matriz C escrita em mN (x1000)")
    A2 = matriz_2_barras(*CONFIGURACOES["C (quase paralelas)"])
    for fator, unidade in ((1.0, "N"), (1000.0, "mN")):
        Ae = fator * A2
        det_2x2 = Ae[0, 0] * Ae[1, 1] - Ae[0, 1] * Ae[1, 0]
        rel.escrever(f"  unidade {unidade:>2}: det = {det_2x2:.4e},  "
                     f"cond_inf = {numero_condicao_inf(Ae, inversa_por_colunas(Ae)):.4e}")


def resolver_e_relatar(rel, A, b, metodo, **opcoes):
    """Resolve com o método indicado e imprime erro, resíduo e trocas."""
    trocas = []
    if metodo == "pivoteamento":
        F = gauss_pivoteamento(A, b, trocas=trocas, **opcoes)
    else:
        F = gauss_ingenua(A, b, **opcoes)
    rel.escrever(f"  F calculado = [{', '.join(f'{v:.10g}' for v in F)}] N")
    if metodo == "pivoteamento":
        rel.escrever(f"  trocas de linha (k, p): {trocas if trocas else 'nenhuma'}")
    rel.escrever(f"  ||F - F_ref||_inf = {norma_inf(F - F_REF_5):.3e} N,  "
                 f"||r||_inf = ||A F - b||_inf = {norma_inf(residuo(A, F, b)):.3e}")
    return F


def main():
    rel = Relatorio()
    rel.titulo("AP06 - Condicionamento e pivoteamento parcial (estrutura de barras)")
    parte1(rel)

    # ---- Partes 2 a 4: sistema 5x5 ----------------------------------------------
    for nome in ("A (afastadas)", "C (quase paralelas)"):
        A = matriz_5_barras(*CONFIGURACOES[nome])
        b = A @ F_REF_5                          # cargas consistentes com F_ref
        rel.subtitulo(f"Partes 2-4 - Sistema 5x5, configuração {nome}")
        rel.escrever("Eliminação INGÊNUA:")
        resolver_e_relatar(rel, A, b, "ingenua")
        rel.escrever("Eliminação com PIVOTEAMENTO PARCIAL (etapas):")
        resolver_e_relatar(rel, A, b, "pivoteamento", mostrar_etapas=True, saida=rel.escrever)

    # ---- item 7: precisão reduzida (arredonda CADA operação) ---------------------
    rel.subtitulo("Item 7 - Configuração C com aritmética de 4 algarismos significativos")
    A_C = matriz_5_barras(*CONFIGURACOES["C (quase paralelas)"])
    b_C = A_C @ F_REF_5
    rel.escrever("Ingênua (4 algarismos):")
    resolver_e_relatar(rel, A_C, b_C, "ingenua", algarismos=4)
    rel.escrever("Pivoteamento (4 algarismos):")
    resolver_e_relatar(rel, A_C, b_C, "pivoteamento", algarismos=4)
    rel.escrever("")
    rel.escrever("Efeito do número de algarismos (configuração C):")
    linhas = []
    for algarismos in (4, 6, 8, 10):
        for metodo, solver in (("ingênua", gauss_ingenua), ("pivoteamento", gauss_pivoteamento)):
            F = solver(A_C, b_C, algarismos=algarismos)
            linhas.append({"algarismos": algarismos, "método": metodo,
                           "||F - F_ref||_inf [N]": norma_inf(F - F_REF_5),
                           "||r||_inf": norma_inf(residuo(A_C, F, b_C))})
    rel.tabela(linhas, casas=4)

    # ---- item 8: condicionamento do 5x5 ----------------------------------------
    rel.subtitulo("Item 8 - Sensibilidade do 5x5 (+0.01 N em b_2 e, como complemento, em b_1)")
    linhas = []
    for nome in ("A (afastadas)", "C (quase paralelas)"):
        A = matriz_5_barras(*CONFIGURACOES[nome])
        b = A @ F_REF_5
        F = gauss_pivoteamento(A, b)
        for indice in (1, 0):                    # b_2 (pedido) e b_1 (complemento)
            b_pert = b.copy()
            b_pert[indice] += PERTURBACAO
            dF = gauss_pivoteamento(A, b_pert) - F
            linhas.append({"config.": nome, "perturbado": f"b_{indice + 1}",
                           "||db||_inf": PERTURBACAO, "||dF||_inf": norma_inf(dF),
                           "amplificação": norma_inf(dF) / PERTURBACAO,
                           "cond_inf(A)": numero_condicao_inf(A, inversa_por_colunas(A))})
    rel.tabela(linhas, casas=6)

    # ---- item 9: FLOPs ----------------------------------------------------------
    rel.subtitulo("Item 9 - FLOPs para n = 5 (configuração C)")
    for nome, solver in (("ingênua", gauss_ingenua), ("pivoteamento", gauss_pivoteamento)):
        flops = novo_contador_flops()
        solver(A_C, b_C, flops=flops)
        rel.escrever(f"  {nome:>13}: {flops['mult_div']} mult/div + {flops['soma_sub']} soma/sub "
                     f"= {flops['mult_div'] + flops['soma_sub']} FLOPs")

    rel.salvar(SAIDA / "resultados.txt")


if __name__ == "__main__":
    main()
