"""
Aula Prática 05 - Eliminação gaussiana ingênua: rede de 5 reatores de mistura
==============================================================================
Balanço de massa da espécie química em cada reator (regime permanente):

    Q01 c01 + Q31 c3        = (Q12 + Q15) c1
    Q12 c1                  = (Q23 + Q24 + Q25) c2
    Q03 c03 + Q23 c2        = (Q31 + Q34) c3
    Q24 c2 + Q34 c3 + Q54 c5 = Q44 c4
    Q15 c1 + Q25 c2         = (Q54 + Q55) c5

-> sistema linear A c = b resolvido por eliminação gaussiana SEM pivoteamento,
   contando as operações de ponto flutuante (FLOPs).

Executar (a partir da raiz do repositório):
    python Aulas_Praticas/AP05_Gauss_Ingenua_Reatores/main.py
"""

import sys
from pathlib import Path

import numpy as np

# ---- acesso à biblioteca utils/ (raiz do repositório) ------------------------
PASTA = Path(__file__).resolve().parent
SAIDA = PASTA / "outputs"
sys.path.insert(0, str(PASTA.parents[1]))
from utils.metricas_erro import norma_inf  # noqa: E402
from utils.saida import Relatorio  # noqa: E402
from utils.sistemas_lineares import (gauss_ingenua, mostrar_sistema,  # noqa: E402
                                     novo_contador_flops, residuo)

# ---- DADOS DE ENTRADA (enunciado) --------------------------------------------
# vazões Qij [m3/min] (i = origem, j = destino; 0 = alimentação externa)
VAZAO = {"01": 5.0, "03": 8.0, "12": 3.0, "15": 3.0, "23": 1.0, "24": 1.0, "25": 1.0,
         "31": 1.0, "34": 8.0, "44": 11.0, "54": 2.0, "55": 2.0}
C01, C03 = 10.0, 20.0     # concentrações das alimentações externas [mg/m3]
VALORES_ESPERADOS = [11.5094, 11.5094, 19.0566, 16.9983, 11.5094]   # para conferência


def montar_sistema(Q):
    """Monta A e b a partir dos balanços (linha i = reator i, coluna j = c_j)."""
    A = np.array([
        [Q["12"] + Q["15"], 0.0, -Q["31"], 0.0, 0.0],
        [-Q["12"], Q["23"] + Q["24"] + Q["25"], 0.0, 0.0, 0.0],
        [0.0, -Q["23"], Q["31"] + Q["34"], 0.0, 0.0],
        [0.0, -Q["24"], -Q["34"], Q["44"], -Q["54"]],
        [-Q["15"], -Q["25"], 0.0, 0.0, Q["54"] + Q["55"]],
    ])
    b = np.array([Q["01"] * C01, 0.0, Q["03"] * C03, 0.0, 0.0])
    return A, b


def flops_teoricos(n):
    """Contagem teórica pelos laços aninhados (eliminação + substituição).

    Eliminação, etapa k (m = n - k linhas abaixo do pivô):
        por linha: 1 divisão (fator) + m mult. em j + 1 mult. em b -> m + 2
                   m subtrações em j + 1 subtração em b              -> m + 1
    Substituição regressiva: n divisões + n(n-1)/2 mult. e n(n-1)/2 subtrações.
    """
    mult_div = sum(m * (m + 2) for m in range(1, n)) + n + n * (n - 1) // 2
    soma_sub = sum(m * (m + 1) for m in range(1, n)) + n * (n - 1) // 2
    return mult_div, soma_sub


def main():
    rel = Relatorio()
    rel.titulo("AP05 - Eliminação gaussiana ingênua: rede de 5 reatores")

    # ---- Parte 1: sistema original --------------------------------------------
    A, b = montar_sistema(VAZAO)
    n = len(b)
    rel.subtitulo("Parte 1 - Sistema original [A | b]")
    mostrar_sistema(A, b, rel.escrever)
    rel.escrever(f"Exemplo: A[4,3] = {A[3, 2]:g} multiplica c3 na 4a equação.")

    # ---- Parte 2: primeira eliminação "no papel" -------------------------------
    rel.subtitulo("Parte 2 - Primeira etapa (conferência manual)")
    m21 = A[1, 0] / A[0, 0]
    rel.escrever(f"m21 = A21/A11 = {A[1, 0]:g}/{A[0, 0]:g} = {m21:g}")
    rel.escrever(f"L2 <- L2 - m21 L1 = {A[1] - m21 * A[0]}  | b2 = {b[1] - m21 * b[0]:g}")

    # ---- Partes 3 e 4: eliminação + substituição, contando FLOPs ---------------
    rel.subtitulo("Parte 3 - Eliminação progressiva (matriz após cada pivô)")
    flops = novo_contador_flops()
    c = gauss_ingenua(A, b, flops=flops, mostrar_etapas=True, saida=rel.escrever)

    rel.subtitulo("Parte 4 - Concentrações (substituição retroativa)")
    rel.tabela([{"reator": i + 1, "c [mg/m3]": c[i], "esperado": VALORES_ESPERADOS[i]}
                for i in range(n)], casas=8)

    # ---- Parte 6: verificação pelo resíduo ------------------------------------
    rel.subtitulo("Parte 6 - Verificação: r = A_orig c - b_orig")
    r = residuo(A, c, b)                    # A e b originais não foram alterados
    rel.escrever("r = [" + ", ".join(f"{valor:.2e}" for valor in r) + "]")
    rel.escrever(f"||r||_inf = {norma_inf(r):.2e}")

    # ---- FLOPs ---------------------------------------------------------------
    rel.subtitulo("Contagem de FLOPs (n = 5)")
    md_teo, ss_teo = flops_teoricos(n)
    rel.escrever(f"Multiplicações/divisões: programa = {flops['mult_div']},  teoria = {md_teo}")
    rel.escrever(f"Somas/subtrações:        programa = {flops['soma_sub']},  teoria = {ss_teo}")
    rel.escrever(f"Total de FLOPs:          programa = {flops['mult_div'] + flops['soma_sub']},"
                 f"  teoria = {md_teo + ss_teo}")
    rel.escrever(f"Estimativa assintótica 2n^3/3 + n^2 = {2 * n ** 3 / 3 + n ** 2:.1f} "
                 "(só vale bem para n grande)")

    # ---- Parte 7: limitação da versão ingênua ----------------------------------
    rel.subtitulo("Parte 7 - Pivô nulo: mesmo sistema com as equações 1 e 3 trocadas")
    ordem = [2, 1, 0, 3, 4]                 # a nova 1a equação tem coeficiente 0 em c1
    try:
        gauss_ingenua(A[ordem], b[ordem])
    except ZeroDivisionError as erro:
        rel.escrever(f"Programa interrompido: {erro}")

    rel.salvar(SAIDA / "resultados.txt")


if __name__ == "__main__":
    main()
