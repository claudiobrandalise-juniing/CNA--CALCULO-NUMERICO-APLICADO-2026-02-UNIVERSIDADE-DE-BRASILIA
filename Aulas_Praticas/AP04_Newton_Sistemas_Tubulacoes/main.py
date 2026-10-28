"""
Aula Prática 04 - Newton-Raphson para sistemas não lineares
=============================================================
Divisão de vazão entre dois ramos paralelos de tubulação (A -> B).

    F1(Q1, Q2) = Q1 + Q2 - Q          = 0   (conservação da vazão)
    F2(Q1, Q2) = R1 Q1^2 - R2 Q2^2    = 0   (mesma perda de carga nos dois ramos)

com  R_i = (f_i L_i/D_i + soma K_i) * 8 / (g pi^2 D_i^4)   [s2/m5]

Jacobiana:
    J = [ dF1/dQ1  dF1/dQ2 ] = [   1         1      ]
        [ dF2/dQ1  dF2/dQ2 ]   [ 2 R1 Q1  -2 R2 Q2  ]

Executar (a partir da raiz do repositório):
    python Aulas_Praticas/AP04_Newton_Sistemas_Tubulacoes/main.py
"""

import math
import sys
from pathlib import Path

import numpy as np

# ---- acesso à biblioteca utils/ (raiz do repositório) ------------------------
PASTA = Path(__file__).resolve().parent
SAIDA = PASTA / "outputs"
sys.path.insert(0, str(PASTA.parents[1]))
from utils.saida import Relatorio  # noqa: E402
from utils.sistemas_nao_lineares import newton_sistemas  # noqa: E402

# ---- DADOS DE ENTRADA (enunciado) --------------------------------------------
Q_TOTAL = 0.050      # vazão total [m3/s]
G = 9.81             # gravidade [m/s2]
# ramo 1: comprimento [m], diâmetro [m], fator de atrito [-], perda localizada [-]
L1, D1, F_ATRITO1, K1 = 30.0, 0.10, 0.020, 1.5
# ramo 2 (com válvula parcialmente aberta, Kv)
L2, D2, F_ATRITO2, K2, KV = 20.0, 0.08, 0.022, 1.0, 8.0
Q_INICIAL = [0.030, 0.020]   # chute (Q1, Q2) [m3/s]
TOL = 1e-8           # critério max(|dQ1|, |dQ2|) < TOL [m3/s]
MAX_ITER = 50


def resistencia(f_atrito, L, D, soma_K):
    """Coeficiente R de h_L = R Q^2 [s2/m5] (perdas distribuída + localizadas)."""
    return (f_atrito * L / D + soma_K) * 8.0 / (G * math.pi ** 2 * D ** 4)


R1 = resistencia(F_ATRITO1, L1, D1, K1)
R2 = resistencia(F_ATRITO2, L2, D2, K2 + KV)


def residuos(Q):
    """Vetor F(Q1, Q2) do sistema não linear."""
    Q1, Q2 = Q
    return np.array([Q1 + Q2 - Q_TOTAL,
                     R1 * Q1 ** 2 - R2 * Q2 ** 2])


def jacobiana(Q):
    """Matriz Jacobiana analítica J_ij = dF_i/dQ_j."""
    Q1, Q2 = Q
    return np.array([[1.0, 1.0],
                     [2.0 * R1 * Q1, -2.0 * R2 * Q2]])


def velocidade(Q, D):
    """U = 4 Q / (pi D^2) [m/s]."""
    return 4.0 * Q / (math.pi * D ** 2)


def main():
    rel = Relatorio()
    rel.titulo("AP04 - Newton-Raphson para sistemas: ramos paralelos de tubulação")

    # ---- Parte 1: coeficientes de resistência -----------------------------------
    rel.subtitulo("Parte 1 - Coeficientes de resistência")
    rel.escrever(f"R1 = {R1:.4f} s2/m5   (f L/D + K = {F_ATRITO1 * L1 / D1 + K1:.2f})")
    rel.escrever(f"R2 = {R2:.4f} s2/m5   (f L/D + K + Kv = {F_ATRITO2 * L2 / D2 + K2 + KV:.2f})")

    # ---- Partes 2 e 3: Newton -------------------------------------------------
    rel.subtitulo("Partes 2 e 3 - Iterações de Newton (sistema linear 2x2 por Cramer)")
    Q, hist = newton_sistemas(residuos, jacobiana, Q_INICIAL, tol=TOL, max_iter=MAX_ITER)
    # renomeia as colunas genéricas (x1, F1, dx1...) para a notação do problema
    nomes = {"x1": "Q1", "x2": "Q2", "dx1": "dQ1", "dx2": "dQ2"}
    rel.tabela([{nomes.get(c, c): v for c, v in linha.items()} for linha in hist], casas=10)
    Q1, Q2 = Q
    rel.escrever(f"\nQ1 = {Q1:.8f} m3/s ({100 * Q1 / Q_TOTAL:.2f} %),  "
                 f"Q2 = {Q2:.8f} m3/s ({100 * Q2 / Q_TOTAL:.2f} %)  em {len(hist)} iterações")

    # ---- Parte 4: verificação física -----------------------------------------
    rel.subtitulo("Parte 4 - Verificação física")
    rel.escrever(f"U1 = {velocidade(Q1, D1):.4f} m/s,  U2 = {velocidade(Q2, D2):.4f} m/s")
    rel.escrever(f"Q1 + Q2 = {Q1 + Q2:.10f} m3/s  (Q = {Q_TOTAL})")
    rel.escrever(f"h_L1 = R1 Q1^2 = {R1 * Q1 ** 2:.8f} m,  h_L2 = R2 Q2^2 = {R2 * Q2 ** 2:.8f} m")
    # solução analítica: R1 Q1^2 = R2 Q2^2 -> Q1/Q2 = sqrt(R2/R1)
    Q1_exato = Q_TOTAL * math.sqrt(R2) / (math.sqrt(R1) + math.sqrt(R2))
    rel.escrever(f"Referência analítica Q1 = Q sqrt(R2)/(sqrt(R1)+sqrt(R2)) = {Q1_exato:.10f}"
                 f"  -> erro = {abs(Q1 - Q1_exato):.1e} m3/s")

    # ---- questão f: outros chutes ---------------------------------------------
    rel.subtitulo("Questão f - Outras estimativas iniciais")
    testes = []
    for chute in ([0.030, 0.020], [0.010, 0.040], [0.100, 0.100], [0.001, 0.001], [0.080, -0.030]):
        Q_t, h_t = newton_sistemas(residuos, jacobiana, chute, tol=TOL, max_iter=MAX_ITER)
        testes.append({"chute (Q1, Q2)": str(chute), "Q1": Q_t[0], "Q2": Q_t[1],
                       "iterações": len(h_t)})
    rel.tabela(testes, casas=8)

    rel.salvar(SAIDA / "resultados.txt")


if __name__ == "__main__":
    main()
