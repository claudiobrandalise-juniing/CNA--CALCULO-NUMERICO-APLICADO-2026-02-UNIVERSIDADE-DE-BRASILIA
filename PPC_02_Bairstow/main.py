"""
PPC 02 - Método de Bairstow para raízes de polinômios
=====================================================
1. Bairstow para polinômios de grau qualquer (reais e complexas)   -> utils/polinomios.py
2. Validação com polinômio de 7a ordem de raízes conhecidas
3. Influência do chute (r0, s0) na convergência
4. Autovalores do sistema massa-mola-amortecedor de 3 GDL (APC 2):
       det(lambda^2 M + lambda C + K) = 0   (polinômio de grau 6)
5. Interpretação física (estabilidade, amortecimento, frequência)
6-7. Varredura de (r0, s0) e fractal de Bairstow                   -> fractal.py

Executar (a partir da raiz do repositório):
    python PPC_02_Bairstow/main.py
"""

import json
import math
import sys
from pathlib import Path

import numpy as np

# ---- acesso à biblioteca utils/ (raiz do repositório) ------------------------
PASTA = Path(__file__).resolve().parent
SAIDA = PASTA / "outputs"
sys.path.insert(0, str(PASTA.parent))
sys.path.insert(0, str(PASTA))
from fractal import desenhar_fractal  # noqa: E402
from utils.polinomios import (avaliar_polinomio, bairstow,  # noqa: E402
                              determinante_polinomial, polinomio_de_raizes)
from utils.saida import Relatorio  # noqa: E402

with open(PASTA / "inputs" / "parametros.json", encoding="utf-8") as arquivo:
    PARAM = json.load(arquivo)
TOL, MAX_ITER = PARAM["tol"], PARAM["max_iter"]


def erro_entre_conjuntos(encontradas, referencia):
    """Maior distância entre cada raiz de referência e a raiz encontrada mais próxima."""
    return max(min(abs(complex(r) - complex(e)) for e in encontradas) for r in referencia)


def matrizes_cadeia(valores):
    """Matriz tridiagonal de uma cadeia parede-elem1-m1-elem2-m2-elem3-m3 (extremo livre).

    Serve para K (rigidezes) e C (amortecimentos): o elemento i liga a massa i-1
    (ou a parede) à massa i.
    """
    n = len(valores)
    matriz = np.zeros((n, n))
    for i in range(n):
        matriz[i, i] = valores[i] + (valores[i + 1] if i + 1 < n else 0.0)
        if i + 1 < n:
            matriz[i, i + 1] = matriz[i + 1, i] = -valores[i + 1]
    return matriz


def polinomio_caracteristico(M, C, K):
    """Coeficientes (ordem crescente) de det(lambda^2 M + lambda C + K).

    Cada elemento da matriz é o polinômio k_ij + c_ij lambda + m_ij lambda^2,
    guardado como [k_ij, c_ij, m_ij]; o determinante é expandido por Laplace.
    """
    n = M.shape[0]
    matriz = [[[K[i, j], C[i, j], M[i, j]] for j in range(n)] for i in range(n)]
    return determinante_polinomial(matriz)


def parte_validacao(rel):
    """Tarefa 2: polinômio de 7a ordem com raízes conhecidas."""
    rel.subtitulo("Tarefa 2 - Validação com polinômio de 7a ordem")
    raizes_ref = [complex(re, im) for re, im in PARAM["validacao"]["raizes"]]
    coef = polinomio_de_raizes(raizes_ref)        # multiplica os monômios (x - r_i)
    rel.escrever("Raízes conhecidas: " + ", ".join(f"{r:g}" for r in raizes_ref))
    rel.escrever("Coeficientes [a0..a7]: " + ", ".join(f"{c:g}" for c in coef))
    encontradas, hist = bairstow(coef, PARAM["validacao"]["r0"], PARAM["validacao"]["s0"],
                                 tol=TOL, max_iter=MAX_ITER)
    rel.tabela(hist, casas=10)
    rel.escrever("Raízes encontradas: " + ", ".join(f"{complex(r):.10g}" for r in encontradas))
    rel.escrever(f"Maior erro em relação às raízes conhecidas: "
                 f"{erro_entre_conjuntos(encontradas, raizes_ref):.2e}")
    return coef, raizes_ref


def parte_convergencia(rel, coef, raizes_ref):
    """Tarefa 3: efeito do chute inicial (r0, s0)."""
    rel.subtitulo("Tarefa 3 - Influência de (r0, s0) na convergência")
    linhas = []
    for r0, s0 in PARAM["chutes_convergencia"]:
        try:
            encontradas, hist = bairstow(coef, r0, s0, tol=TOL, max_iter=MAX_ITER)
            linhas.append({"r0": r0, "s0": s0, "iter. 1o fator": hist[0]["iteracoes"],
                           "iter. totais": sum(h["iteracoes"] for h in hist),
                           "1o par encontrado": f"{complex(hist[0]['raiz_1']):.4g}, "
                                                f"{complex(hist[0]['raiz_2']):.4g}",
                           "erro máx": erro_entre_conjuntos(encontradas, raizes_ref)})
        except (ArithmeticError, ZeroDivisionError) as erro:
            linhas.append({"r0": r0, "s0": s0, "iter. 1o fator": "-", "iter. totais": "-",
                           "1o par encontrado": f"falhou: {erro}", "erro máx": "-"})
    rel.tabela(linhas, casas=4)


def parte_apc2(rel):
    """Tarefas 4 e 5: autovalores do sistema de 3 GDL e interpretação física."""
    sistema = PARAM["sistema_apc2"]
    M = np.diag(sistema["massas_kg"])
    K = matrizes_cadeia(sistema["rigidezes_N_m"])
    C = matrizes_cadeia(sistema["amortecimentos_Ns_m"])
    rel.subtitulo("Tarefa 4 - Polinômio característico do sistema da APC 2")
    rel.matriz(M, "M [kg]")
    rel.matriz(C, "C [N s/m]")
    rel.matriz(K, "K [N/m]")
    coef = polinomio_caracteristico(M, C, K)
    coef_apc = PARAM["polinomio_apc2_crescente"]
    rel.escrever("P(lambda) = " + " + ".join(f"{coef[i]:.10g} l^{i}" for i in range(len(coef) - 1, -1, -1)))
    rel.escrever(f"Confere com o polinômio deduzido à mão na APC 2? "
                 f"{'sim' if np.allclose(coef, coef_apc, rtol=1e-12) else 'NÃO'}")
    autovalores, hist = bairstow(coef, -1.0, -1.0, tol=TOL, max_iter=MAX_ITER)
    rel.tabela(hist, casas=10)

    # ---- tarefa 5: interpretação -------------------------------------------------
    rel.subtitulo("Tarefa 5 - Interpretação física dos autovalores lambda = sigma +- i omega_d")
    # frequências naturais NÃO amortecidas: det(K - omega^2 M) = 0, polinômio em omega^2
    n = M.shape[0]
    coef_sem_amort = determinante_polinomial([[[K[i, j], -M[i, j]] for j in range(n)]
                                              for i in range(n)])
    omega2, _ = bairstow(coef_sem_amort, -1.0, -1.0, tol=TOL, max_iter=MAX_ITER)
    omegas_sem_amort = sorted(math.sqrt(complex(w).real) for w in omega2)
    pares = sorted([z for z in autovalores if complex(z).imag > 0], key=lambda z: complex(z).imag)
    linhas = []
    for modo, (lam, omega_0) in enumerate(zip(pares, omegas_sem_amort), start=1):
        lam = complex(lam)
        omega_n = abs(lam)
        linhas.append({"modo": modo, "sigma = Re [1/s]": lam.real, "omega_d [rad/s]": lam.imag,
                       "f_d [Hz]": lam.imag / (2 * math.pi), "omega_n [rad/s]": omega_n,
                       "zeta [%]": -100 * lam.real / omega_n, "tau = 1/|sigma| [s]": -1 / lam.real,
                       "omega sem amort. [rad/s]": omega_0,
                       "|P(lambda)|/|a6 lambda^6|": abs(avaliar_polinomio(coef, lam)) / abs(coef[-1] * lam ** 6)})
    rel.tabela(linhas, casas=6)
    estavel = all(complex(z).real < 0 for z in autovalores)
    rel.escrever(f"Todas as partes reais negativas -> sistema ESTÁVEL? {'sim' if estavel else 'não'}")
    soma_raizes = sum(complex(z) for z in autovalores).real
    rel.escrever(f"Vieta: soma das raízes = {soma_raizes:.10f}, -a5/a6 = {-coef[5] / coef[6]:.10f}")
    return coef


def main():
    rel = Relatorio()
    rel.titulo("PPC 02 - Método de Bairstow")
    coef_validacao, raizes_ref = parte_validacao(rel)
    parte_convergencia(rel, coef_validacao, raizes_ref)
    coef_apc2 = parte_apc2(rel)

    # ---- tarefas 6 e 7: fractais --------------------------------------------------------
    rel.subtitulo("Tarefas 6 e 7 - Fractais de Bairstow")
    desenhar_fractal(rel, coef_validacao, raizes_ref, PARAM["fractal"],
                     "Polinômio de validação (7a ordem)", SAIDA / "fractal_validacao")
    # o polinômio da APC 2 tem raízes ~ 10-35 rad/s; com lambda = omega_ref * mu as
    # raízes em mu ficam O(1) e cabem na janela (r0, s0) em [-3, 3]
    omega_ref = PARAM["omega_ref_rad_s"]
    coef_escalado = [coef_apc2[i] * omega_ref ** i for i in range(len(coef_apc2))]
    coef_escalado = [c / coef_escalado[-1] for c in coef_escalado]      # mônico
    raizes_escaladas, _ = bairstow(coef_escalado, -1.0, -1.0, tol=TOL, max_iter=MAX_ITER)
    rel.escrever(f"Polinômio da APC 2 reescalado (mu = lambda/{omega_ref:g}): "
                 + ", ".join(f"{c:.6g}" for c in coef_escalado))
    desenhar_fractal(rel, coef_escalado, raizes_escaladas, PARAM["fractal"],
                     f"Polinômio da APC 2 (λ = {omega_ref:g} μ)", SAIDA / "fractal_apc2")

    rel.salvar(SAIDA / "resultados.txt")


if __name__ == "__main__":
    main()
