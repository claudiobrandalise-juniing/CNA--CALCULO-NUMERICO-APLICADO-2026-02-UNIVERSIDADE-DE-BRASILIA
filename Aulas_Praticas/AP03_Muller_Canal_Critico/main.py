"""
Aula Prática 03 - Método de Müller: profundidade crítica em canal trapezoidal
===========================================================================
Escoamento crítico (Fr = 1):  Q^2 T / (g A^3) = 1, com
    A(H) = bH + zH^2   (área molhada)      T(H) = b + 2zH   (largura da superfície)

Multiplicando por g A^3 e expandindo, a profundidade crítica é raiz de

    f(H) = g z^3 H^6 + 3 g b z^2 H^5 + 3 g b^2 z H^4 + g b^3 H^3 - 2 Q^2 z H - Q^2 b = 0

Executar (a partir da raiz do repositório):
    python Aulas_Praticas/AP03_Muller_Canal_Critico/main.py
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
from utils.polinomios import avaliar_polinomio, bairstow  # noqa: E402
from utils.raizes import muller  # noqa: E402
from utils.saida import Relatorio  # noqa: E402

# ---- DADOS DE ENTRADA (enunciado) --------------------------------------------
B_FUNDO = 4.0        # largura de fundo b [m]
Z_TALUDE = 1.5       # inclinação lateral z:1 [-]
G = 9.81             # gravidade [m/s2]
Q_TESTE = 20.0       # vazão do primeiro teste [m3/s]
H_INICIAIS = (0.5, 1.0, 1.5)   # três aproximações iniciais H0, H1, H2 [m]
TOL = 1e-8           # critério |H3 - H2| < TOL [m]
MAX_ITER = 100
VAZOES = [5, 10, 15, 20, 25, 30, 35, 40]   # varredura da Parte 5 [m3/s]


def coeficientes(Q, b=B_FUNDO, z=Z_TALUDE, g=G):
    """Coeficientes [a0, a1, ..., a6] (ordem crescente) do polinômio crítico."""
    return [-Q ** 2 * b,            # a0
            -2.0 * Q ** 2 * z,      # a1
            0.0,                    # a2 (nulo: não há termo em H^2)
            g * b ** 3,             # a3
            3.0 * g * b ** 2 * z,   # a4
            3.0 * g * b * z ** 2,   # a5
            g * z ** 3]             # a6


def profundidade_critica(Q, H_iniciais=H_INICIAIS):
    """Aplica Müller ao polinômio da vazão Q. Retorna (Hc, historico)."""
    coef = coeficientes(Q)
    return muller(lambda H: avaliar_polinomio(coef, H), *H_iniciais,
                  tol=TOL, max_iter=MAX_ITER, criterio="passo")


def froude_quadrado(H, Q):
    """Fr^2 = Q^2 T / (g A^3): deve dar 1 na profundidade crítica."""
    area = B_FUNDO * H + Z_TALUDE * H ** 2
    largura = B_FUNDO + 2.0 * Z_TALUDE * H
    return Q ** 2 * largura / (G * area ** 3)


def main():
    rel = Relatorio()
    configurar_estilo()
    rel.titulo("AP03 - Profundidade crítica em canal trapezoidal (Müller)")

    # ---- Parte 1: coeficientes do polinômio para Q = 20 m3/s ----------------
    coef = coeficientes(Q_TESTE)
    rel.subtitulo("Parte 1 - Polinômio para Q = 20 m3/s (conferir com o enunciado)")
    rel.escrever("f(H) = " + " + ".join(f"({coef[i]:.10g}) H^{i}" for i in range(6, -1, -1)))
    rel.escrever("Esperado: 33.10875 H^6 + 264.87 H^5 + 706.32 H^4 + 627.84 H^3 - 1200 H - 1600")

    # ---- Partes 2 e 3: Müller --------------------------------------------------
    rel.subtitulo("Partes 2 e 3 - Iterações de Müller (H0, H1, H2) = (0.5, 1.0, 1.5) m")
    Hc, hist = profundidade_critica(Q_TESTE)
    rel.tabela(hist, casas=12)
    rel.escrever(f"\nProfundidade crítica Hc = {Hc:.8f} m em {len(hist)} iterações")

    # ---- Parte 4: verificação física ----------------------------------------
    rel.subtitulo("Parte 4 - Verificação física")
    Fr2 = froude_quadrado(Hc, Q_TESTE)
    rel.escrever(f"A_c = {B_FUNDO * Hc + Z_TALUDE * Hc ** 2:.6f} m2,  "
                 f"T_c = {B_FUNDO + 2 * Z_TALUDE * Hc:.6f} m")
    rel.escrever(f"Fr^2 = {Fr2:.12f}  ->  Fr = {np.sqrt(Fr2):.12f}")
    real_positiva = isinstance(Hc, float) and Hc > 0
    rel.escrever(f"Raiz real e positiva (fisicamente admissível)? {'sim' if real_positiva else 'não'}")

    # ---- todas as raízes (Bairstow) para discutir a questão d ----------------
    rel.subtitulo("Todas as 6 raízes do polinômio (Bairstow, para comparação)")
    todas, _ = bairstow(coef, -1.0, -1.0)
    for raiz in sorted(todas, key=lambda r: (complex(r).real, complex(r).imag)):
        tipo = "real" if isinstance(raiz, float) else "complexa"
        rel.escrever(f"  {raiz:.6f}  ({tipo}{', ADMISSÍVEL' if tipo == 'real' and raiz > 0 else ''})")

    # ---- questão f: outros pontos iniciais ------------------------------------
    rel.subtitulo("Questão f - Outras aproximações iniciais (Q = 20 m3/s)")
    testes = []
    for iniciais in ((0.5, 1.0, 1.5), (0.1, 0.2, 0.3), (2.0, 3.0, 4.0), (-3.0, -2.5, -2.0)):
        H_teste, h_teste = profundidade_critica(Q_TESTE, iniciais)
        testes.append({"H0, H1, H2": str(iniciais), "raiz": H_teste, "iterações": len(h_teste)})
    rel.tabela(testes, casas=10)

    # ---- Parte 5: varredura de vazões -----------------------------------------
    rel.subtitulo("Parte 5 - Profundidade crítica para diferentes vazões")
    resultados = []
    for Q in VAZOES:
        Hc_Q, h_Q = profundidade_critica(Q)
        resultados.append({"Q [m3/s]": Q, "Hc [m]": Hc_Q, "Fr": np.sqrt(froude_quadrado(Hc_Q, Q)),
                           "iterações": len(h_Q)})
    rel.tabela(resultados, casas=8)
    tabela_Q_Hc = np.array([[r["Q [m3/s]"], r["Hc [m]"]] for r in resultados])
    SAIDA.mkdir(parents=True, exist_ok=True)
    np.savetxt(SAIDA / "Hc_vs_Q.txt", tabela_Q_Hc, fmt="%10.6f", header="Q[m3/s]   Hc[m]")

    # ---- gráficos -------------------------------------------------------------
    fig, ax = plt.subplots()
    ax.plot(tabela_Q_Hc[:, 0], tabela_Q_Hc[:, 1], "o-")
    ax.set(xlabel="Q [m³/s]", ylabel="Hc [m]",
           title="Profundidade crítica x vazão (b = 4 m, z = 1,5)")
    salvar_figura(fig, SAIDA / "Hc_vs_Q.png")

    H_graf = np.linspace(0, 2.5, 300)
    fig, ax = plt.subplots()
    ax.plot(H_graf, avaliar_polinomio(coef, H_graf))
    ax.axhline(0, color="k", lw=0.8)
    ax.plot(Hc, 0, "o", label=f"Hc = {Hc:.4f} m")
    ax.set(xlabel="H [m]", ylabel="f(H)", title="Polinômio crítico (Q = 20 m³/s)")
    ax.legend()
    salvar_figura(fig, SAIDA / "polinomio.png")

    rel.salvar(SAIDA / "resultados.txt")


if __name__ == "__main__":
    main()
