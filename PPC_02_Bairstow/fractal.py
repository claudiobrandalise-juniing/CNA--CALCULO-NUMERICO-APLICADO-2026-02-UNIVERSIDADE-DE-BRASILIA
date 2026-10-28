"""
fractal.py (PPC 02)
===================
Varredura sistemática do chute (r0, s0) do método de Bairstow e desenho do
"fractal de Bairstow".

Para cada ponto (r0, s0) de uma grade, roda-se o Bairstow (só o 1o fator
quadrático) e registra-se:
    - o número de iterações até convergir (cor do mapa); preto = divergiu;
    - para qual fator quadrático x^2 - r x - s ele convergiu (bacias de atração).

O cálculo pesado está em utils.polinomios.mapa_convergencia_bairstow, que
aplica a MESMA correção de Bairstow a todos os pontos da grade de uma vez.
"""

import itertools

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

from utils.padrao_grafico import salvar_figura
from utils.polinomios import mapa_convergencia_bairstow


def _texto_raiz(z):
    """Raiz em texto curto: '2' para reais, '1+2j' para complexas."""
    return f"{z.real:.3g}" if abs(z.imag) < 1e-12 else f"{z:.3g}"


def fatores_candidatos(raizes):
    """Todos os fatores quadráticos REAIS x^2 - r x - s formados por 2 raízes.

    Para um par (x1, x2):  r = x1 + x2  e  s = -x1 x2. O par só gera um fator
    real se as duas forem reais ou se forem complexas conjugadas.
    """
    fatores = []
    for x1, x2 in itertools.combinations([complex(z) for z in raizes], 2):
        r, s = x1 + x2, -x1 * x2
        if abs(r.imag) < 1e-8 and abs(s.imag) < 1e-8:
            fatores.append((r.real, s.real, f"{_texto_raiz(x1)} e {_texto_raiz(x2)}"))
    return fatores


def desenhar_fractal(rel, coef, raizes, config, titulo, prefixo_arquivo):
    """Gera o mapa de iterações e o mapa de bacias para um polinômio."""
    r0 = np.linspace(config["r_min"], config["r_max"], config["resolucao"])
    s0 = np.linspace(config["s_min"], config["s_max"], config["resolucao"])
    R0, S0 = np.meshgrid(r0, s0)          # R0 varia nas colunas (eixo x), S0 nas linhas
    iteracoes, r_final, s_final = mapa_convergencia_bairstow(
        coef, R0, S0, tol=config["tol"], max_iter=config["max_iter"])

    # ---- bacias: fator candidato mais próximo do (r, s) final -----------------
    fatores = fatores_candidatos(raizes)
    distancias = np.stack([np.hypot(r_final - r, s_final - s) for r, s, _ in fatores])
    bacia = np.argmin(np.nan_to_num(distancias, nan=np.inf), axis=0)
    convergiu = iteracoes > 0
    extensao = (config["r_min"], config["r_max"], config["s_min"], config["s_max"])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    mapa_cores = plt.get_cmap("viridis").copy()
    mapa_cores.set_bad("black")            # divergências em preto
    # escala logarítmica: a maioria converge rápido e as bordas das bacias demoram
    imagem = ax1.imshow(np.ma.masked_where(~convergiu, iteracoes), origin="lower",
                        extent=extensao, cmap=mapa_cores, aspect="auto",
                        norm=LogNorm(vmin=1, vmax=config["max_iter"]))
    fig.colorbar(imagem, ax=ax1, label="iterações até convergir")
    ax1.set(xlabel="r0", ylabel="s0", title=f"{titulo}\nnúmero de iterações (preto = divergiu)")
    ax1.grid(False)

    paleta = plt.get_cmap("tab20", len(fatores)).copy()
    paleta.set_bad("black")
    imagem = ax2.imshow(np.ma.masked_where(~convergiu, bacia), origin="lower", extent=extensao,
                        cmap=paleta, vmin=-0.5, vmax=len(fatores) - 0.5, aspect="auto")
    usados = sorted(set(bacia[convergiu].tolist()))
    barra = fig.colorbar(imagem, ax=ax2, ticks=usados)
    barra.ax.set_yticklabels([fatores[i][2] for i in usados], fontsize=7)
    ax2.set(xlabel="r0", ylabel="s0", title="bacias de atração (par de raízes encontrado)")
    ax2.grid(False)
    salvar_figura(fig, f"{prefixo_arquivo}.png")

    # ---- resumo numérico ------------------------------------------------------------
    total = iteracoes.size
    rel.escrever(f"{titulo}: grade {config['resolucao']}x{config['resolucao']}, "
                 f"convergiu em {100 * convergiu.sum() / total:.1f} % dos chutes, "
                 f"iterações (mediana) = {np.median(iteracoes[convergiu]):.0f}, "
                 f"máx = {iteracoes[convergiu].max()}")
    for i in usados:
        rel.escrever(f"   bacia do par {fatores[i][2]}: {100 * np.sum(convergiu & (bacia == i)) / total:.1f} %"
                     " da janela")
