"""
MODELO EM BRANCO PARA A PROVA
=============================
Copie este arquivo para a raiz do repositório (ou para esta pasta), renomeie
e preencha as 5 seções. Ele já roda com um problema-exemplo trivial
(processador com convecção pura, que tem solução exata) só para mostrar
o formato.

    1. DADOS        -> todas as constantes com unidade no comentário
    2. MODELO       -> uma função por equação física (balanço de energia, etc.)
    3. SOLUÇÃO      -> chamar o método da biblioteca utils/
    4. VERIFICAÇÃO  -> resíduo ~ 0, balanço de energia, caso limite com
                       solução exata, comparação entre dois métodos
    5. RESULTADOS   -> tabela (Relatorio) + gráfico salvo em arquivo

CHECKLIST RÁPIDO (veja também GUIA_DA_PROVA.md na raiz):
    [ ] arrays com dtype=float (np.array([...], dtype=float))
    [ ] temperaturas de radiação em KELVIN
    [ ] índices do Python começam em 0 (nó 1 da apostila = índice 0)
    [ ] critério de parada E limite de iterações
    [ ] imprimir o número de iterações e o erro final
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# ---- acesso à biblioteca utils/ ------------------------------------------------------------
# se este arquivo estiver na RAIZ do repositório, troque .parent.parent por .parent
PASTA = Path(__file__).resolve().parent
sys.path.insert(0, str(PASTA.parent))
from utils.padrao_grafico import configurar_estilo, salvar_figura  # noqa: E402
from utils.raizes import bissecao, newton_raphson  # noqa: E402
from utils.saida import Relatorio  # noqa: E402

# =============================================================================================
# 1. DADOS
# =============================================================================================
POTENCIA = 65.0          # potência dissipada [W]
H_A = 2.0                # condutância convectiva h*A [W/K]
T_INF = 25.0             # temperatura do ar [°C]
TOL = 1e-8               # tolerância do método


# =============================================================================================
# 2. MODELO
# =============================================================================================
def residuo(T):
    """Balanço de energia: calor que sai - calor gerado [W] (zero no equilíbrio)."""
    return H_A * (T - T_INF) - POTENCIA


def derivada(T):
    """d(residuo)/dT [W/K]."""
    return H_A


# =============================================================================================
# 3, 4 e 5. SOLUÇÃO, VERIFICAÇÃO E RESULTADOS
# =============================================================================================
def main():
    rel = Relatorio()
    configurar_estilo()
    rel.titulo("PROVA - <título do problema>")

    # ---- 3. solução (dois métodos para poder comparar) --------------------------------------
    T_bis, hist_bis = bissecao(residuo, T_INF, T_INF + 200.0, tol=TOL, criterio="funcao")
    T_nr, hist_nr = newton_raphson(residuo, derivada, T_INF + 10.0, tol=TOL)
    rel.tabela(hist_nr, casas=10)

    # ---- 4. verificação -----------------------------------------------------------------------
    T_exata = T_INF + POTENCIA / H_A                 # caso com solução analítica
    rel.escrever(f"Bisseção: T = {T_bis:.8f} °C ({len(hist_bis)} it.) | "
                 f"Newton: T = {T_nr:.8f} °C ({len(hist_nr)} it.)")
    rel.escrever(f"Exata: T = {T_exata:.8f} °C | resíduo f(T) = {residuo(T_nr):.2e} W")

    # ---- 5. resultados --------------------------------------------------------------------------
    T_graf = np.linspace(T_INF, T_INF + 100.0, 200)
    fig, ax = plt.subplots()
    ax.plot(T_graf, residuo(T_graf))
    ax.axhline(0.0, color="k", lw=0.8)
    ax.plot(T_nr, 0.0, "o", label=f"T = {T_nr:.3f} °C")
    ax.set(xlabel="T [°C]", ylabel="f(T) [W]", title="Resíduo do balanço de energia")
    ax.legend()
    salvar_figura(fig, PASTA / "outputs" / "modelo_grafico.png")
    rel.salvar(PASTA / "outputs" / "modelo_resultados.txt")


if __name__ == "__main__":
    main()
