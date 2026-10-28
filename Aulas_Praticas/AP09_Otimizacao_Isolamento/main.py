"""
Aula Prática 09 - Otimização unidimensional: espessura ótima de isolamento
==========================================================================
Tubulação de raio r1 isolada com espessura e (raio externo r2 = r1 + e):

    R'_cond(e) = ln((r1 + e)/r1) / (2 pi k)        [m K / W]
    R'_conv(e) = 1 / (2 pi h (r1 + e))              [m K / W]
    q'(e)      = (Ti - T_inf) / (R'_cond + R'_conv)  [W/m]

Custos no período de operação t:
    C_E(e) = pE * L q'(e) t / 1000                   (energia perdida, R$)
    C_I(e) = cI * pi * L * ((r1 + e)^2 - r1^2)       (isolante instalado, R$)
    C(e)   = C_E(e) + C_I(e)     -> minimizar em 0 <= e <= 0.12 m

Executar (a partir da raiz do repositório):
    python Aulas_Praticas/AP09_Otimizacao_Isolamento/main.py
"""

import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# ---- acesso à biblioteca utils/ (raiz do repositório) ------------------------
PASTA = Path(__file__).resolve().parent
SAIDA = PASTA / "outputs"
sys.path.insert(0, str(PASTA.parents[1]))
from utils.otimizacao import FuncaoContada, interpolacao_quadratica, razao_aurea  # noqa: E402
from utils.padrao_grafico import configurar_estilo, salvar_figura  # noqa: E402
from utils.saida import Relatorio  # noqa: E402

# ---- DADOS DE ENTRADA (Tabela 1 do enunciado) ----------------------------------
R1 = 0.040          # raio externo da tubulação sem isolamento [m]
L_TUBO = 20.0       # comprimento [m]
T_I = 180.0         # temperatura na interface tubo-isolante [°C]
T_INF = 40.0        # temperatura ambiente [°C]
K_ISOL = 0.045      # condutividade do isolante [W/(m K)]
H_EXT = 12.0        # convecção externa [W/(m2 K)]
TEMPO = 30000.0     # tempo acumulado de operação [h]
P_ENERGIA = 0.70    # preço da energia térmica [R$/kWh]
C_ISOLANTE = 9000.0  # custo instalado do isolante [R$/m3]
E_MIN, E_MAX = 0.0, 0.12    # intervalo admissível da espessura [m]
TOL = 1e-6          # tolerância em espessura [m]
MAX_ITER = 200
PONTOS_QUADRATICA = (0.04, 0.07, 0.11)    # tripla inicial da interpolação [m]
PRECOS = [0.35, 0.70, 1.05, 1.40]          # estudo paramétrico [R$/kWh]


def perda_por_metro(e):
    """q'(e) [W/m]: resistências de condução (isolante) e convecção em série."""
    r2 = R1 + e
    resistencia_cond = math.log(r2 / R1) / (2.0 * math.pi * K_ISOL)
    resistencia_conv = 1.0 / (2.0 * math.pi * H_EXT * r2)
    return (T_I - T_INF) / (resistencia_cond + resistencia_conv)


def custo_energia(e, preco=P_ENERGIA):
    """C_E(e) [R$]: (L q' [W]) * t [h] / 1000 = energia em kWh, vezes o preço."""
    return preco * L_TUBO * perda_por_metro(e) * TEMPO / 1000.0


def custo_isolante(e):
    """C_I(e) [R$]: volume do anel de isolante vezes o custo por m3."""
    return C_ISOLANTE * math.pi * L_TUBO * ((R1 + e) ** 2 - R1 ** 2)


def custo_total(e, preco=P_ENERGIA):
    """Função objetivo C(e) = C_E + C_I [R$]."""
    return custo_energia(e, preco) + custo_isolante(e)


def otimizar(rel, preco):
    """Roda razão áurea e interpolação quadrática para um preço de energia."""
    resultados = {}
    # --- razão áurea -------------------------------------------------------------
    C_aurea = FuncaoContada(lambda e: custo_total(e, preco))
    e_aurea, hist_aurea = razao_aurea(C_aurea, E_MIN, E_MAX, tol=TOL, max_iter=MAX_ITER)
    resultados["Razão áurea"] = (e_aurea, len(hist_aurea), C_aurea.n_avaliacoes,
                                 "|b - a| < tol" if hist_aurea[-1]["largura"] < TOL else "max_iter")
    # --- interpolação quadrática ------------------------------------------------
    C_quad = FuncaoContada(lambda e: custo_total(e, preco))
    e_quad, hist_quad = interpolacao_quadratica(C_quad, *PONTOS_QUADRATICA, tol=TOL,
                                                max_iter=MAX_ITER)
    protecoes = sum(1 for h in hist_quad if h["passo"] != "parabólico")
    resultados["Interp. quadrática"] = (e_quad, len(hist_quad), C_quad.n_avaliacoes,
                                        f"x3 - x1 < tol ({protecoes} passos de proteção)")
    rel.subtitulo("Parte 2 - Razão áurea (últimas iterações)")
    rel.tabela(hist_aurea[-5:], casas=10)
    rel.subtitulo("Parte 3 - Interpolação quadrática (todas as iterações)")
    rel.tabela(hist_quad, casas=10)
    return resultados


def main():
    rel = Relatorio()
    configurar_estilo()
    rel.titulo("AP09 - Espessura ótima de isolamento térmico")

    # ---- Parte 1: exploração da função objetivo ---------------------------------
    rel.subtitulo("Parte 1 - Exploração de C(e)")
    rel.tabela([{"e [m]": e, "q' [W/m]": perda_por_metro(e), "C_E [R$]": custo_energia(e),
                 "C_I [R$]": custo_isolante(e), "C [R$]": custo_total(e)}
                for e in np.arange(0.0, 0.1201, 0.02)], casas=8)
    f1, f2, f3 = (custo_total(e) for e in PONTOS_QUADRATICA)
    rel.escrever(f"\nTripla da interpolação {PONTOS_QUADRATICA}: C = ({f1:.2f}, {f2:.2f}, {f3:.2f})"
                 f" -> f2 < f1 e f2 < f3? {'sim' if f2 < f1 and f2 < f3 else 'NÃO'}")

    # ---- Partes 2 e 3: os dois métodos ----------------------------------------------
    resultados = otimizar(rel, P_ENERGIA)

    # ---- Parte 3, item 9: comparação -------------------------------------------------
    rel.subtitulo("Item 9 - Comparação dos métodos (mesma tolerância em e)")
    rel.tabela([{"resultado": nome, "e ótima [m]": e, "C mínimo [R$]": custo_total(e),
                 "iterações": it, "avaliações de C": av, "parada": situacao}
                for nome, (e, it, av, situacao) in resultados.items()], casas=10)

    # ---- sensibilidade aos pontos iniciais da interpolação (questão d) ---------------
    rel.subtitulo("Questão d - Interpolação quadrática com outras triplas iniciais")
    linhas = []
    for tripla in ((0.04, 0.07, 0.11), (0.01, 0.05, 0.12), (0.03, 0.06, 0.09)):
        C_teste = FuncaoContada(custo_total)
        e_t, h_t = interpolacao_quadratica(C_teste, *tripla, tol=TOL, max_iter=MAX_ITER)
        linhas.append({"tripla [m]": str(tripla), "e ótima [m]": e_t, "iterações": len(h_t),
                       "avaliações": C_teste.n_avaliacoes})
    rel.tabela(linhas, casas=10)

    # ---- Parte 4: estudo paramétrico do preço da energia -----------------------------
    rel.subtitulo("Parte 4 - Efeito do preço da energia (razão áurea em [0, 0.12] m)")
    parametrico = []
    for preco in PRECOS:
        # a razão áurea só precisa do intervalo admissível (não depende de uma
        # tripla que cerque o mínimo, que muda de lugar com o preço)
        e_otima, _ = razao_aurea(lambda e: custo_total(e, preco), E_MIN, E_MAX,
                                 tol=TOL, max_iter=MAX_ITER)
        parametrico.append({"pE [R$/kWh]": preco, "e ótima [m]": e_otima,
                            "C_E [R$]": custo_energia(e_otima, preco),
                            "C_I [R$]": custo_isolante(e_otima),
                            "C total [R$]": custo_total(e_otima, preco)})
    rel.tabela(parametrico, casas=8)

    # ---- Parte 5: raio crítico ---------------------------------------------------------
    rel.subtitulo("Parte 5 - Raio crítico do isolamento cilíndrico")
    raio_critico = K_ISOL / H_EXT
    rel.escrever(f"r_c = k/h = {raio_critico * 1000:.2f} mm   (r1 = {R1 * 1000:.0f} mm)")
    rel.escrever(f"r1 > r_c -> qualquer isolante REDUZ a perda: q'(0) = {perda_por_metro(0):.2f} W/m,"
                 f" q'(0.01) = {perda_por_metro(0.01):.2f} W/m")

    # ---- gráficos -------------------------------------------------------------------------
    e_graf = np.linspace(E_MIN, E_MAX, 300)
    e_otima = resultados["Razão áurea"][0]
    fig, ax = plt.subplots()
    ax.plot(e_graf * 1000, [custo_energia(e) for e in e_graf], label="C_E (energia)")
    ax.plot(e_graf * 1000, [custo_isolante(e) for e in e_graf], label="C_I (isolante)")
    ax.plot(e_graf * 1000, [custo_total(e) for e in e_graf], "k", label="C total")
    ax.plot(e_otima * 1000, custo_total(e_otima), "o", color="C3",
            label=f"ótimo: e = {e_otima * 1000:.1f} mm")
    ax.set(xlabel="espessura e [mm]", ylabel="custo [R$]", title="Compromisso energia x investimento")
    ax.legend()
    salvar_figura(fig, SAIDA / "custos.png")

    fig, ax = plt.subplots()
    ax.plot([p["pE [R$/kWh]"] for p in parametrico], [p["e ótima [m]"] * 1000 for p in parametrico],
            "o-")
    ax.set(xlabel="preço da energia pE [R$/kWh]", ylabel="espessura ótima [mm]",
           title="Espessura ótima x preço da energia")
    salvar_figura(fig, SAIDA / "e_otima_vs_preco.png")

    rel.salvar(SAIDA / "resultados.txt")


if __name__ == "__main__":
    main()
