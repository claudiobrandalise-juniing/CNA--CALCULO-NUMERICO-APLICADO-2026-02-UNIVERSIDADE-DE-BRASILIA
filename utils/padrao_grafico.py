"""
padrao_grafico.py
=================
Formatação padronizada dos gráficos do repositório (Matplotlib).

Uso típico:
    from utils.padrao_grafico import configurar_estilo, salvar_figura
    configurar_estilo()
    fig, ax = plt.subplots()
    ...
    salvar_figura(fig, "outputs/meu_grafico.png")

Os scripts salvam as figuras em arquivo em vez de abrir janelas, para que
qualquer pessoa consiga reproduzir as imagens só rodando o programa.
"""

from pathlib import Path

import matplotlib.pyplot as plt


def configurar_estilo():
    """Aplica um estilo único (fontes, grade, tamanho) a todas as figuras."""
    plt.rcParams.update({
        "figure.figsize": (7.0, 4.5),   # tamanho padrão em polegadas
        "figure.dpi": 110,              # resolução na tela
        "savefig.dpi": 150,             # resolução do arquivo salvo
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "axes.grid": True,              # grade de fundo ligada
        "grid.alpha": 0.3,
        "axes.spines.top": False,       # remove bordas de cima e da direita
        "axes.spines.right": False,
        "legend.frameon": False,
        "lines.linewidth": 1.6,
    })


def salvar_figura(fig, caminho):
    """Salva a figura (criando a pasta, se preciso) e libera a memória."""
    caminho = Path(caminho)
    caminho.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(caminho, bbox_inches="tight")
    plt.close(fig)
    print(f"[figura salva em {caminho}]")
