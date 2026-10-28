"""
saida.py
========
Rotinas de saída de resultados.

A classe Relatorio imprime cada linha no console e, ao mesmo tempo, guarda tudo
numa lista. No final do programa basta chamar salvar() para gravar exatamente o
que apareceu na tela num arquivo .txt dentro de outputs/. Isso garante que os
resultados do repositório sejam reprodutíveis (exigência do Manual).

Todos os métodos iterativos da biblioteca devolvem um "historico": uma lista de
dicionários, um por iteração, sempre com as mesmas chaves. O método tabela()
transforma esse histórico numa tabela alinhada.
"""

from pathlib import Path

import numpy as np


def formatar_valor(valor, casas=10):
    """Converte um número em texto curto e legível para tabelas.

    Parâmetros
    ----------
    valor : int, float, complex, bool, str ou array
    casas : int
        Algarismos significativos usados para números reais.
    """
    # bool precisa vir antes de int, porque em Python True também é int
    if isinstance(valor, (bool, np.bool_)):
        return "sim" if valor else "não"
    # inteiros (contadores, número da iteração) aparecem sem casas decimais
    if isinstance(valor, (int, np.integer)):
        return f"{int(valor):d}"
    # complexos aparecem como a+bj
    if isinstance(valor, (complex, np.complexfloating)):
        return f"{valor.real:.{casas}g}{valor.imag:+.{casas}g}j"
    # reais: notação compacta com 'casas' algarismos significativos
    if isinstance(valor, (float, np.floating)):
        return f"{float(valor):.{casas}g}"
    # vetores aparecem entre colchetes, separados por vírgula
    if isinstance(valor, (list, tuple, np.ndarray)):
        return "[" + ", ".join(formatar_valor(v, casas) for v in np.ravel(valor)) + "]"
    # qualquer outra coisa (texto, None) é convertida diretamente
    return str(valor)


class Relatorio:
    """Imprime no console e acumula o texto para salvar em arquivo.

    Exemplo
    -------
        rel = Relatorio()
        rel.titulo("Método da bisseção")
        rel.escrever(f"Raiz = {raiz:.6f}")
        rel.tabela(historico)
        rel.salvar("outputs/resultados.txt")
    """

    def __init__(self):
        # lista com todas as linhas já escritas, na ordem em que apareceram
        self.linhas = []

    def escrever(self, texto=""):
        """Imprime uma linha e guarda uma cópia dela."""
        print(texto)
        self.linhas.append(str(texto))

    def titulo(self, texto):
        """Escreve um título destacado entre duas linhas de '='."""
        barra = "=" * 78
        self.escrever("")
        self.escrever(barra)
        self.escrever(texto)
        self.escrever(barra)

    def subtitulo(self, texto):
        """Escreve um subtítulo sublinhado com '-'."""
        self.escrever("")
        self.escrever(texto)
        self.escrever("-" * len(texto))

    def tabela(self, historico, casas=10, colunas=None):
        """Imprime uma lista de dicionários como tabela alinhada.

        Parâmetros
        ----------
        historico : list[dict]
            Uma linha por dicionário; as chaves viram o cabeçalho.
        casas : int
            Algarismos significativos dos números reais.
        colunas : list[str] ou None
            Se informado, mostra apenas essas chaves (nessa ordem).
        """
        # tabela vazia: avisa e sai sem erro
        if not historico:
            self.escrever("(histórico vazio)")
            return
        # se o usuário não escolheu colunas, usa todas as chaves que aparecem
        # (na ordem em que surgem; linhas sem alguma chave ficam em branco)
        if colunas is None:
            colunas = []
            for linha in historico:
                colunas += [c for c in linha if c not in colunas]
        # converte todas as células em texto antes de medir a largura
        celulas = [[formatar_valor(linha.get(c, ""), casas) for c in colunas]
                   for linha in historico]
        # largura de cada coluna = maior texto entre o cabeçalho e as células
        larguras = [max(len(c), *(len(linha[j]) for linha in celulas))
                    for j, c in enumerate(colunas)]
        # cabeçalho, linha separadora e uma linha de texto por iteração
        self.escrever("  ".join(c.rjust(larguras[j]) for j, c in enumerate(colunas)))
        self.escrever("  ".join("-" * larguras[j] for j in range(len(colunas))))
        for linha in celulas:
            self.escrever("  ".join(linha[j].rjust(larguras[j]) for j in range(len(colunas))))

    def matriz(self, M, nome="", casas=6):
        """Imprime uma matriz (ou vetor) com colunas alinhadas."""
        if nome:
            self.escrever(f"{nome} =")
        M = np.atleast_2d(M)
        textos = [[formatar_valor(v, casas) for v in linha] for linha in M]
        largura = max(len(t) for linha in textos for t in linha)
        for linha in textos:
            self.escrever("  [ " + "  ".join(t.rjust(largura) for t in linha) + " ]")

    def salvar(self, caminho):
        """Grava todas as linhas acumuladas num arquivo de texto (UTF-8)."""
        caminho = Path(caminho)
        # cria a pasta outputs/ caso ela ainda não exista
        caminho.parent.mkdir(parents=True, exist_ok=True)
        caminho.write_text("\n".join(self.linhas) + "\n", encoding="utf-8")
        print(f"\n[resultados salvos em {caminho}]")
