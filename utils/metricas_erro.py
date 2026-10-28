"""
metricas_erro.py
================
Funções globais de erro usadas para validar as soluções numéricas.

Definições (Chapra & Canale, cap. 3):
    erro verdadeiro          E_t   = valor_exato - aproximacao
    erro relativo verdadeiro eps_t = |E_t / valor_exato| * 100 %
    erro relativo aproximado eps_a = |(atual - anterior) / atual| * 100 %
"""

import numpy as np


def erro_absoluto(aproximacao, referencia):
    """|aproximacao - referencia| (funciona para escalares e vetores)."""
    return np.abs(np.asarray(aproximacao, dtype=float) - np.asarray(referencia, dtype=float))


def erro_relativo_verdadeiro(aproximacao, exato):
    """Erro relativo verdadeiro em %: |(exato - aprox)/exato| * 100."""
    return np.abs((exato - aproximacao) / exato) * 100.0


def erro_relativo_aproximado(atual, anterior):
    """Erro relativo aproximado em %: |(atual - anterior)/atual| * 100.

    É o critério de parada clássico do Chapra quando não se conhece a resposta
    exata. Se 'atual' for zero, devolve o erro absoluto para não dividir por 0.
    """
    if atual == 0:
        return abs(atual - anterior)
    return abs((atual - anterior) / atual) * 100.0


def rmse(aproximacao, referencia):
    """Raiz do erro quadrático médio entre dois vetores do mesmo tamanho."""
    diferenca = np.asarray(aproximacao, dtype=float) - np.asarray(referencia, dtype=float)
    return float(np.sqrt(np.mean(diferenca ** 2)))


def norma_l2(vetor):
    """Norma euclidiana ||v||_2 = sqrt(soma v_i^2)."""
    v = np.asarray(vetor, dtype=float)
    return float(np.sqrt(np.sum(v ** 2)))


def norma_inf(vetor):
    """Norma infinito de um vetor: ||v||_inf = max |v_i|."""
    return float(np.max(np.abs(np.asarray(vetor, dtype=float))))


def norma_inf_matriz(A):
    """Norma infinito de uma matriz: maior soma absoluta de linha."""
    A = np.asarray(A, dtype=float)
    return float(np.max(np.sum(np.abs(A), axis=1)))


def numero_condicao_inf(A, A_inversa):
    """Número de condição cond(A) = ||A||_inf * ||A^-1||_inf.

    Interpretação: um erro relativo de 1 parte em 10^p nos dados pode virar um
    erro relativo cond(A) vezes maior na solução. cond(A) ~ 10^q significa que
    se perdem cerca de q algarismos significativos.
    """
    return norma_inf_matriz(A) * norma_inf_matriz(A_inversa)
