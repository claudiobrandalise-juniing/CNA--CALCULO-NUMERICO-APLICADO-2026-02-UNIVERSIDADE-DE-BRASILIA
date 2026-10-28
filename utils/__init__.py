"""
Biblioteca de métodos numéricos da disciplina ENM0227 - Cálculo Numérico Aplicado.

Cada módulo reúne os métodos de um capítulo do curso. Todos foram escritos
"na mão", usando o NumPy apenas para guardar vetores e matrizes (como exige o
Manual de Boas Práticas). Nenhuma função aqui chama um resolvedor pronto.

Módulos:
    raizes                 -> zeros de funções de uma variável (bisseção, falsa
                              posição, falsa posição modificada, ponto fixo,
                              Newton-Raphson, secante, Müller, busca de intervalos)
    polinomios             -> Horner, Briot-Ruffini, Bairstow, operações com
                              polinômios e determinante de matriz de polinômios
    sistemas_nao_lineares  -> Newton-Raphson para sistemas e Jacobiana numérica
    sistemas_lineares      -> métodos diretos (Gauss ingênua, pivoteamento,
                              Gauss-Jordan, Cramer, LU, Crout, inversa, Thomas,
                              Cholesky)
    sistemas_iterativos    -> Jacobi e Gauss-Seidel (com relaxação) para Ax = b
    otimizacao             -> razão áurea, interpolação quadrática, Newton 1D
    derivadas              -> derivadas numéricas por diferenças finitas
    edo                    -> Euler, Heun e Runge-Kutta de 4a ordem (sistemas)
    edp                    -> Laplace/Poisson 2D por Gauss-Seidel (Liebmann) e
                              condução transiente 1D (explícito, implícito,
                              Crank-Nicolson) com algoritmo de Thomas
    metricas_erro          -> erros, RMSE, normas e número de condição
    padrao_grafico         -> estilo único dos gráficos do repositório
    saida                  -> impressão de tabelas e gravação de resultados .txt

Como importar a partir de uma pasta qualquer do repositório:

    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # sobe até a raiz
    from utils.raizes import bissecao
"""
