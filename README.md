# CNA--CALCULO-NUMERICO-APLICADO-2026-02-UNIVERSIDADE-DE-BRASILIA
Aluno: Claudio Junior Franco Brandalise escrito na matricula 211064816, professor Rafael Gabler Gontijo

# Cálculo Numérico Aplicado — CNA

Repositório de práticas computacionais da disciplina de Cálculo Numérico Aplicado: uma apostila digital dos métodos numéricos implementados ao longo do semestre.

> **Para a prova:** comece pelo [GUIA_DA_PROVA.md](GUIA_DA_PROVA.md), que tem o mapa problema → método → função, as fórmulas e as armadilhas, e pelos [exemplos de temperatura em processadores](Exemplos_Processador/README.md).

---

## 1. Identificação

| Campo | Descrição |
|---|---|
| **Instituição** | Universidade de Brasília (UnB) |
| **Unidade** | Faculdade de Tecnologia (FT) |
| **Departamento** | Departamento de Engenharia Mecânica|
| **Disciplina** | Cálculo Numérico Aplicado — ENM0227 |
| **Semestre** | 2026/2 |
| **Docente** | Prof. Rafael Gabler Gontijo |
| **Discente** | Claudio Junior Franco Brandalise |
| **Matrícula** | 211064816 |

---

## 2. Apresentação

Boa parte dos problemas de engenharia não admite solução analítica fechada. Escoamentos
governados por equações não lineares, sistemas estruturais com muitos graus de liberdade e
fenômenos transientes de transferência de calor só se tornam tratáveis quando reformulados
em termos discretos e resolvidos por aproximações sucessivas. O Cálculo Numérico é a
disciplina que fornece esse ferramental, e mais importante, que fornece os critérios para
julgar quando a aproximação obtida é confiável.

A motivação central deste repositório é que implementar um método numérico é diferente de
utilizá-lo. Ao escrever o algoritmo desde a formulação matemática até a convergência, o
comportamento do erro deixa de ser um resultado tabelado e passa a ser uma consequência
observável de escolhas de discretização, de condicionamento e de aritmética de ponto
flutuante. É essa compreensão que separa o engenheiro que interpreta um resultado de
simulação daquele que apenas o reporta.

Cada prática reunida aqui documenta um método, sua implementação, os dados de entrada e
saída correspondentes e a validação contra uma solução de referência.

---

## 3. Topologia do Repositório

```
.
├── README.md                    # Este arquivo: metadados, topologia e bibliografia geral
├── GUIA_DA_PROVA.md             # Consulta rápida: qual método usar, fórmulas, armadilhas
├── requirements.txt             # Dependências (numpy, matplotlib)
├── executar_tudo.py             # Roda testes + todos os programas e regenera os outputs/
│
├── utils/                       # BIBLIOTECA de métodos (usada por todas as pastas)
│   ├── raizes.py                #   bisseção, falsa posição (+modificada), ponto fixo,
│   │                            #   Newton-Raphson, secante, Müller, varredura de sinais
│   ├── polinomios.py            #   Horner, Briot-Ruffini, Bairstow (+ mapa/fractal),
│   │                            #   produto e determinante de polinômios
│   ├── sistemas_lineares.py     #   Gauss ingênua/pivoteamento, Gauss-Jordan, Cramer,
│   │                            #   LU, Crout, inversa, Thomas, Cholesky, contagem de FLOPs
│   ├── sistemas_iterativos.py   #   Jacobi, Gauss-Seidel com relaxação, dominância diagonal
│   ├── sistemas_nao_lineares.py #   Newton para sistemas, Jacobiana numérica
│   ├── otimizacao.py            #   razão áurea, interpolação quadrática, Newton 1D
│   ├── derivadas.py             #   diferenças finitas progressiva e central
│   ├── edo.py                   #   Euler, Heun, Runge-Kutta 4 (sistemas)
│   ├── edp.py                   #   Laplace/Poisson 2D (Liebmann/SOR); calor 1D transiente
│   ├── metricas_erro.py         #   erros, RMSE, normas, número de condição
│   ├── padrao_grafico.py        #   estilo padronizado dos gráficos
│   └── saida.py                 #   tabelas de iterações e gravação de resultados
│
├── testes/testar_tudo.py        # Valida cada método da biblioteca (50 testes)
│
├── APC_01/  APC_02/  APC_03/    # Atividades para casa: deduções e conferências
├── PPC_01_RungeKutta_Sedimentacao/   # RK4: esfera sedimentando (Stokes x Oseen)
├── PPC_02_Bairstow/                  # Bairstow, autovalores de 3 GDL, fractal
├── PPC_03_Conducao_Transiente/       # FD implícito + Thomas: combustível nuclear
│   ├── README.md                #   resumo, dicionário, I/O, execução, validação
│   ├── main.py                  #   programa principal
│   ├── inputs/parametros.json   #   dados de entrada
│   └── outputs/                 #   resultados.txt, gráficos .png, animação .gif
│
├── Aulas_Praticas/              # 9 atividades feitas em sala (AP01 ... AP09)
└── Exemplos_Processador/        # 8 problemas de temperatura em processadores
                                 # no formato de prova + modelo_em_branco.py
```

**Organização adotada.** Cada prática ocupa um diretório isolado e autocontido, de modo que
seu conteúdo possa ser executado e auditado independentemente das demais. Dentro de cada
diretório, os dados de entrada e os artefatos gerados ficam separados em `inputs/` e
`outputs/`, o que torna explícito o que alimenta o modelo e o que dele resulta. O `README.md`
interno de cada prática descreve o problema modelado, o dicionário de variáveis, as
dependências, a especificação de I/O, os procedimentos de execução e a validação
metodológica.

O diretório `utils/` concentra as rotinas reaproveitadas por mais de uma prática, evitando a
duplicação de lógica entre os scripts (princípio DRY do Manual de Boas Práticas). Todos os
métodos foram escritos a partir das equações. O NumPy é usado apenas para guardar vetores e
matrizes, e nenhum resolvedor pronto (`numpy.linalg.solve`, `numpy.roots`, `scipy`) é chamado.

### Conteúdo por atividade

| Atividade | Método | Problema |
|---|---|---|
| [PPC 01](PPC_01_RungeKutta_Sedimentacao) / [APC 01](APC_01) | Runge-Kutta de 4ª ordem | sedimentação de esfera em baixo Reynolds |
| [PPC 02](PPC_02_Bairstow) / [APC 02](APC_02) | Bairstow | autovalores de sistema massa-mola-amortecedor; fractal |
| [PPC 03](PPC_03_Conducao_Transiente) / [APC 03](APC_03) | diferenças finitas implícitas + Thomas | condução transiente com geração (combustível nuclear) |
| [AP01](Aulas_Praticas/AP01_Bissecao_FalsaPosicao) | bisseção e falsa posição | profundidade em canal (Manning) |
| [AP02](Aulas_Praticas/AP02_Newton_Secante) | Newton-Raphson e secante | temperatura de equilíbrio (convecção + radiação) |
| [AP03](Aulas_Praticas/AP03_Muller_Canal_Critico) | Müller | profundidade crítica em canal trapezoidal |
| [AP04](Aulas_Praticas/AP04_Newton_Sistemas_Tubulacoes) | Newton para sistemas | ramos paralelos de tubulação |
| [AP05](Aulas_Praticas/AP05_Gauss_Ingenua_Reatores) | Gauss ingênua + FLOPs | rede de reatores |
| [AP06](Aulas_Praticas/AP06_Pivoteamento_Barras) | pivoteamento e condicionamento | treliça de barras |
| [AP07](Aulas_Praticas/AP07_Inversao_LU_Molas) | inversão por LU | matriz de flexibilidade |
| [AP08](Aulas_Praticas/AP08_Gauss_Seidel_Placa) | Gauss-Seidel | placa 2D |
| [AP09](Aulas_Praticas/AP09_Otimizacao_Isolamento) | razão áurea e interpolação quadrática | isolamento térmico ótimo |
| [Exemplos](Exemplos_Processador) | todos | temperatura em processadores |

---

## 4. Dependências e execução

- **Python 3.10+**, `numpy` e `matplotlib` (`pip install -r requirements.txt`).
- Cada programa roda a partir da raiz do repositório, por exemplo `python PPC_02_Bairstow/main.py`.
- Para conferir tudo de uma vez (testes + todos os programas, ~1 min): `python executar_tudo.py`.
- Para validar só a biblioteca: `python testes/testar_tudo.py`. Ela compara cada método com soluções conhecidas ou fabricadas, sem usar resolvedores prontos.

---

## 5. Bibliografia Geral

RUGGIERO, M. A. G.; LOPES, V. L. R. **Cálculo numérico: aspectos teóricos e computacionais**.
2. ed. São Paulo: Pearson Makron Books, 1996.

CHAPRA, S. C.; CANALE, R. P. **Métodos numéricos para engenharia**. 7. ed. Porto Alegre:
AMGH, 2016.

BURDEN, R. L.; FAIRES, J. D. **Análise numérica**. São Paulo: Cengage Learning, 2011.

GONTIJO, R. G. **Notas de aula do curso de Cálculo Numérico Aplicado**. Brasília: UnB, 2026.

INCROPERA, F. P.; DEWITT, D. P.; BERGMAN, T. L.; LAVINE, A. S. **Fundamentos de transferência
de calor e de massa**. 7. ed. Rio de Janeiro: LTC, 2014.

HARRIS, C. R. et al. Array programming with NumPy. **Nature**, v. 585, p. 357–362, 2020.

HUNTER, J. D. Matplotlib: a 2D graphics environment. **Computing in Science & Engineering**,
v. 9, n. 3, p. 90–95, 2007.
