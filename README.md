# CNA--CALCULO-NUMERICO-APLICADO-2026-02-UNIVERSIDADE-DE-BRASILIA
Aluno: Claudio Junior Franco Brandalise escrito na matricula 211064816, professor Rafael Gabler Gontijo

# Cálculo Numérico Aplicado — CNA

Repositório de práticas computacionais da disciplina de Cálculo Numérico Aplicado.

---

## 1. Identificação

| Campo | Descrição |
|---|---|
| **Instituição** | Universidade de Brasília (UnB) |
| **Unidade** | Faculdade de Tecnologia (FT) |
| **Departamento** | Departamento de Engenharia Mecânica (ENM) |
| **Disciplina** | Cálculo Numérico Aplicado — *[código da disciplina]* |
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
├── README.md                      # Este arquivo: metadados, topologia e bibliografia geral
│
├── PPC_01/                        # Prática 01 — [tema da prática]
│   ├── README.md                  # Resumo, dicionário de variáveis, I/O e validação
│   ├── main.py                    # Script principal: implementação do algoritmo
│   ├── inputs/
│   │   └── matriz_coeficientes.csv    # Dados de entrada e parâmetros iniciais
│   └── outputs/
│       ├── resultados_analise.txt     # Saídas numéricas e logs de processamento
│       └── erro_convergencia.png      # Gráfico comparativo com a solução de referência
│
├── PPC_02/                        # Prática 02 — [tema da prática]
│   ├── README.md                  # Especificações exclusivas da PPC 02
│   ├── solver_integracao.py       # Script principal da segunda prática
│   ├── inputs/
│   │   └── parametros_iniciais.json   # Condições de contorno e PVI
│   └── outputs/
│       ├── serie_temporal.csv         # Dados exportados após o processamento
│       └── retrato_fase.pdf           # Visualização gráfica padronizada
│
└── utils/                         # Módulos compartilhados entre as práticas
    ├── metricas_erro.py           # Funções globais de erro (RMSE, norma L2)
    └── padrao_grafico.py          # Formatação padronizada de plotagem
```

**Organização adotada.** Cada prática ocupa um diretório isolado e autocontido, de modo que
seu conteúdo possa ser executado e auditado independentemente das demais. Dentro de cada
diretório, os dados de entrada e os artefatos gerados ficam separados em `inputs/` e
`outputs/`, o que torna explícito o que alimenta o modelo e o que dele resulta. O `README.md`
interno de cada prática descreve o problema modelado, o dicionário de variáveis, as
dependências, a especificação de I/O, os procedimentos de execução e a validação
metodológica.

O diretório `utils/` concentra rotinas reaproveitadas por mais de uma prática, evitando a
duplicação de lógica entre os scripts.

---

## 4. Bibliografia Geral

RUGGIERO, M. A. G.; LOPES, V. L. R. **Cálculo numérico: aspectos teóricos e computacionais**.
2. ed. São Paulo: Pearson Makron Books, 1996.

CHAPRA, S. C.; CANALE, R. P. **Métodos numéricos para engenharia**. 7. ed. Porto Alegre:
AMGH, 2016.

BURDEN, R. L.; FAIRES, J. D. **Análise numérica**. São Paulo: Cengage Learning, 2011.

HARRIS, C. R. et al. Array programming with NumPy. **Nature**, v. 585, p. 357–362, 2020.

HUNTER, J. D. Matplotlib: a 2D graphics environment. **Computing in Science & Engineering**,
v. 9, n. 3, p. 90–95, 2007.