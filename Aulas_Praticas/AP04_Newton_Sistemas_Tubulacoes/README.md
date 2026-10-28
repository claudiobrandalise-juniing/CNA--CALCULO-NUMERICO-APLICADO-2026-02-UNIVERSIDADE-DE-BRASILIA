# AP 04 — Newton-Raphson para sistemas não lineares

**Problema:** quanto da vazão total $Q$ passa por cada um de dois ramos paralelos de tubulação entre os pontos A e B.

## Resumo operacional

Duas condições físicas precisam valer ao mesmo tempo:

1. **Conservação da vazão** na bifurcação: $Q_1 + Q_2 = Q$.
2. **Mesma perda de carga** pelos dois caminhos, já que ambos ligam A a B: $h_{L,1} = h_{L,2}$.

A perda de carga de cada ramo usa Darcy-Weisbach para a parcela distribuída e coeficientes $K$ para as perdas localizadas. Com $U = 4Q/(\pi D^2)$, ela fica:

$$h_{L,i} = \left(f_i\frac{L_i}{D_i} + \sum K_i\right)\frac{U_i^2}{2g} = R_i Q_i^2, \qquad R_i = \left(f_i\frac{L_i}{D_i} + \sum K_i\right)\frac{8}{g\pi^2 D_i^4}$$

Isso gera o sistema não linear:

$$F_1 = Q_1 + Q_2 - Q = 0, \qquad F_2 = R_1 Q_1^2 - R_2 Q_2^2 = 0$$

## Algoritmo (Newton-Raphson multivariável)

Em cada iteração, o sistema é linearizado em torno de $(Q_1^{(k)}, Q_2^{(k)})$:

$$\underbrace{\begin{bmatrix} 1 & 1 \\ 2R_1Q_1 & -2R_2Q_2 \end{bmatrix}}_{J}\begin{bmatrix}\Delta Q_1 \\ \Delta Q_2\end{bmatrix} = -\begin{bmatrix} F_1 \\ F_2 \end{bmatrix}$$

O sistema 2×2 é resolvido pela **regra de Cramer**. Com $D = ad - bc$, as correções são:

$$\Delta Q_1 = \frac{r_1 d - b r_2}{D}, \qquad \Delta Q_2 = \frac{a r_2 - r_1 c}{D}$$

Depois, atualiza-se $Q \leftarrow Q + \Delta Q$, até que $\max(\lvert\Delta Q_1\rvert, \lvert\Delta Q_2\rvert) < 10^{-8}$ m³/s.

A função usada é `newton_sistemas`, de [`utils/sistemas_nao_lineares.py`](../../utils/sistemas_nao_lineares.py). Para $n = 2$ ela usa Cramer; para $n > 2$, Gauss com pivoteamento.

## Dicionário de variáveis

| Variável | Significado | Unidade | Tipo |
|---|---|---|---|
| `Q_TOTAL` | vazão total | m³/s | float |
| `L1`, `D1`, `F_ATRITO1`, `K1` | comprimento, diâmetro, fator de atrito e perda localizada do ramo 1 | m, m, –, – | float |
| `L2`, `D2`, `F_ATRITO2`, `K2`, `KV` | idem para o ramo 2, mais a válvula | m, m, –, –, – | float |
| `R1`, `R2` | coeficientes de resistência ($h_L = RQ^2$) | s²/m⁵ | float |
| `Q_INICIAL` | chute $(Q_1^{(0)}, Q_2^{(0)})$ | m³/s | list |
| `TOL`, `MAX_ITER` | critério e limite de iterações | m³/s, – | float, int |
| `Q` | vazões convergidas $(Q_1, Q_2)$ | m³/s | array |

## Dependências, I/O e execução

- **Dependências:** `numpy`.
- **Saída:** `outputs/resultados.txt`.

```bash
python Aulas_Praticas/AP04_Newton_Sistemas_Tubulacoes/main.py
```

## Resultados e validação

Os coeficientes de resistência são $R_1 = 6197{,}01$ s²/m⁵ e $R_2 = 29250{,}23$ s²/m⁵.

| Grandeza | Valor |
|---|---|
| $Q_1$ | 0,03423990 m³/s (68,5 %) |
| $Q_2$ | 0,01576010 m³/s (31,5 %) |
| $U_1$, $U_2$ | 4,360 m/s e 3,135 m/s |
| $h_{L,1} = h_{L,2}$ | 7,265197 m |
| Iterações | 4 |

A validação tem duas partes:

- **Conferência física:** $Q_1 + Q_2 = 0{,}0500000000$ e $h_{L,1} = h_{L,2}$ até a 8ª casa.
- **Conferência analítica:** de $R_1Q_1^2 = R_2Q_2^2$ sai $Q_1/Q_2 = \sqrt{R_2/R_1}$, ou seja, $Q_1 = Q\sqrt{R_2}/(\sqrt{R_1}+\sqrt{R_2})$. O erro do Newton em relação a esse valor é de $7\times10^{-18}$ m³/s.

## Respostas às questões

**a) Qual ramo leva mais vazão?**
O ramo 1, com 68,5 %. Ele tem resistência $R_1$ cerca de 4,7 vezes menor: diâmetro maior (o $D^4$ no denominador pesa muito) e não tem válvula. Como $Q_1/Q_2 = \sqrt{R_2/R_1} \approx 2{,}17$, a vazão "prefere" o caminho de menor resistência.

**b) Número de iterações?**
4, partindo de (0,030; 0,020).

**c) Conservação e perdas iguais?**
Sim. $F_1$ fica no zero de máquina ($7\times10^{-18}$). $F_2$ cai para $3\times10^{-8}$ m na última iteração avaliada, e as perdas recalculadas coincidem em 7,265197 m.

**d) Por que o Newton transforma o sistema não linear em linear a cada iteração?**
Porque ele troca cada $F_i$ pela série de Taylor de 1ª ordem em torno do ponto atual (o plano tangente). Assim, as incógnitas passam a ser as correções $\Delta Q$, que aparecem linearmente. Resolver um sistema linear é fácil e direto; o não linear, não.

**e) Papel da Jacobiana?**
Ela é a "derivada" do sistema: contém as inclinações de cada equação em relação a cada incógnita, $J_{ij} = \partial F_i / \partial Q_j$. Faz o papel do $f'(x)$ do Newton 1D e define o plano tangente usado na linearização. Se $J$ for singular ($\det J = 0$), o passo não existe.

**f) Outros chutes iniciais?**

| Chute $(Q_1, Q_2)$ | $Q_1$ | $Q_2$ | Iterações |
|---|---|---|---|
| (0,030; 0,020) | 0,034240 | 0,015760 | 4 |
| (0,010; 0,040) | 0,034240 | 0,015760 | 5 |
| (0,100; 0,100) | 0,034240 | 0,015760 | 6 |
| (0,001; 0,001) | 0,034240 | 0,015760 | 5 |
| (0,080; −0,030) | 0,092641 | **−0,042641** | 5 |

Chutes razoáveis convergem para a mesma solução, mudando só o número de iterações. Mas o último chute converge para uma **solução espúria**, com $Q_2 < 0$ (escoamento invertido no ramo 2). Ela satisfaz as equações porque $R_2 Q_2^2$ não enxerga o sinal de $Q_2$. Fisicamente, a perda de carga troca de sinal com o sentido do escoamento. A forma correta seria $h_L = R\,Q\lvert Q\rvert$, que eliminaria essa raiz. **Lição:** sempre conferir se a raiz encontrada tem sentido físico.

## Bibliografia específica

- CHAPRA, S. C.; CANALE, R. P. *Métodos numéricos para engenharia*. 7. ed. AMGH, 2016. Seção 6.6 (sistemas não lineares).
- FOX, R. W.; McDONALD, A. T.; PRITCHARD, P. J. *Introdução à mecânica dos fluidos*. 8. ed. LTC, 2014. Cap. 8 (perdas de carga).
