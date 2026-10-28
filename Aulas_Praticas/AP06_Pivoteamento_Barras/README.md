# AP 06 — Condicionamento e pivoteamento parcial em uma estrutura de barras

## Resumo operacional

Barras articuladas só transmitem esforço axial ($F > 0$ é tração). O equilíbrio de forças em cada nó gera um sistema linear $AF = b$. A atividade separa dois fenômenos que costumam ser confundidos:

- **Condicionamento**, que é uma propriedade do *problema*. Quando duas barras ficam quase paralelas, pequenas variações nas cargas geram grandes variações nos esforços. Isso acontece com qualquer algoritmo.
- **Estabilidade**, que é uma propriedade do *algoritmo*. A eliminação sem pivoteamento pode dividir por pivôs muito pequenos e amplificar o erro de arredondamento.

**Sistema 2×2 (nó A):**

$$\begin{bmatrix}\sin\theta_1 & \sin\theta_2\\ \cos\theta_1 & \cos\theta_2\end{bmatrix}\begin{bmatrix}F_1\\F_2\end{bmatrix} = -\begin{bmatrix}P_x\\P_y\end{bmatrix}, \qquad \det A_2 = \sin\theta_1\cos\theta_2 - \sin\theta_2\cos\theta_1 = \sin(\theta_1 - \theta_2)$$

**Sistema 5×5 (nós A, B, C):** é escrito na ordem do enunciado, de propósito com o pivô inicial $\sin\theta_1$, minúsculo quando $\theta_1 \to 0$.

## Algoritmo: eliminação com pivoteamento parcial

Em cada etapa $k$:

1. Procura $p = \arg\max_{i \ge k} \lvert a_{ik}\rvert$.
2. Se $p \ne k$, troca as linhas $k$ e $p$ de $A$ e de $b$.
3. Se o maior $\lvert a_{pk}\rvert$ ainda for "zero", para (matriz singular).
4. Elimina a coluna $k$ normalmente.

O critério de pivô nulo é **relativo à escala**: $\lvert a_{kk}\rvert \le n\,\varepsilon_{maq}\max\lvert a_{ij}\rvert$. Ele funciona igual se as cargas estiverem em N ou em mN, ao contrário de um limiar fixo.

As funções vêm de [`utils/sistemas_lineares.py`](../../utils/sistemas_lineares.py): `gauss_ingenua` e `gauss_pivoteamento`. O parâmetro `trocas` registra as trocas de linha, e `algarismos=4` arredonda **cada operação** para simular uma aritmética de 4 dígitos (imprimir menos casas não simula precisão reduzida).

## Dicionário de variáveis

| Variável | Significado | Unidade | Tipo |
|---|---|---|---|
| `CONFIGURACOES` | pares $(\theta_1, \theta_2)$ das configurações A, B e C | graus | dict |
| `PHI` | ângulo $\varphi$ da barra 4 | graus | float |
| `F_REF_2`, `F_REF_5` | esforços de referência (usados só para gerar $b$ e conferir) | N | array |
| `PERTURBACAO` | acréscimo aplicado a uma carga | N | float |
| `A`, `b` | matriz de equilíbrio e vetor de cargas | –, N | array |
| `trocas` | lista de trocas de linha $(k, p)$ | – | list |
| `cond_inf` | número de condição $\lVert A\rVert_\infty\lVert A^{-1}\rVert_\infty$ | – | float |

## Dependências, I/O e execução

- **Dependências:** `numpy`.
- **Saída:** `outputs/resultados.txt`, com as etapas da eliminação, os erros, os resíduos, a sensibilidade e os FLOPs.

```bash
python Aulas_Praticas/AP06_Pivoteamento_Barras/main.py
```

## Resultados

### Parte 1: sistema 2×2 (perturbação de +0,01 N)

| Config. | $\lvert\det\rvert$ | $\Delta F$ ($p_y$) | $\Delta F$ ($p_x$) | cond∞ |
|---|---|---|---|---|
| A (15°, 65°) | 0,766 | ~0,01 N | ~0,01 N | 2,4 |
| B (10°, 10,1°) | 1,75×10⁻³ | ~1 N | ~5,6 N | 1309 |
| C (0,01°, 0,02°) | 1,75×10⁻⁴ | 0,02 N | **57 N** | 11463 |

O resultado mais interessante está na configuração C: perturbar $p_y$ quase não muda nada, mas perturbar $p_x$ amplifica o erro cerca de 5700 vezes. As barras estão quase **verticais**, então uma carga vertical é resistida "ao longo" delas com facilidade. Já uma carga horizontal precisa ser equilibrada pela pequena diferença entre os ângulos. O número de condição mede o **pior caso**; uma perturbação específica pode ou não excitar a direção mal condicionada.

O determinante também depende da escala. Escrevendo a mesma matriz em mN (×1000), o $\det$ muda $10^6$ vezes e o cond∞ continua igual.

### Partes 2 a 4: sistema 5×5

**Precisão dupla:**

| Config. | Método | $\lVert F - F_{ref}\rVert_\infty$ | $\lVert r\rVert_\infty$ |
|---|---|---|---|
| A | ingênua | 4×10⁻¹⁴ | 3×10⁻¹⁴ |
| A | pivoteamento (troca 1↔2) | 7×10⁻¹⁵ | 1×10⁻¹⁴ |
| C | ingênua | 6×10⁻¹¹ | 5×10⁻¹² |
| C | pivoteamento (troca 1↔2) | 5×10⁻¹¹ | 4×10⁻¹⁵ |

**Aritmética de $d$ algarismos (configuração C):**

| $d$ | Ingênua: erro / resíduo | Pivoteamento: erro / resíduo |
|---|---|---|
| 4 | 20 N / 5,4 | 22,7 N / 0,004 |
| 6 | 0,27 N / 0,27 | 0,21 N / 1×10⁻⁴ |
| 8 | 1,7×10⁻³ N / 1,7×10⁻³ | 2,1×10⁻³ N / 1×10⁻⁶ |

O pivoteamento sempre produz um **resíduo** muito menor (o algoritmo fica estável), mas o **erro** em $F$ continua grande. Com cond∞ ≈ 4,6×10⁴, perdem-se cerca de 4,7 algarismos, e com 4 dígitos não sobra nenhum.

### Itens 8 e 9

Perturbar $b_2$ amplifica 1,2× (A) e 2× (C). Perturbar $b_1$ amplifica 1,3× (A) e **5730×** (C). As duas versões fazem os mesmos **115 FLOPs** (65 mult./div. e 50 soma/sub.). O pivoteamento só acrescenta comparações e trocas, que não são contadas.

## Respostas às questões

**a) Relação entre a diferença angular e o determinante?**
$\det A_2 = \sin(\theta_1 - \theta_2)$, que vai a zero quando as barras ficam paralelas. Mas **o determinante sozinho não mede condicionamento**, porque depende das unidades e da escala das equações (multiplicar a matriz por 1000 multiplica o det 2×2 por $10^6$). A medida correta é $\text{cond}(A) = \lVert A\rVert\,\lVert A^{-1}\rVert$, que não muda com a escala.

**b) Quais linhas foram pivô em cada cenário?**
Nos dois cenários, na etapa $k = 1$ a linha 2 (equilíbrio vertical em A) virou pivô, porque $\lvert\cos\theta_1\rvert > \lvert\sin\theta_1\rvert$ (0,966 contra 0,259 em A, e ≈1 contra 1,7×10⁻⁴ em C). Nas outras etapas não houve troca: abaixo da diagonal só havia zeros. Em C, o 2º pivô continua sendo $1{,}7\times10^{-4}$, porque o mau condicionamento está no próprio problema.

**c) O pivoteamento alterou o sistema físico?**
Não. Ele só muda a **ordem** em que as equações são processadas. As equações, a solução e o significado físico continuam os mesmos.

**d) Um resíduo pequeno garante forças próximas das exatas?**
Não. Com 4 algarismos e pivoteamento, $\lVert r\rVert_\infty = 0{,}004$, mas o erro em $F$ é de 22,7 N. A relação é $\frac{\lVert\Delta F\rVert}{\lVert F\rVert} \le \text{cond}(A)\,\frac{\lVert r\rVert}{\lVert b\rVert}$. Num problema mal condicionado, um resíduo pequeno é compatível com um erro grande.

**e) Estabilidade do algoritmo x condicionamento do problema?**
*Estabilidade* é uma propriedade do método: o quanto ele amplifica os erros de arredondamento. O pivoteamento torna a eliminação estável, e isso aparece como resíduo pequeno. *Condicionamento* é uma propriedade dos dados ($A$): o quanto a solução exata muda quando $A$ ou $b$ mudam, seja qual for o algoritmo. Um bom algoritmo não conserta um problema mal condicionado.

**f) E se $\theta_1 = \theta_2$?**
Então $\det = 0$ e a matriz fica **singular**: as duas colunas são iguais. As barras colineares não resistem a uma carga perpendicular a elas, e o nó é um mecanismo. Para cargas exatamente na direção das barras existem infinitas soluções (estrutura hiperestática nessa direção). O programa detecta isso como "pivô nulo / matriz singular".

## Bibliografia específica

- CHAPRA, S. C.; CANALE, R. P. *Métodos numéricos para engenharia*. 7. ed. AMGH, 2016. Seções 9.3–9.4 (pivoteamento) e 10.3 (condicionamento).
- HIBBELER, R. C. *Estática: mecânica para engenharia*. 12. ed. Pearson, 2011. Cap. 6 (treliças).
