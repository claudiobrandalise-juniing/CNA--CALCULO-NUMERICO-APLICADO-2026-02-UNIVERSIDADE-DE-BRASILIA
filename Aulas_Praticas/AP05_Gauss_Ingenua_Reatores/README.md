# AP 05 — Eliminação gaussiana ingênua: rede de reatores de mistura

## Resumo operacional

Cinco reatores perfeitamente agitados, ligados por correntes $Q_{ij}$ (da origem $i$ para o destino $j$), operam em regime permanente. O balanço de massa da espécie química em cada reator diz que **entra = sai**, e isso gera cinco equações **lineares** nas concentrações $c_1, \dots, c_5$:

$$\begin{bmatrix} 6 & 0 & -1 & 0 & 0 \\ -3 & 3 & 0 & 0 & 0 \\ 0 & -1 & 9 & 0 & 0 \\ 0 & -1 & -8 & 11 & -2 \\ -3 & -1 & 0 & 0 & 4 \end{bmatrix}\begin{bmatrix} c_1 \\ c_2 \\ c_3 \\ c_4 \\ c_5 \end{bmatrix} = \begin{bmatrix} 50 \\ 0 \\ 160 \\ 0 \\ 0 \end{bmatrix}$$

Cada **linha** é o balanço de um reator e cada **coluna** é uma incógnita: $A_{ij}$ é o coeficiente de $c_j$ na equação $i$. O sistema é resolvido por eliminação gaussiana **sem pivoteamento**, e o próprio programa conta os FLOPs.

> Este código evoluiu do `gausstwo.py` original do repositório, que já estava correto. A lógica agora fica na biblioteca (`utils/sistemas_lineares.py`) e passou a imprimir a matriz depois de **cada** pivô, como o enunciado pede.

## Algoritmo

**Eliminação progressiva** (os índices seguem a apostila):

```
para k = 1 .. n-1            (k: coluna do pivô = linha pivô)
    para i = k+1 .. n        (i: linha onde a_ik será zerado)
        fator = a_ik / a_kk                      (multiplicador m_ik)
        para j = k+1 .. n    (j: colunas atualizadas da linha i)
            a_ij = a_ij - fator * a_kj
        b_i = b_i - fator * b_k
        a_ik = 0
```

**Substituição regressiva**, de baixo para cima:

$$x_n = b_n / a_{nn}, \qquad x_i = \frac{b_i - \sum_{j>i} a_{ij}x_j}{a_{ii}}$$

Se algum $\lvert a_{kk}\rvert$ ficar abaixo de $n\,\varepsilon_{maq}\max\lvert a_{ij}\rvert$, o programa **para** com uma mensagem, porque a versão ingênua não troca linhas.

### Contagem de FLOPs (teoria)

Na etapa $k$ há $m = n - k$ linhas abaixo do pivô. Cada linha custa:

- **Multiplicações e divisões:** 1 divisão (fator), mais $m$ multiplicações nas colunas $j$, mais 1 multiplicação no $b$, ou seja, $m + 2$.
- **Somas e subtrações:** $m$ subtrações nas colunas $j$, mais 1 no $b$, ou seja, $m + 1$.

Somando tudo:

| | Eliminação | Substituição | Total ($n = 5$) |
|---|---|---|---|
| mult./div. | $\sum_{m=1}^{n-1} m(m+2) = 50$ | $n + n(n-1)/2 = 15$ | **65** |
| soma/sub. | $\sum_{m=1}^{n-1} m(m+1) = 40$ | $n(n-1)/2 = 10$ | **50** |

Total: **115 FLOPs**. Para $n$ grande, o custo é $\approx 2n^3/3$.

## Dicionário de variáveis

| Variável | Significado | Unidade | Tipo |
|---|---|---|---|
| `VAZAO` | dicionário com as vazões $Q_{ij}$ | m³/min | dict |
| `C01`, `C03` | concentrações das alimentações externas | mg/m³ | float |
| `A`, `b` | matriz dos coeficientes e vetor independente | m³/min, mg/min | array 5×5, array 5 |
| `c` | concentrações calculadas | mg/m³ | array |
| `flops` | contadores `mult_div` e `soma_sub` | – | dict |
| `r` | resíduo $A_{orig}c - b_{orig}$ | mg/min | array |

## Dependências, I/O e execução

- **Dependências:** `numpy`.
- **Saída:** `outputs/resultados.txt`, com o sistema, a matriz após cada pivô, as concentrações, o resíduo, os FLOPs e o teste de pivô nulo.

```bash
python Aulas_Praticas/AP05_Gauss_Ingenua_Reatores/main.py
```

## Resultados e validação

| Reator | $c$ [mg/m³] | Esperado |
|---|---|---|
| 1 | 11,509434 | 11,5094 |
| 2 | 11,509434 | 11,5094 |
| 3 | 19,056604 | 19,0566 |
| 4 | 16,998285 | 16,9983 |
| 5 | 11,509434 | 11,5094 |

O resíduo é $\lvert\lvert r\rvert\rvert_\infty = 7{,}1\times10^{-15}$, no nível do arredondamento de máquina. Os FLOPs contados pelo programa (65 + 50 = 115) batem **exatamente** com a análise dos laços.

## Respostas às questões

**a) Significado físico de $c_1, \dots, c_5$?**
A concentração (uniforme, porque a mistura é perfeita) da espécie química dentro de cada reator. Como cada reator é perfeitamente agitado, essa também é a concentração da corrente que sai dele.

**b) Por que o problema gera um sistema linear?**
O fluxo de massa de cada corrente é vazão vezes concentração, $Q\,c$. As vazões são conhecidas e as concentrações entram só elevadas a 1, sem produtos entre si. Em regime permanente, a soma das entradas é igual à soma das saídas em cada reator, o que dá $n$ equações lineares acopladas.

**c) Função do índice $k$?**
Indica o pivô: a coluna que está sendo zerada e a linha usada para zerá-la. Ele fica fixo enquanto toda a coluna abaixo de $a_{kk}$ é eliminada.

**d) Função do índice $i$?**
Percorre as linhas abaixo do pivô ($i = k+1..n$), em que o elemento $a_{ik}$ será eliminado.

**e) Função do índice $j$?**
Percorre as colunas da linha $i$ que precisam ser atualizadas ($j = k+1..n$) pela operação $L_i \leftarrow L_i - m_{ik}L_k$.

**f) Por que laços aninhados?**
Cada pivô ($k$) atua sobre várias linhas ($i$), e cada linha tem vários elementos ($j$): são três níveis de repetição. É daí que vem o custo $\propto n^3$.

**g) O que é o multiplicador $m_{ik}$?**
É o fator $a_{ik}/a_{kk}$ que, multiplicado pela linha pivô e subtraído da linha $i$, zera o elemento $a_{ik}$. Na decomposição LU, ele é exatamente o elemento $\ell_{ik}$ da matriz $L$.

**h) Por que a substituição começa pela última equação?**
Depois da triangularização, a última equação tem uma única incógnita ($c_n$). Cada equação acima só depende de incógnitas que já foram calculadas abaixo dela.

**i) O resíduo foi pequeno?**
Sim, $\lvert\lvert r\rvert\rvert_\infty \approx 7\times10^{-15}$, compatível com a precisão dupla ($\varepsilon \approx 2{,}2\times10^{-16}$).

**j) E se um pivô fosse zero?**
O multiplicador $a_{ik}/a_{kk}$ seria uma divisão por zero e o algoritmo não poderia continuar. O programa demonstra isso: trocando a ordem das equações 1 e 3, o primeiro pivô vira $a_{11} = 0$, e a execução é interrompida com aviso. A solução é o **pivoteamento** (AP 06). Pivôs só *pequenos* também são perigosos, porque amplificam os erros de arredondamento.

## Bibliografia específica

- CHAPRA, S. C.; CANALE, R. P. *Métodos numéricos para engenharia*. 7. ed. AMGH, 2016. Cap. 9 (eliminação gaussiana) e seção 12.1 (reatores).
