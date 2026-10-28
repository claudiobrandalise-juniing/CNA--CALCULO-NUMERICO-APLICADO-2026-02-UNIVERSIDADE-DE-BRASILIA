# APC 03 — Sistemas lineares feitos à mão (LU, Thomas, Cholesky, Gauss-Seidel)

> A APC é entregue em PDF (LaTeX). Esta pasta tem o script [`verificar_apc3.py`](verificar_apc3.py), que refaz **cada passo** das contas pedidas para você conferir o caderno. A saída completa está em [`outputs/verificacao.txt`](outputs/verificacao.txt).
>
> **Atenção:** o PDF do enunciado termina com um item "1. Mostre que a substituição da solução proposta..." que foi copiado por engano da APC 2 e não pertence a esta atividade.

```bash
python APC_03/verificar_apc3.py
```

## Item 1 — Inversa por LU

$$A = \begin{bmatrix}3 & -0{,}1 & -0{,}2\\ 0{,}1 & 7 & -0{,}3\\ 0{,}3 & -0{,}2 & 10\end{bmatrix} \;\Rightarrow\; L = \begin{bmatrix}1&0&0\\0{,}0333333&1&0\\0{,}1&-0{,}0271300&1\end{bmatrix},\quad U = \begin{bmatrix}3&-0{,}1&-0{,}2\\0&7{,}0033333&-0{,}2933333\\0&0&10{,}012042\end{bmatrix}$$

Para cada coluna, resolve-se $Ly = e_j$ (progressiva) e depois $Ux = y$ (regressiva):

$$A^{-1} = \begin{bmatrix}0{,}332489 & 0{,}004944 & 0{,}006798\\ -0{,}005182 & 0{,}142903 & 0{,}004183\\ -0{,}010078 & 0{,}002710 & 0{,}099880\end{bmatrix}, \qquad \max\lvert AA^{-1} - I\rvert = 1{,}1\times10^{-16}$$

## Item 2 — Tridiagonal: LU completa x Thomas

O sistema tem $f_k = 2{,}04$, $e_k = g_k = -1$ e $b = [40{,}8;\ 0{,}8;\ 0{,}8;\ 200{,}8]$.

**Tabela dos laços de Thomas:**

| Laço | k | $e_k$ | $f_k$ | $r_k$ | $x_k$ |
|---|---|---|---|---|---|
| decomposição | 2 | −0,490196 | 1,549804 | | |
| decomposição | 3 | −0,645243 | 1,394757 | | |
| decomposição | 4 | −0,716971 | 1,323029 | | |
| progressiva | 2 | | | 20,8 | |
| progressiva | 3 | | | 14,221053 | |
| progressiva | 4 | | | 210,996078 | |
| regressiva | 4 | | | | 159,4795 |
| regressiva | 3 | | | | 124,5382 |
| regressiva | 2 | | | | 93,7785 |
| regressiva | 1 | | | | 65,9698 |

Os dois caminhos dão a **mesma** solução, $x = [65{,}970;\ 93{,}778;\ 124{,}538;\ 159{,}480]$. Os $e_k$ do Thomas são exatamente os multiplicadores $\ell_{k,k-1}$ da LU, e os $f_k$ são a diagonal de $U$: conceitualmente é o mesmo método.

| Método | mult./div. | soma/sub. | Total |
|---|---|---|---|
| LU completa (decomposição + 2 substituições) | 40 | 26 | 66 |
| Thomas | 16 | 9 | **25** |

No caso geral, o Thomas faz $8n - 7$ operações (linear em $n$), contra $\sim 2n^3/3$ da LU cheia, que opera inclusive sobre os zeros. Ele também guarda só $3n - 2$ números, em vez de $n^2$.

## Item 3 — Cholesky (exemplo, troque pelos seus números)

Escolhendo $b$ para que a solução seja $x = (1, 2, 3)$:

$$A = \begin{bmatrix}6&15&55\\15&55&225\\55&225&979\end{bmatrix},\quad b = \begin{bmatrix}201\\800\\3442\end{bmatrix}$$

Aplicando as recorrências $\ell_{ki} = (a_{ki} - \sum_{j<i}\ell_{ij}\ell_{kj})/\ell_{ii}$ e $\ell_{kk} = \sqrt{a_{kk} - \sum_{j<k}\ell_{kj}^2}$:

$$L = \begin{bmatrix}2{,}449490&0&0\\6{,}123724&4{,}183300&0\\22{,}453656&20{,}916501&6{,}110101\end{bmatrix}$$

Resolvendo $Ly = b$ chega-se a $y = [82{,}0579;\ 71{,}1161;\ 18{,}3303]$. Resolvendo $L^Tx = y$ chega-se a $x = [1;\ 2;\ 3]$.

> Na apostila, a eq. (40) traz $a_{22} = \ell_{21}^2\ell_{22}^2$. O correto é $a_{22} = \ell_{21}^2 + \ell_{22}^2$, como na eq. (39).

## Item 4 — Gauss-Seidel com erro relativo de 5 %

Na ordem original, a matriz **não** é diagonalmente dominante. Reordenando para (eq. 3, eq. 1, eq. 2):

$$\begin{cases}-8x_1 + x_2 - 2x_3 = -20\\ 2x_1 - 6x_2 - x_3 = -38\\ -3x_1 - x_2 + 7x_3 = -34\end{cases}$$

Agora $\lvert a_{ii}\rvert > \sum_{j\ne i}\lvert a_{ij}\rvert$ em todas as linhas, o que garante a convergência.

| $\lambda$ | Iterações até $\varepsilon_a < 5\%$ | Resultado |
|---|---|---|
| 1,0 (sem relaxação) | 3 | (4,0047; 7,9917; −1,9992) |
| 1,2 (sobre-relaxação) | 6 | (4,0122; 8,0273; −1,9790) |

A solução exata é $(4, 8, -2)$. Aqui a sobre-relaxação **piorou**. O sistema já converge muito rápido com Gauss-Seidel puro, e $\lambda = 1{,}2$ faz a solução "passar do ponto" e oscilar em torno da resposta (veja $x_1$: 3,0, depois 4,55, depois 3,78, depois 4,08). O SOR ajuda em sistemas de convergência **lenta**, como as malhas grandes de Laplace.
