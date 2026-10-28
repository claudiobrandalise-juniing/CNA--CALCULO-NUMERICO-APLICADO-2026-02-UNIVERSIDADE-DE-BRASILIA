# APC 02 — Polinômio característico de um sistema massa-mola-amortecedor

> As APCs são entregues em PDF (LaTeX). Este arquivo guarda a dedução para consulta. O polinômio obtido aqui é resolvido no [PPC 02](../PPC_02_Bairstow/README.md), que também confere automaticamente os coeficientes.

## 1. Do sistema de EDOs ao polinômio de ordem 2N

A equação do movimento é:

$$M\ddot{x} + C\dot{x} + Kx = 0$$

Propondo $x(t) = \phi\,e^{\lambda t}$, temos $\dot x = \lambda\phi e^{\lambda t}$ e $\ddot x = \lambda^2\phi e^{\lambda t}$. Como $e^{\lambda t} \ne 0$:

$$(\lambda^2 M + \lambda C + K)\,\phi = 0$$

Esse é o **problema de autovalores quadrático**. Para existir uma solução não trivial ($\phi \ne 0$), a matriz precisa ser singular:

$$P(\lambda) = \det(\lambda^2 M + \lambda C + K) = 0$$

Cada elemento da matriz é um polinômio de grau 2 em $\lambda$. O determinante de uma matriz $N \times N$ soma produtos de $N$ elementos, então tem grau $2N$. O coeficiente de $\lambda^{2N}$ é $\det M = \prod m_i \ne 0$, logo o grau é **exatamente** $2N$: são $2N$ autovalores.

## 2. Significado físico das raízes

Como os coeficientes são reais, as raízes complexas aparecem em **pares conjugados** $\lambda = \sigma \pm i\omega_d$. Cada par corresponde a um modo:

$$x(t) \propto e^{\sigma t}\cos(\omega_d t + \varphi)$$

| Raiz | Comportamento |
|---|---|
| $\sigma < 0$ | o modo decai: estável, amortecido |
| $\sigma > 0$ | o modo cresce: instável |
| $\sigma = 0$ ($\lambda$ imaginário puro) | oscilação sem amortecimento |
| $\lambda$ real negativo | modo superamortecido (decai sem oscilar) |

Para cada par: frequência amortecida $\omega_d = \lvert\text{Im}\,\lambda\rvert$; frequência natural $\omega_n = \lvert\lambda\rvert$; fator de amortecimento $\zeta = -\sigma/\omega_n$; constante de tempo $\tau = 1/\lvert\sigma\rvert$.

O sistema é **estável** se, e somente se, todas as raízes têm parte real negativa.

## 3. Sistema proposto (N = 3)

A cadeia é: parede – ($k_1$, $c_1$) – $m_1$ – ($k_2$, $c_2$) – $m_2$ – ($k_3$, $c_3$) – $m_3$, com extremo livre.

| | 1 | 2 | 3 |
|---|---|---|---|
| $m_i$ [kg] | 2,0 | 1,5 | 1,0 |
| $k_i$ [N/m] | 800 | 600 | 400 |
| $c_i$ [N s/m] | 4 | 3 | 2 |

A escolha $c_i = k_i/200$ dá **amortecimento proporcional**, $C = \beta K$ com $\beta = 0{,}005$ s. Isso permite uma validação analítica no PPC: os modos amortecidos têm as mesmas $\omega_n$ do sistema sem amortecimento, e $\zeta_i = \beta\,\omega_{n,i}/2$.

## 4. Montagem explícita e expansão do determinante

As matrizes são:

$$M = \begin{bmatrix}2&0&0\\0&1{,}5&0\\0&0&1\end{bmatrix},\quad C = \begin{bmatrix}7&-3&0\\-3&5&-2\\0&-2&2\end{bmatrix},\quad K = \begin{bmatrix}1400&-600&0\\-600&1000&-400\\0&-400&400\end{bmatrix}$$

Portanto:

$$\lambda^2M + \lambda C + K = \begin{bmatrix} 2\lambda^2+7\lambda+1400 & -3\lambda-600 & 0\\ -3\lambda-600 & 1{,}5\lambda^2+5\lambda+1000 & -2\lambda-400 \\ 0 & -2\lambda-400 & \lambda^2+2\lambda+400\end{bmatrix}$$

Expandindo pela 1ª linha, e usando $A_{13} = 0$:

$$P = A_{11}\,(A_{22}A_{33} - A_{23}^2) - A_{12}^2\,A_{33}$$

Os termos intermediários são:

- $A_{22}A_{33} = 1{,}5\lambda^4 + 8\lambda^3 + 1610\lambda^2 + 4000\lambda + 400\,000$
- $A_{23}^2 = 4\lambda^2 + 1600\lambda + 160\,000$
- $A_{22}A_{33} - A_{23}^2 = 1{,}5\lambda^4 + 8\lambda^3 + 1606\lambda^2 + 2400\lambda + 240\,000$
- $A_{11}(\cdots) = 3\lambda^6 + 26{,}5\lambda^5 + 5368\lambda^4 + 27\,242\lambda^3 + 2\,745\,200\lambda^2 + 5\,040\,000\lambda + 336\,000\,000$
- $A_{12}^2A_{33} = 9\lambda^4 + 3618\lambda^3 + 370\,800\lambda^2 + 2\,160\,000\lambda + 144\,000\,000$

## 5. Forma padrão (pronta para o PPC 2)

$$\boxed{P(\lambda) = 3\lambda^6 + 26{,}5\lambda^5 + 5359\lambda^4 + 23\,624\lambda^3 + 2\,374\,400\lambda^2 + 2\,880\,000\lambda + 192\,000\,000}$$

| $a_6$ | $a_5$ | $a_4$ | $a_3$ | $a_2$ | $a_1$ | $a_0$ |
|---|---|---|---|---|---|---|
| 3 | 26,5 | 5359 | 23 624 | 2 374 400 | 2 880 000 | 192 000 000 |

Há algumas conferências rápidas, que valem para qualquer cadeia:

- $a_6 = m_1m_2m_3 = 3$.
- $a_0 = \det K = k_1k_2k_3 = 1{,}92\times10^8$.
- $a_1 = c_1k_2k_3 + k_1c_2k_3 + k_1k_2c_3 = 2{,}88\times10^6$, porque $\det(\lambda C + K) = \prod(k_i + \lambda c_i)$.
- $a_5 = m_1m_2c_3 + m_1m_3(c_2+c_3) + m_2m_3(c_1+c_2) = 26{,}5$.

No código, os coeficientes ficam em ordem **crescente**: `[192000000, 2880000, 2374400, 23624, 5359, 26.5, 3]` (arquivo `PPC_02_Bairstow/inputs/parametros.json`).
