# Guia da prova — consulta rápida

> Prova prática: escrever um programa durante a prova, consultando apenas os próprios códigos. Tema anunciado: **temperatura em processadores**, com todo o conteúdo (zeros de funções, sistemas lineares e não lineares, otimização, EDOs e EDPs).

**Sumário:** [1. Qual método usar](#1-qual-método-usar) · [2. Como importar](#2-como-importar-a-biblioteca) · [3. Funções da biblioteca](#3-funções-da-biblioteca-assinaturas) · [4. Fórmulas](#4-fórmulas-de-cada-método) · [5. Convergência e estabilidade](#5-convergência-e-estabilidade) · [6. Modelagem térmica](#6-modelagem-térmica-de-processadores) · [7. Armadilhas](#7-armadilhas-que-custam-pontos) · [8. Roteiro](#8-roteiro-de-prova-em-8-passos)

---

## 1. Qual método usar

| O problema pede… | Forma matemática | Função | Exemplo pronto |
|---|---|---|---|
| uma temperatura de equilíbrio (radiação $T^4$, fuga exponencial…) | $f(x) = 0$ | `bissecao`, `falsa_posicao`, `newton_raphson`, `secante`, `muller` | [E01](Exemplos_Processador/E01_temperatura_equilibrio_dissipador.py), [AP01](Aulas_Praticas/AP01_Bissecao_FalsaPosicao), [AP02](Aulas_Praticas/AP02_Newton_Secante) |
| raízes de um polinômio / constantes de tempo / frequências | $P(\lambda) = 0$ | `bairstow`, `muller`, `determinante_polinomial` | [E08](Exemplos_Processador/E08_constantes_tempo_bairstow.py), [PPC 2](PPC_02_Bairstow), [AP03](Aulas_Praticas/AP03_Muller_Canal_Critico) |
| temperaturas de vários nós (rede de resistências) | $A\,x = b$ | `gauss_pivoteamento`, `decomposicao_lu`/`resolver_lu`, `cholesky`, `gauss_seidel` | [E02](Exemplos_Processador/E02_rede_resistencias_multinucleo.py), [AP05](Aulas_Praticas/AP05_Gauss_Ingenua_Reatores), [AP06](Aulas_Praticas/AP06_Pivoteamento_Barras) |
| "quanto sobe cada nó se…" / muitos casos de carga | $A^{-1}$ | `inversa_lu` | [E02](Exemplos_Processador/E02_rede_resistencias_multinucleo.py), [AP07](Aulas_Praticas/AP07_Inversao_LU_Molas) |
| sistema **tridiagonal** (1D, diferenças finitas) | $e, f, g$ | `thomas` | [PPC 3](PPC_03_Conducao_Transiente), [APC 3](APC_03) |
| várias equações **não lineares** acopladas | $F(x) = 0$ | `newton_sistemas` | [E06](Exemplos_Processador/E06_fuga_termica_newton_sistemas.py), [AP04](Aulas_Praticas/AP04_Newton_Sistemas_Tubulacoes) |
| mapa de temperatura 2D em regime permanente | Laplace/Poisson | `gauss_seidel_laplace_2d` | [E03](Exemplos_Processador/E03_mapa_termico_spreader_2d.py), [AP08](Aulas_Praticas/AP08_Gauss_Seidel_Placa) |
| evolução no tempo de blocos com T uniforme | $\dot y = f(t, y)$ | `rk4` (`euler`, `heun`) | [E04](Exemplos_Processador/E04_transiente_rc_rk4.py), [PPC 1](PPC_01_RungeKutta_Sedimentacao) |
| evolução no tempo **e** no espaço (1D) | $T_t = \alpha T_{xx} + \dot q/\rho c$ | `conducao_transiente_1d` | [E05](Exemplos_Processador/E05_transiente_die_implicito.py), [PPC 3](PPC_03_Conducao_Transiente) |
| "melhor", "ótimo", "mínimo custo/energia" | $\min f(x)$ | `razao_aurea`, `interpolacao_quadratica`, `newton_otimizacao` | [E07](Exemplos_Processador/E07_velocidade_otima_ventoinha.py), [AP09](Aulas_Praticas/AP09_Otimizacao_Isolamento) |

---

## 2. Como importar a biblioteca

**Arquivo novo na raiz do repositório:**

```python
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))          # raiz do repo
import numpy as np
import matplotlib.pyplot as plt
from utils.raizes import bissecao, falsa_posicao, newton_raphson, secante, muller
from utils.polinomios import bairstow, avaliar_polinomio, determinante_polinomial
from utils.sistemas_lineares import gauss_pivoteamento, decomposicao_lu, resolver_lu, inversa_lu, thomas, cholesky, resolver_cholesky
from utils.sistemas_iterativos import gauss_seidel, jacobi, dominancia_diagonal, reordenar_para_dominancia
from utils.sistemas_nao_lineares import newton_sistemas
from utils.otimizacao import razao_aurea, interpolacao_quadratica, newton_otimizacao, FuncaoContada
from utils.edo import rk4
from utils.edp import gauss_seidel_laplace_2d, conducao_transiente_1d
from utils.saida import Relatorio
```

**Arquivo numa subpasta** (por exemplo, `Exemplos_Processador/`): use `.parent.parent` no lugar de `.parent`.

**Se não der para importar** e só der para consultar, copie o **módulo inteiro** de `utils/`, que leva junto as funções auxiliares (as que começam com `_`). Os módulos dependem só do `numpy`, com estas exceções:

- `polinomios`, `sistemas_nao_lineares` e `edp` usam `sistemas_lineares`;
- `otimizacao` usa `derivadas`.

---

## 3. Funções da biblioteca (assinaturas)

Todos os métodos iterativos devolvem `(resultado, historico)`. O `historico` é uma lista de dicionários e pode ser impresso com `Relatorio().tabela(historico)`.

```python
# ---------- raízes (utils/raizes.py) ----------  criterio: "funcao" |f|, "passo" |dx|, "relativo" %
x, h = bissecao(f, xl, xu, tol=1e-6, max_iter=200, criterio="funcao")
x, h = falsa_posicao(f, xl, xu, ...)              # e falsa_posicao_modificada(...)
x, h = newton_raphson(f, df, x0, tol=1e-6, max_iter=100, criterio="passo")
x, h = secante(f, x_anterior, x_atual, ...)
x, h = muller(f, x0, x1, x2, tol=1e-8)            # acha raízes complexas também
x, h = ponto_fixo(g, x0)                          # x = g(x)
intervalos = varredura_mudanca_sinal(f, a, b, n_subintervalos=1000)   # "método gráfico"

# ---------- polinômios (utils/polinomios.py) ----------  coef em ordem CRESCENTE [a0, a1, ..., an]
raizes, h = bairstow(coef, r0, s0, tol=1e-10)
y = avaliar_polinomio(coef, x)                    # Horner
coef = polinomio_de_raizes([r1, r2, ...])
coef = determinante_polinomial(M)                 # M[i][j] = [c0, c1, c2] (polinômio em lambda)

# ---------- sistemas lineares (utils/sistemas_lineares.py) ----------
x = gauss_ingenua(A, b, flops=None, mostrar_etapas=False)
x = gauss_pivoteamento(A, b, trocas=None, algarismos=None)   # algarismos=4 -> precisão reduzida
x = gauss_jordan(A, b);  x = cramer(A, b);  d = determinante(A)
L, U = decomposicao_lu(A);  x = resolver_lu(L, U, b);  Ainv = inversa_lu(A)
L, U = decomposicao_crout(A)
x = thomas(e, f, g, r)                            # e[0] e g[n-1] não são usados
L = cholesky(A);  x = resolver_cholesky(L, b)
r = residuo(A, x, b);  flops = novo_contador_flops()   # passe flops=flops para contar

# ---------- iterativos (utils/sistemas_iterativos.py) ----------
x, h = gauss_seidel(A, b, x0=None, tol=1e-6, relaxacao=1.0, criterio="relativo")  # ou "absoluto"
x, h = jacobi(A, b, ...)
A2, b2, ordem = reordenar_para_dominancia(A, b)

# ---------- não lineares (utils/sistemas_nao_lineares.py) ----------
x, h = newton_sistemas(F, J, x0, tol=1e-8)        # J=None -> Jacobiana numérica

# ---------- otimização (utils/otimizacao.py) ----------  MINIMIZA (para máximo use -f)
x, h = razao_aurea(f, a, b, tol=1e-6)
x, h = interpolacao_quadratica(f, x1, x2, x3, tol=1e-6)   # exige f(x2) < f(x1), f(x3)
x, h = newton_otimizacao(f, x0)                    # derivadas numéricas se não passar df, d2f
fc = FuncaoContada(f);  ...;  fc.n_avaliacoes

# ---------- EDO (utils/edo.py) ----------
t, Y = rk4(f, t0, y0, h, n_passos)                 # f(t, y) -> array; Y[:, i] = variável i

# ---------- EDP (utils/edp.py) ----------  T[i, j]: i = x, j = y, j = 0 é a borda de BAIXO
T, erros = gauss_seidel_laplace_2d(T_com_bordas, tol=1e-6, relaxacao=1.0, fonte=S)  # S = q''' dx²/k
x, t, T = conducao_transiente_1d(L, n_nos, alfa, k, dt, n_passos, T_inicial,
                                 contorno_esq=("adiabatico",),          # ou ("simetria",)
                                 contorno_dir=("conveccao", h, T_inf),  # ou ("temperatura", Tf), ("fluxo", q)
                                 geracao=0.0, beta=0.0)                  # 0 implícito, 0.5 CN, 1 explícito
```

---

## 4. Fórmulas de cada método

### Zeros de funções

| Método | Atualização | Ordem | Precisa |
|---|---|---|---|
| Bisseção | $x_r = (x_l + x_u)/2$ | linear (erro ÷2 por passo) | intervalo com troca de sinal |
| Falsa posição | $x_r = x_u - f(x_u)\dfrac{x_l - x_u}{f(x_l) - f(x_u)}$ | linear (pode travar um lado) | intervalo |
| Newton-Raphson | $x_{k+1} = x_k - f/f'$ | 2 | $f'$ e chute bom |
| Secante | $x_{k+1} = x_k - f_k\dfrac{x_k - x_{k-1}}{f_k - f_{k-1}}$ | 1,62 | 2 chutes |
| Müller | parábola por 3 pontos, $x_3 = x_2 - 2c/(b \pm \sqrt{b^2 - 4ac})$ | 1,84 | 3 chutes |

- Iterações da bisseção para um erro $E$: $n = \log_2[(x_u - x_l)/E]$.
- Erro relativo aproximado: $\varepsilon_a = \lvert(x_{novo} - x_{velho})/x_{novo}\rvert \times 100\,\%$.

### Bairstow (fator $x^2 - rx - s$)

- Divisão sintética: $b_n = a_n$, $b_{n-1} = a_{n-1} + rb_n$, $b_i = a_i + rb_{i+1} + sb_{i+2}$.
- Derivadas: $c_i$ com a mesma recorrência sobre os $b_i$.
- Correção (Cramer): $c_2\Delta r + c_3\Delta s = -b_1$ e $c_1\Delta r + c_2\Delta s = -b_0$.
- Raízes do fator: $x = \left(r \pm \sqrt{r^2 + 4s}\right)/2$.

### Sistemas lineares

| Método | Custo | Quando usar |
|---|---|---|
| Gauss (com pivoteamento) | $\sim 2n^3/3$ | sistema cheio, uma vez |
| LU + substituições | $2n^3/3$ uma vez + $2n^2$ por lado direito | muitos $b$ diferentes; inversa |
| Cholesky ($A = LL^T$) | $\sim n^3/3$ | $A$ simétrica positiva definida (matriz de condutâncias!) |
| Thomas | $\sim 8n$ | tridiagonal (FD 1D) |
| Gauss-Seidel / SOR | por iteração $\sim n^2$ (cheia) ou $\sim n$ (malha) | grande e esparsa, diagonal dominante |

- Gauss-Seidel: $x_i = \left(b_i - \sum_{j\ne i} a_{ij}x_j\right)/a_{ii}$, usando os valores **já atualizados**.
- Relaxação: $x_i \leftarrow \lambda x_i^{novo} + (1 - \lambda)x_i^{velho}$.
- Número de condição: $\text{cond}(A) = \lVert A\rVert\,\lVert A^{-1}\rVert$. Perdem-se $\approx \log_{10}\text{cond}$ algarismos.

### Newton para sistemas

$$J(x^{(k)})\,\Delta x = -F(x^{(k)}), \qquad x^{(k+1)} = x^{(k)} + \Delta x, \qquad J_{ij} = \partial F_i/\partial x_j$$

Para 2×2, use Cramer: $D = ad - bc$, $\Delta x_1 = (r_1 d - b r_2)/D$, $\Delta x_2 = (a r_2 - r_1 c)/D$.

### Otimização 1D

- Razão áurea: $\tau = 0{,}618$, $x_1 = b - \tau(b-a)$, $x_2 = a + \tau(b-a)$. Uma avaliação nova por iteração.
- Interpolação quadrática: $d_1 = \dfrac{f_2 - f_1}{x_2 - x_1}$, $A = \dfrac{\frac{f_3 - f_2}{x_3 - x_2} - d_1}{x_3 - x_1}$, $B = d_1 - A(x_1 + x_2)$, $x^* = -B/(2A)$.
- Newton: $x_{k+1} = x_k - f'(x_k)/f''(x_k)$. Se $f'' > 0$ é mínimo.

### EDOs

- RK4: $k_1 = f(t, y)$, $k_2 = f(t + \tfrac h2, y + \tfrac h2 k_1)$, $k_3 = f(t + \tfrac h2, y + \tfrac h2 k_2)$, $k_4 = f(t + h, y + hk_3)$.
- Passo: $y_{novo} = y + \tfrac h6(k_1 + 2k_2 + 2k_3 + k_4)$. Erro $O(h^4)$: dividir $h$ por 2 divide o erro por 16.
- EDO de ordem 2 vira sistema: $y_1 = y$, $y_2 = y'$.

### Diferenças finitas (calor)

| Nó | Equação |
|---|---|
| 2D interno, permanente | $T_{ij} = \left(T_{i\pm1,j} + T_{i,j\pm1} + \dot q\,\Delta x^2/k\right)/4$ |
| 1D interno, implícito | $-Fo\,T_{m-1} + (1 + 2Fo)\,T_m - Fo\,T_{m+1} = T_m^p + A$ |
| simetria/adiabático (meia célula) | $(1 + 2Fo)\,T_1 - 2Fo\,T_2 = T_1^p + A$ |
| convecção (meia célula) | $-2Fo\,T_{N-1} + (1 + 2Fo + 2Fo\,Bi_{\Delta x})\,T_N = T_N^p + 2Fo\,Bi_{\Delta x}T_\infty + A$ |
| fluxo $q''$ entrando | lado direito do adiabático $+\,2Fo\,q''\Delta x/k$ |

Aqui $Fo = \alpha\Delta t/\Delta x^2$, $A = \dot q\Delta t/(\rho c_p)$ e **$Bi_{\Delta x} = h\Delta x/k$**, que não é o mesmo que $hL/k$.

---

## 5. Convergência e estabilidade

| Método | Converge/estável se… | Se falhar… |
|---|---|---|
| Bisseção / falsa posição | $f$ contínua e $f(x_l)f(x_u) < 0$ (sempre converge) | conferir o sinal; usar `varredura_mudanca_sinal` |
| Newton | chute perto da raiz, $f' \ne 0$ | começar com bisseção e refinar com Newton |
| Ponto fixo | $\lvert g'(x)\rvert < 1$ perto da raiz | reescrever $x = g(x)$ de outro jeito |
| Gauss-Seidel/Jacobi | diagonal **estritamente** dominante, **ou** $A$ simétrica positiva definida | reordenar as equações; sub-relaxar ($\lambda < 1$) |
| SOR (Laplace, malha $N \times N$) | $0 < \lambda < 2$; ótimo $\lambda \approx 2/(1 + \sin(\pi/N))$ | – |
| Explícito (FD) | $Fo \le 1/2$; no nó de convecção, $Fo\,(1 + Bi_{\Delta x}) \le 1/2$ | usar implícito |
| Implícito / Crank-Nicolson | sempre estáveis (o CN pode oscilar com $Fo$ grande) | – |
| RK4 | $h\,\lvert\lambda_{max}\rvert < 2{,}785$ (autovalor real) | diminuir $h$; problema rígido pede implícito |
| Euler explícito | $h\,\lvert\lambda_{max}\rvert < 2$ | – |

---

## 6. Modelagem térmica de processadores

**Caminho do calor:** junção (die) → TIM1 → IHS/spreader → TIM2 → base do dissipador → aletas → ar.

### Resistências térmicas [K/W]

| Fenômeno | Resistência |
|---|---|
| condução plana | $R = L/(kA)$ |
| condução radial (cilindro) | $R = \ln(r_2/r_1)/(2\pi k L)$ |
| convecção | $R = 1/(hA)$ |
| radiação (linearizada) | $R \approx 1/(4\varepsilon\sigma A T_m^3)$, com $T_m$ em **kelvin** |
| em série | $R = \sum R_i$ |
| em paralelo | $1/R = \sum 1/R_i$ |
| contato / pasta | $R = R''/A$ (R'' tipicamente 5 a 50 mm² K/W) |

### Balanço de energia num nó $i$

$$C_i\frac{dT_i}{dt} = P_i + \sum_j G_{ij}(T_j - T_i), \qquad G = 1/R, \qquad C = \rho\,c_p\,V$$

- Em **regime permanente** isso vira $GT = P$ (sistema linear, simétrico e positivo definido). Use Cholesky ou Gauss-Seidel.
- Com **radiação ou fuga** ($P$ depende de $T$), o sistema é não linear e pede Newton para sistemas.
- Capacitância concentrada vale se $Bi = hL_c/k < 0{,}1$, com $L_c = V/A_s$. Se valer, $T(t) = T_\infty + (T_0 - T_\infty)e^{-t/\tau}$, com $\tau = \rho c V/(hA)$.

### Potência do chip

- Dinâmica: $P_{din} \approx \alpha\,C\,V^2 f$.
- Fuga: $P_{fuga} \approx P_0\,e^{\gamma(T - T_{ref})}$, que dobra a cada ~30–40 K. É o mecanismo da fuga térmica.
- Limites típicos: $T_{j,max} \approx 95$–$105$ °C; o throttling reduz o clock perto do limite.

### Propriedades típicas

| Material | $k$ [W/(m K)] | $\rho$ [kg/m³] | $c_p$ [J/(kg K)] |
|---|---|---|---|
| silício (~80 °C) | 120–150 | 2330 | 705 |
| cobre | 390–400 | 8960 | 385 |
| alumínio | 200–240 | 2700 | 900 |
| pasta térmica | 1–10 | – | – |
| ar (convecção) | $h$ ≈ 5–25 (natural), 25–250 (forçada) | | |
| água (convecção) | $h$ ≈ 500–10 000 | | |

---

## 7. Armadilhas que custam pontos

1. **Matriz de inteiros.** `np.full((n, n), 20)` e `np.array([[1, 2], [3, 4]])` são **int**, e toda conta é truncada. Use `dtype=float` ou `20.0`. Foi esse o erro do "PPC surpresa".
2. **Radiação em °C.** $T^4$ só funciona em kelvin: $T[K] = T[°C] + 273{,}15$.
3. **Potência em Python é `**`, não `^`.** `x^2` é XOR!
4. **Índices começam em 0.** O nó 1 da apostila é `T[0]`, e `range(1, n)` não inclui `n`.
5. **Cópia x referência.** `B = A` não copia: alterar `B` altera `A`. Use `B = A.copy()`. As funções de `utils/` já trabalham em cópias.
6. **Biot da malha x Biot da placa** nas equações de FD: $h\Delta x/k$ nos nós, $hL/k$ na solução exata.
7. **Graus x radianos.** `math.sin` usa radianos; converta com `math.radians(theta)`.
8. **Tolerância sem limite de iterações** gera laço infinito. Sempre use `max_iter`.
9. **Raiz sem sentido físico.** Temperatura abaixo do ambiente, vazão negativa ou profundidade complexa indicam que é preciso mudar o chute (AP03, AP04, E06).
10. **Unidades:** mm contra m, W contra kW, h contra s. Escreva a unidade em **todo** comentário de dado.
11. **`np.arange(0, 1, 0.1)` não inclui o 1.** Prefira `np.linspace(0, 1, 11)`.
12. **Resolvedores prontos são proibidos** (`np.linalg.solve`, `np.roots`, `scipy.optimize`…). Use as funções de `utils/`.

---

## 8. Roteiro de prova em 8 passos

1. Desenhar o esquema e marcar os nós e as temperaturas conhecidas.
2. Escrever o balanço de energia de cada nó ou volume.
3. Classificar o problema e escolher o método (tabela da seção 1).
4. Copiar o exemplo mais próximo ou o [`modelo_em_branco.py`](Exemplos_Processador/modelo_em_branco.py).
5. Trocar dados e funções, mantendo `dtype=float` e as unidades.
6. Rodar e conferir o histórico: o erro está caindo?
7. **Verificar:** resíduo, balanço de energia, caso limite com solução exata, sentido físico.
8. Apresentar tabela, resultado com unidade e gráfico, e comentar em 2–3 frases.
