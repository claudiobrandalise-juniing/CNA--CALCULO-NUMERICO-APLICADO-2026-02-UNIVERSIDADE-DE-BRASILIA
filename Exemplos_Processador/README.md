# Exemplos — Temperatura em processadores

Problemas no **formato de prova** sobre o gerenciamento térmico de um processador, cobrindo todos os grupos de métodos da disciplina. Cada script tem:

- um **enunciado** no cabeçalho;
- os **dados** com unidades;
- o **modelo** (uma função por equação física);
- a **solução** com a biblioteca [`utils/`](../utils);
- uma **verificação**: solução exata, balanço de energia ou comparação entre métodos;
- um bloco **"COMO ADAPTAR"**.

Para começar um programa novo na prova, use o [`modelo_em_branco.py`](modelo_em_branco.py).

```bash
python Exemplos_Processador/E01_temperatura_equilibrio_dissipador.py   # (e assim por diante)
```

As saídas (`.txt` e `.png`) vão para `Exemplos_Processador/outputs/`.

## Mapa dos exemplos

| Exemplo | Situação física | Método(s) | Verificação |
|---|---|---|---|
| [E01](E01_temperatura_equilibrio_dissipador.py) | Temperatura do dissipador com convecção + radiação ($T^4$); potência máxima para $T_j \le 95$ °C | bisseção, falsa posição, Newton, secante, Müller | 5 métodos concordam; resíduo ~10⁻¹¹ W |
| [E02](E02_rede_resistencias_multinucleo.py) | Rede de resistências térmicas de 4 núcleos + IHS + dissipador | Gauss com pivoteamento, LU, Cholesky, Gauss-Seidel/SOR, inversa por LU | balanço de energia (90 W = 90 W); 4 métodos concordam |
| [E03](E03_mapa_termico_spreader_2d.py) | Mapa 2D do spreader de cobre com 4 núcleos (Poisson) | Gauss-Seidel (Liebmann) com fonte e SOR | energia que sai pelas bordas = 90 W; resíduo das equações |
| [E04](E04_transiente_rc_rk4.py) | Aquecimento transiente (modelo RC de 2 nós), carga variável e *thermal throttling* | RK4 para sistema de EDOs | solução analítica (autovalores 2×2); ordem 4; limite de estabilidade |
| [E05](E05_transiente_die_implicito.py) | Pico de temperatura no die durante um "turbo" de 0,1 s | diferenças finitas implícitas + Thomas, por trechos | permanente analítico (perfil linear); explícito x implícito |
| [E06](E06_fuga_termica_newton_sistemas.py) | Fuga térmica: potência de fuga exponencial em $T_j$ | Newton 2×2 (Cramer) e 3×3 (Jacobiana numérica + Gauss) | duas soluções (estável/instável) e ponto crítico |
| [E07](E07_velocidade_otima_ventoinha.py) | Velocidade ótima da ventoinha (ventoinha x fuga), com restrição $T_j \le 85$ °C | razão áurea, interpolação quadrática; Newton e bisseção por dentro | dois otimizadores concordam; restrição ativa |
| [E08](E08_constantes_tempo_bairstow.py) | Constantes de tempo da rede RC: $\det(\lambda C + G) = 0$ | Bairstow, Müller, RK4 | Vieta; decaimento simulado = $\lambda_{lento}$ |

## Principais resultados

- **E01:** dissipador a 60,6 °C e junção a 89,1 °C. A potência máxima para $T_j \le 95$ °C é **105,4 W**. Cerca de 14 % do calor sai por radiação.
- **E02:** o núcleo 1 (30 W) é o mais quente, com 85,6 °C. Um turbo de +10 W no núcleo 2 aquece ele mesmo em 6,4 K e os vizinhos em ~5 K, porque o calor se espalha. O Gauss-Seidel leva 373 iterações, contra 115 com SOR ($\lambda = 1{,}5$).
- **E03:** o ponto mais quente fica a 77,3 °C, dentro do núcleo de 30 W. O SOR com $\lambda \approx 1{,}85$ (o ótimo teórico) cai de **1840 para 132 iterações**.
- **E04:** as constantes de tempo são 4,9 s (chip) e 132 s (dissipador). Sem throttling o chip chega a 100 °C; com throttling ele "segura" em 95 °C. O RK4 fica instável para $h > 13{,}5$ s.
- **E05:** o turbo leva a face ativa de 65,5 para 76,9 °C em ~0,1 s ($\tau \approx 26$ ms). O implícito usa 40 passos onde o explícito precisaria de ~4400, e o explícito com passo 20 % acima do limite explode.
- **E06:** com a ventoinha nominal, $T_j = 88{,}6$ °C (estável); existe outra solução, instável, a 206,7 °C. Abaixo de $hA_{crit} = 1{,}358$ W/K não há equilíbrio: é a **fuga térmica**.
- **E07:** o ótimo é $v = 2{,}95$ m/s (96,2 W no total). Com $T_j \le 85$ °C a restrição fica ativa, e o ótimo passa a $v = 3{,}71$ m/s.
- **E08:** as constantes de tempo são 0,024 s, 1,67 s e 132 s. A razão de 5400 entre a mais lenta e a mais rápida é o que se chama **problema rígido** (stiff).

## Receita para a prova

1. **Ler o problema e escrever o balanço de energia** (o que entra = o que sai + o que acumula).
2. **Classificar o problema e escolher o método:**

   | O que aparece | Método |
   |---|---|
   | uma equação em uma incógnita | raízes (E01) |
   | várias equações **lineares** | sistemas lineares (E02) |
   | uma malha 2D em regime permanente | Gauss-Seidel/SOR (E03) |
   | derivada no tempo, sem espaço | RK4 (E04) |
   | derivada no tempo e no espaço (1D) | implícito + Thomas (E05) |
   | várias equações **não lineares** | Newton para sistemas (E06) |
   | "melhor", "ótimo", "mínimo" | otimização (E07) |
   | polinômio ou constantes de tempo | Bairstow/Müller (E08) |

3. **Copiar o exemplo mais parecido** (ou o `modelo_em_branco.py`) e trocar dados e equações.
4. **Verificar sempre:** resíduo, balanço de energia, um caso limite com solução exata e o sentido físico (temperatura acima do ambiente? raiz positiva?).
5. **Imprimir:** tabela de iterações, resultado com unidade e um gráfico.

O [GUIA_DA_PROVA.md](../GUIA_DA_PROVA.md) tem as fórmulas e as armadilhas mais comuns.
