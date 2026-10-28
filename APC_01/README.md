# APC 01 — Adimensionalização e planejamento do RK4

> As APCs são entregues em PDF (LaTeX). Este arquivo guarda a **dedução** e o **planejamento** para consulta. A implementação está no [PPC 01](../PPC_01_RungeKutta_Sedimentacao/README.md).

## 1. Adimensionalização da equação do movimento

A equação dimensional é:

$$m_p\frac{dv_z}{dt} = -6\pi\eta a v_z - \frac94\pi\rho_f a^2 v_z^2 + \frac{4\pi a^3}{3}\Delta\rho\,g$$

As escalas usadas são a velocidade de Stokes $U_s$ e o tempo convectivo $a/U_s$:

$$v_z = U_s v^*, \qquad t = \frac{a}{U_s}t^*, \qquad U_s = \frac{2\Delta\rho\,g\,a^2}{9\eta}$$

$U_s$ é a velocidade em que o empuxo líquido equilibra o arrasto de Stokes.

Substituindo, como $\frac{dv_z}{dt} = \frac{U_s^2}{a}\frac{dv^*}{dt^*}$, fica:

$$m_p\frac{U_s^2}{a}\frac{dv^*}{dt^*} = -6\pi\eta a U_s v^* - \frac94\pi\rho_f a^2 U_s^2 {v^*}^2 + \frac{4\pi a^3}{3}\Delta\rho\,g$$

Dividindo tudo por $6\pi\eta a U_s$ (a escala da força de Stokes), cada termo vira um grupo adimensional:

| Termo | Resultado |
|---|---|
| inercial | $\dfrac{m_p U_s}{6\pi\eta a^2}\dfrac{dv^*}{dt^*} = St\,\dfrac{dv^*}{dt^*}$ |
| Stokes | $-v^*$ |
| Oseen | $-\dfrac{9\pi\rho_f a^2 U_s^2}{4\cdot6\pi\eta a U_s}{v^*}^2 = -\dfrac38\dfrac{\rho_f U_s a}{\eta}{v^*}^2 = -\dfrac38 Re_s{v^*}^2$ |
| empuxo | $\dfrac{4\pi a^3\Delta\rho g/3}{6\pi\eta a U_s} = \dfrac{2a^2\Delta\rho g}{9\eta U_s} = 1$ (pela definição de $U_s$) |

O resultado é:

$$\boxed{St\,\frac{dv^*}{dt^*} = 1 - v^* - \frac38 Re_s\,{v^*}^2}, \qquad St = \frac{m_p U_s}{6\pi\eta a^2}, \qquad Re_s = \frac{\rho_f U_s a}{\eta}$$

A forma adimensional depende só de $St$ (inércia da partícula frente ao arrasto viscoso) e de $Re_s$ (inércia do fluido). Para $Re_s \to 0$, ela volta à equação de Stokes $St\,\dot v^* = 1 - v^*$.

**Verificação da solução de Riccati.** Faça $v = v_\infty + 1/w$, em que $v_\infty$ zera o lado direito. A equação vira uma EDO linear em $w$:

$$w' = -Pw + Q, \qquad P = -\frac{1}{St}\left(1 + \frac34 Re\,v_\infty\right), \qquad Q = \frac{3Re}{8St}$$

A solução é $w = Q/P + (w_0 - Q/P)e^{-Pt}$, com $w_0 = -1/v_\infty$, porque $v(0) = 0$. Isso reproduz a expressão do enunciado.

## 2. Planejamento do algoritmo

**Variáveis principais:**

- parâmetros físicos: $St$ e $Re$;
- condição inicial: $v_0 = 0$;
- discretização: passo $h$, instante final $t_f$ e número de passos $N = t_f/h$;
- vetores de saída: $t[0..N]$ e $v[0..N]$.

**Fluxo lógico:**

```
ler parâmetros (St, Re, v0, h, t_final)            <- inputs/parametros.json
definir f(t, v) = (1 - v - 3/8 Re v^2) / St
N = t_final / h ;  t[0] = 0 ;  v[0] = v0
para i = 0 .. N-1:                                 <- avanço no tempo
    k1 = f(t[i],       v[i])
    k2 = f(t[i] + h/2, v[i] + h k1/2)
    k3 = f(t[i] + h/2, v[i] + h k2/2)
    k4 = f(t[i] + h,   v[i] + h k3)
    v[i+1] = v[i] + h (k1 + 2k2 + 2k3 + k4)/6
    t[i+1] = t[i] + h
calcular a solução exata nos mesmos t[i]            <- validação
erro = |v - v_exata| ;  RMSE
salvar tabelas (.txt) e gráficos (.png)            <- outputs/
```

**Entradas:** arquivo JSON com $St$, $Re$, $h$ e $t_f$, para mudar os parâmetros sem tocar no código.

**Armazenamento:** vetores pré-alocados com $N+1$ posições (não é preciso guardar $k_1..k_4$).

**Saídas:** gráficos de $v(t)$, do erro em função de $h$ (log-log, inclinação 4) e do efeito de $Re$, mais um arquivo `.txt` com as tabelas de erro.
