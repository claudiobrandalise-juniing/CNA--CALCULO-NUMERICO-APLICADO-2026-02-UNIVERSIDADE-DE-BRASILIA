# Aulas práticas de programação

Atividades feitas em sala (laboratório LCCC). Cada pasta tem:

- um `main.py`, que segue o enunciado item por item;
- um `README.md` com o modelo, o algoritmo, o dicionário de variáveis, a validação e **as respostas das questões de discussão**;
- uma pasta `outputs/` com os resultados (`resultados.txt`) e os gráficos.

| Pasta | Método | Problema físico |
|---|---|---|
| [AP01_Bissecao_FalsaPosicao](AP01_Bissecao_FalsaPosicao) | bisseção, falsa posição (e a modificada) | profundidade normal em canal retangular (Manning) |
| [AP02_Newton_Secante](AP02_Newton_Secante) | Newton-Raphson, secante | temperatura de equilíbrio com convecção + radiação |
| [AP03_Muller_Canal_Critico](AP03_Muller_Canal_Critico) | Müller (polinômio de 6ª ordem) | profundidade crítica em canal trapezoidal |
| [AP04_Newton_Sistemas_Tubulacoes](AP04_Newton_Sistemas_Tubulacoes) | Newton para sistemas (Cramer 2×2) | divisão de vazão em ramos paralelos |
| [AP05_Gauss_Ingenua_Reatores](AP05_Gauss_Ingenua_Reatores) | eliminação gaussiana ingênua + FLOPs | rede de 5 reatores de mistura |
| [AP06_Pivoteamento_Barras](AP06_Pivoteamento_Barras) | pivoteamento parcial, condicionamento | treliça com barras quase paralelas |
| [AP07_Inversao_LU_Molas](AP07_Inversao_LU_Molas) | inversão por decomposição LU | matriz de flexibilidade de molas |
| [AP08_Gauss_Seidel_Placa](AP08_Gauss_Seidel_Placa) | Gauss-Seidel (Liebmann) | condução 2D numa placa |
| [AP09_Otimizacao_Isolamento](AP09_Otimizacao_Isolamento) | razão áurea, interpolação quadrática | espessura ótima de isolamento térmico |

Para executar qualquer uma, a partir da raiz do repositório:

```bash
python Aulas_Praticas/AP01_Bissecao_FalsaPosicao/main.py
```
