#PPC SURPRESA

import numpy as np
import matplotlib.pyplot as plt
#definição de parametros

n=5  #numero de nós para cada lado
L=1  #comprimento na placa
tamanho = n + 2  #numero de posições da matriz 

def atualizar_ponto(T, i, j):
    
    return (T[i-1, j] + T[i+1, j] + T[i, j-1] + T[i, j+1]) / 4

iter = 0
iter_max=100000
tol = 1e-6
erro = 1500

T = np.full((tamanho, tamanho), 20) # Inicializando a matriz T com 20 °C em todas as posições
T[n+1, 1:n+1] = 100                 #Atualizando a base com valor de 100 °C
print(T)

while erro > tol and iter < iter_max:
        erro = 0.0
        
        for j in range(1, n + 1):
            for i in range(1, n + 1):
            
                antigo = T[i][j]
            
                T[i][j] = atualizar_ponto(T, i, j)
            
                diferenca = abs(T[i][j] - antigo)
                erro = max(erro, diferenca)
            
        iter += 1


print(f"\nConvergência atingida após {iter} iterações.")
print(f"Erro final: {erro:.5e}")
print("")

print(f"\ntemperatura no ponto [4,4]: {T[3, 3]}")  # Imprimindo o valor da temperatura no ponto (4,4)
print("Matriz T final:")
print("")


print(T)

##Plotando a matriz de temperatura
plt.imshow(T, cmap='rainbow', interpolation='nearest') #testado com "hot", "coolwarm", "rainbow", "viridis", "plasma", "inferno", "magma", "cividis"
plt.colorbar(label='Temperatura (°C)')
plt.title('Distribuição de temperatura na placa')
plt.xlabel('Em x/i')
plt.ylabel('Em y/j')
plt.show()

#devido ao curto tempo, não foi possivel implementar o codigo de comparação para n = 10 e n=20 para
# comprara o numero de iterções e o erro final, mas a implementação do código para n=10 e n=20 seria
# similar, apenas alterando o valor de n e ajustando a matriz T de acordo. para n=5 as iterações foram de
#14, para n=10 foi de 24 e para n=20 foi de 24 também, para plotar a temperatura no meio da placa em qualquer
# n, basta acessar o ponto T[n//2, n//2] e imprimir o valor.

#as temperaturas na borda não são atualizadas, pois são condições de contorno fixas. A temperatura na borda 
# superior é de 20 °C, enquanto a borda inferior é de 100 °C. As bordas laterais também são mantidas em 20 °C,
#o numero de iterações muda quando n é maior pois a precisão desejada é mantida constante, mas o 
# número de pontos na malha aumenta, o que pode levar a mais iterações para atingir a convergência.