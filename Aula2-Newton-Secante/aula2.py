
#bibliotecas utilizadas
import numpy as np
import matplotlib.pyplot as plt

#Inputs

A= 0.01                 #Area em m**2
h=10                    # W/(m**2K)
Ep=0.8                  
T8=300                  #temperatura em K
P=20                    #Watts
Sigma= 5.7*10**-8       #W/((m**2)*(K**2))
tol= 10**-6
T0= 350                 #chute inicial
K= 100                  #iterações



##PARTE 1 METODO NEWTON RAPHSON:
def t(T):
    return h*A*(T-T8)+ Ep*Sigma*A*((T**4)-(T8**4))-P     #Função proposta

def dt(T):
    return h*A + 4*Ep*Sigma*A*T**3                      #Derivada da  função proposta

vetRes = []
it = 0
err = 1

while(err>tol and it < K):
    
 vetRes.append(Ts) 
 it 







