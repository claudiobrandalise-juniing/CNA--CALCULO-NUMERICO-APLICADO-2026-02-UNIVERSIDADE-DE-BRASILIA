#bibliotecas utilizadas
import numpy as np
import matplotlib.pyplot as plt


# exercicio: escreva um programa de computador que calcule, utilizando o metodo da bissecção a raiz da seguinte equação algébrica:

##INPUTS
a=0
b=1.5
m= (a+b)/2
tol= 10**(-4)
Ft=1

def f(x):
    f=((0.47140*((20*x)**(5/3)))/((20+2*x)**(2/3))) -5 
    return f

# vetX = np.arange(0,1.5, step = 0.1)
# plt.plot(vetX,f(vetX))
# plt.savefig("teste1.png")

Fa= f(a)
Fb= f(b)
Fm= f(m)

while abs(Fm) >= tol:
    m= (a+b)/2
    if Fm


