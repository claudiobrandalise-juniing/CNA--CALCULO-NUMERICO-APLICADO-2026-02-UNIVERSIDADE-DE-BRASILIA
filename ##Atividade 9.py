##Atividade 9

from asyncio import log

import numpy as np
import matplotlib.pyplot as plt

#definição de parametros

r1=0.040
L=20
Ti=180
Tinf=40
k=0.045
h=12
t=30000
Pe=0.7
Ci=9000
e=float(input("Qual valor de e?"))


def Rcond(e):
    return (np.log((r1+e)/r1)/(2*np.pi*k))

def Rconv(e):
    return 1/(2*np.pi*h*(r1+e))

def Qlinha(e):
    return (Ti-Tinf)/(Rcond(e)+Rconv(e))

Qe=L*Qlinha(e)
print(Qe)
print(Qlinha(e))



        
    
    