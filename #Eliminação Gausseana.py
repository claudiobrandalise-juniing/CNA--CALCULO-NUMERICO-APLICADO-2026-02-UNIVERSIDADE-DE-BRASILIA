#Eliminação Gausseana
import numpy as np
#Inputs

c1=11.5094
c2=11.5094
c3=19.0566
c4=16.9983
c5=11.5094

#matrizes

A= np.array([[6,0,-1,0,0],[-3,3,0,0,0],[0,-1,9,0,0],[0,-1,-8,11,-2],[-3,-1,0,0,4]],dtype=(float))
 
c=np.array([[c1],[c2],[c3],[c4],[c5]],dtype=(float))

b= np.linalg.solve(A, c)
print(b)

