#from math import *
import numpy as np
import sympy as sym
from sympy import *
from sympy import diff

'''
fceSwitch == 0 -> Matyas function x[0.5,0.5], N=50
fceSwitch == 1 -> Eggholder function x[100,100], N=4
'''
fceSwitch=1

'''
konjug.grad
'''

x=sym.Symbol("x")
y=sym.Symbol("y")

if fceSwitch==0:
    xzad = [0.5,0.5]
    N=50
if fceSwitch==1:
    xzad = [100,100]
    N=2


def grad(xz1):
    if fceSwitch==0:
        fce=(0.26*(x**2+y**2)-0.48*x*y)
    if fceSwitch==1:
        fce=(-(y+47)*sin(sqrt(abs((x/2)+(y+47))))-x*sin(sqrt(abs(x-(y+47)))))
    xzad=xz1
    
    ga=diff(fce,x)
    gb=diff(fce,y)
    gaa=ga.subs(x,xzad[0])
    gaa2=gaa.subs(y,xzad[0])
    
    gbb=gb.subs(x,xzad[1])
    gbb2=gbb.subs(y,xzad[1])
    
    grad=[gaa2,gbb2]
    return grad
#print(grad(xzad))

#lambda odhad
lamb=0.3

matrZ=[[-1.0,0.0],[0.0,-1.0]]


x1=np.zeros((N+1,2))
x1=x1.tolist()
x1[0]=xzad

d_xc=np.zeros((N+1,2))
xc=np.zeros((N+1,2))

xc=xc.tolist()
for i in range(N):
    xc[i]=(grad(x1[i]))
    
    if i==0: #init
        H=0
        d_xc=(d_xc.tolist())
        d_xc[i]=(-np.dot(lamb,grad(x1[i]))) #dot arr
    else:
        H=(np.dot(xc[i],xc[i])/np.dot(xc[i-1],xc[i-1]))
        d_xc[i]=(lamb*np.dot(matrZ,grad(x1[i]))+H*d_xc[i-1])

    x1[i+1]=x1[i]+d_xc[i] 
    
    x1=np.array(x1)   
print(x1[:,0],'\n')
print(x1[:,1],'\n')

