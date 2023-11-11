import math
import numpy as np
import matplotlib.pyplot as plt
import sympy as smp
from sympy import *
from mpl_toolkits import mplot3d
import plotly.graph_objects as go


# Task A: Lagrangian polynomials and interpolation

def lagrangian(j, xp, nodes):
    Lj = 1
    for i in [i for i in range(len(nodes)) if i != j]:
        Lj *= (xp - nodes[i]) / (nodes[j] - nodes[i])
    return Lj


def lagrInterp(xn, yn, x):
    pn = 0
    for i in range(len(xn)):
        pn += yn[i] * lagrangian(i, x, xn)
    return pn


#'''

xl = np.arange(0, 3 + 0.05, 0.05)
yl = []
for i in range(len(xl)):
    yl.append(math.sin(xl[i]))

xn1 = np.linspace(1, 2, 2)
yn1 = []
yl1 = []
for i in range(len(xn1)):
    yn1.append(math.sin(xn1[i]))
for i in range(len(xl)):
    yl1.append(lagrInterp(xn1, yn1, xl[i]))

xn2 = np.linspace(1, 2, 3)
yn2 = []
yl2 = []
for i in range(len(xn2)):
    yn2.append(math.sin(xn2[i]))
for i in range(len(xl)):
    yl2.append(lagrInterp(xn2, yn2, xl[i]))

xn3 = np.linspace(1, 2, 4)
yn3 = []
yl3 = []
for i in range(len(xn3)):
    yn3.append(math.sin(xn3[i]))
for i in range(len(xl)):
    yl3.append(lagrInterp(xn3, yn3, xl[i]))


plt.plot(xl, yl)
#plt.plot(xl, yl1)
#plt.plot(xl, yl2)
plt.plot(xl, yl3)

print("Cubic interp at x = 0.8 -->", lagrInterp(xn3, yn3, 0.8))

#plt.show()

#'''
#'''


# Task B: Newton interpolation


def NewtDivDiff(xn, yn):
    l = len(xn)
    if l == 1:
        return yn[0]
    else:
        left = NewtDivDiff(xn[1:l], yn[1:l])
        right = NewtDivDiff(xn[0:l - 1], yn[0:l - 1])
        return (left - right) / (xn[-1] - xn[0])


def NewtonInterp(xn, yn, x):
    k = len(xn) - 1
    p = yn[0]
    for i in range(1, k + 1):
        pi = 1
        for j in range(0, i):
            pi = pi * (x - xn[j])
        p += NewtDivDiff(xn[0:i + 1],
                         yn[0:i + 1]) * pi  # slice the first i+1 items i.e. 0 to i, i+1 excluded in slicing
    return p


'''

xl = np.arange(0, 3 + 0.05, 0.05)
yl = []
for i in range(len(xl)):
    yl.append(math.sin(xl[i]))

xn1 = np.linspace(1, 2, 2)
yn1 = []
yl1 = []
for i in range(len(xn1)):
    yn1.append(math.sin(xn1[i]))
for i in range(len(xl)):
    yl1.append(NewtonInterp(xn1, yn1, xl[i]))

xn2 = np.linspace(1, 2, 3)
yn2 = []
yl2 = []
for i in range(len(xn2)):
    yn2.append(math.sin(xn2[i]))
for i in range(len(xl)):
    yl2.append(NewtonInterp(xn2, yn2, xl[i]))

xn3 = np.linspace(1, 2, 3)
yn3 = []
yl3 = []
for i in range(len(xn3)):
    yn3.append(math.sin(xn3[i]))
for i in range(len(xl)):
    yl3.append(NewtonInterp(xn3, yn3, xl[i]))

plt.plot(xl, yl)
plt.plot(xl, yl1)
plt.plot(xl, yl2)
plt.plot(xl, yl3)
plt.show()

'''
#'''

xl = np.arange(-1, 1 + 0.02, 0.02)
yl = []
for i in range(len(xl)):
    yl.append(1 / (1 + 25 * xl[i] ** 2))

#'''
#'''

n = 11 # order
xn1 = np.linspace(-1, 1, n+1)
yn1 = []
yl1 = []
for i in range(len(xn1)):
    yn1.append(1/(1+25*xn1[i]**2))
for i in range(len(xl)):
    yl1.append(NewtonInterp(xn1, yn1, xl[i]))

#print(xn1)
#print(yn1)


print("Newton interp at x = 0.8 -->", NewtonInterp(xn1, yn1, 0.8))

plt.plot(xl, yl)
plt.plot(xl, yl1)
#plt.show()

#'''
'''

xn = np.zeros((15, 14))
yn = np.zeros((15, 14))
ynl = np.zeros((15, 101))

for i in range(0, 14):
    for j in range(i+1):
        xn[i][j] = np.linspace(-1, 1, i+1)[j]
        yn[i][j] = (1 / (1 + 25 * xn[i][j] ** 2))
    for j in range(len(xl)):
        ynl[i][j] = NewtonInterp(xn[i][0:i+1], yn[i][0:i+1], xl[j])

'''
'''

plt.plot(xl, yl)
for i in range(1,14):
    plt.plot(xl, ynl[i])
#plt.show()

'''
#'''

# Task C: Splines

x = symbols("x", real=True)
f = 1 / (1 + 25 * x ** 2)
fprime = diff(f, x)
a = fprime.subs(x, -1).evalf()
b = fprime.subs(x, -1).evalf()

xn = np.linspace(-1, 1, 11)
yn = []
for i in range(len(xn)):
    yn.append(f.subs(x, xn[i]).evalf())

#'''

def Splines(x, xn, yn, a, b):
    n = len(xn)
    h = xn[1] - xn[0]
    A = np.zeros((n, n))
    d = np.zeros(n)

    A[0][0] = 1
    A[-1][-1] = 1
    d[0] = a
    d[-1] = b
    for i in range(1, n - 1):
        A[i][i - 1:i + 2] = [1, 4, 1]
        d[i] = 3 * (yn[i + 1] - yn[i - 1]) / h
    A = np.linalg.inv(A)
    v = np.matmul(A, d)

    if x < xn[-1]:  # < instead of != to deal with floating point error (1.000000018 instead of 1)
        j = floor((x - xn[0]) / h)  # j = range index, indicates between which two nodes the point x lies
    else:
        j = n - 2

    acoeff = yn[j]
    bcoeff = v[j]
    ccoeff = 3 * (yn[j + 1] - yn[j]) / (xn[j + 1] - xn[j]) ** 2 - (v[j + 1] + 2 * v[j]) / (xn[j + 1] - xn[j])
    dcoeff = -2 * (yn[j + 1] - yn[j]) / (xn[j + 1] - xn[j]) ** 3 + (v[j + 1] + v[j]) / (xn[j + 1] - xn[j]) ** 2
    q = acoeff + bcoeff * (x - xn[j]) + ccoeff * (x - xn[j]) ** 2 + dcoeff * (x - xn[j]) ** 3
    return q

#'''

ynl = []
for i in range(len(xl)):
    ynl.append(Splines(xl[i], xn, yn, a, b))

print("Splines interp at x = 0.8 -->", Splines(0.8, xn, yn, a, b))

plt.plot(xl, yl)
plt.plot(xl, ynl)
#plt.show()

#'''
