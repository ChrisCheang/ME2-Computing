import math
import numpy as np
import matplotlib.pyplot as plt
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
#plt.plot(xl, yl3)

#plt.show()

'''


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
        p += NewtDivDiff(xn[0:i + 1], yn[0:i + 1]) * pi  # slice the first i+1 items i.e. 0 to i, i+1 excluded in slicing
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

xl = np.arange(-1, 1 + 0.02, 0.02)
yl = []
for i in range(len(xl)):
    yl.append(1/(1+25*xl[i]**2))

'''
xn1 = np.linspace(-1, 1, 7)
yn1 = []
yl1 = []
for i in range(len(xn1)):
    yn1.append(1/(1+25*xn1[i]**2))
for i in range(len(xl)):
    yl1.append(NewtonInterp(xn1, yn1, xl[i]))

print(xn1)
print(yn1)

plt.plot(xl, yl)
plt.plot(xl, yl1)
plt.show()
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

#'''

plt.plot(xl, yl)
for i in range(1,14):
    plt.plot(xl, ynl[i])
plt.show()

#'''
