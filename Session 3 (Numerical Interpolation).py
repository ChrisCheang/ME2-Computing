import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import plotly.graph_objects as go


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


xl = np.arange(0, 3+0.05, 0.05)
yl= []
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
plt.plot(xl, yl1)
plt.plot(xl, yl2)
plt.plot(xl, yl3)
print(xl)

plt.show()
