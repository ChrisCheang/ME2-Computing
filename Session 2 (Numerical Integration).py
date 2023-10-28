import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import plotly.graph_objects as go
import sympy as smp
from sympy import *


# Task 1 Trapezium rule for functions with equidistant nodes


x = symbols("x", real=True)
y = symbols("y", real=True)


def trapzeqd(x_a, x_b, n, f):
    # n = number of subintervals, n+1 = number of nodes
    h = (x_b - x_a) / n   #n has to be an integer
    func = smp.sympify(f)
    nodes = np.linspace(x_a, x_b, n+1)
    I = 0
    for element in range(1,n):
        I += func.subs(x, nodes[element]).evalf()
    I += 0.5 * (func.subs(x, x_a).evalf() + func.subs(x, x_b).evalf())
    I = I*h
    return I


#print("5 nodes: ", trapzeqd(0, 2, n=4, f=1/smp.sqrt(x**17.10+2023)))
#print("11 nodes: ", trapzeqd(0, 2, n=10, f=1/smp.sqrt(x**17.10+2023)))

b_list = [10, 100, 1000, 10000]

I_list = []
#for i in range(4):
    #I_list.append(trapzeqd(0, b_list[i], n=4, f=1/smp.sqrt(x**17.10+2023)))
#plt.plot(b_list, I_list)
#.show()

I_list2 = []
#for i in range(4):
    #I_list2.append(trapzeqd(0, b_list[i], n=int(b_list[i]/0.5), f=1/smp.sqrt(x**17.10+2023)))
#plt.plot(b_list, I_list2)
#plt.show()


# Task 2 Numerical integration of diverging improper integrals


b_list = [10, 100, 1000, 10000]
f2 = 1/smp.sqrt(x**1.10+2023)

I_list = []
#for i in range(4):
    #I_list.append(trapzeqd(0, b_list[i], n=4, f=f2))
#plt.plot(b_list, I_list)
#plt.show()

I_list2 = []
#for i in range(4):
    #I_list2.append(trapzeqd(0, b_list[i], n=int(b_list[i]/0.5), f=f2))
#plt.plot(b_list, I_list2)
#plt.show()


# Task 3 Trapezium rule for functions with non-equidistant nodes


def trapz(x_list, f):
    n = len(x_list) - 1 # number of subdivisions
    func = smp.sympify(f)
    I = 0
    for s in range(n):
        y_n = func.subs(x, x_list[s]).evalf()
        y_np1 = func.subs(x, x_list[s + 1]).evalf()
        I += 0.5 * (y_n + y_np1) * (x_list[s + 1] - x_list[s])
    return I


#print(trapz([0, 0.5, 1, 1.5, 2], f=1/smp.sqrt(x**17.10+2023)))


# Task 4 The river Thames basin in London


