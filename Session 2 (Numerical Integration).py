import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import plotly.graph_objects as go
import sympy as smp
from sympy import *
from scipy import integrate


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

print("4 nodes: ", trapzeqd(0, 2, n=3, f=1/smp.sqrt(x**17.10+2023)))

b_list = [10, 100, 1000, 10000]

I_list = []
#for i in range(4):
    #I_list.append(trapzeqd(0, b_list[i], n=4, f=1/smp.sqrt(x**17.10+2023)))
#plt.plot(b_list, I_list)
#plt.show()

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


def trapz(x_list, y_list):
    n = len(x_list) - 1 # number of subdivisions
    I = 0
    for s in range(n):
        I += 0.5 * (y_list[s + 1] + y_list[s]) * (x_list[s + 1] - x_list[s])
    return I


# Task 4 The river Thames basin in London


t = open("Thames.txt", "r")
t_read = t.readlines()
t.close()
thames_np = np.zeros((72, 4))
thames = []
for item in t_read:
    r = item.rstrip().split(",")
    row = []
    for j in range(len(r)):
        row.append(float(r[j]))
    thames.append(row)
for i in range(72):
    for j in range(4):
        thames_np[i][j] = thames[i][j]

x_n = thames_np[:, 0]
y_n = thames_np[:, 1]
x_s = thames_np[:, 2]
y_s = thames_np[:, 3]

plt.plot(x_n, y_n)
plt.plot(x_s, y_s)
plt.axis('equal')
plt.show()

thames_area = trapz(x_n, y_n) - trapz(x_s, y_s)
print("Thames Area = ", thames_area)


# Task 5 Multiple integrals (with given analytical function): volume of the dome of
# the Royal Albert Hall

# QUESTION: how to discretise with constant dx dy when limits are not integers?

d = 0.5
h = 25
a = 67
b = 56


def z(x, y):
    f = np.sqrt(h ** 2 * (1 - x ** 2 / a ** 2 - y ** 2 / b ** 2))
    if x ** 2 / a ** 2 + y ** 2 / b ** 2 < 1:
        return f
    else:
        return 0


def G(x_i):
    x_b = math.sqrt(b**2 * (1 - x_i**2/a**2))
    x_a = - x_b
    #n = math.floor((math.floor(x_b) - math.floor(x_a)) / d)
    #return trapzeqd(x_a, x_b, n, smp.sqrt(h**2 * (1 - x_i**2/a**2 - x**2/b**2)))
    nodes = np.arange(x_a, x_b, d)
    I = 0
    for element in range(1, len(nodes)):
        I += z(x_i, nodes[element])
    I += 0.5 * (z(x_i, x_a) + z(x_i, x_b))
    I = I * d
    return I


def Rvol(x_a, x_b):
    # n = number of subintervals, n+1 = number of nodes
    n = int((x_b - x_a) / d)  #n has to be an integer
    nodes = np.linspace(x_a, x_b, n+1)
    I = 0
    for element in range(1,n):
        I += G(nodes[element])
    I += 0.5 * (G(x_a) + G(x_b))
    I = I*d
    return I


#print(Rvol(-a, a))


# scipy method (DO THIS)


f = lambda y, x: np.sqrt(h ** 2 * (1 - x ** 2 / a ** 2 - y ** 2 / b ** 2))
V = integrate.dblquad(f, -a, a, lambda x: - np.sqrt(b**2 * (1 - x**2/a**2)), lambda x: np.sqrt(b**2 * (1 - x**2/a**2)))
print("Volume = ", V[0])

xr = np.arange(-a, a, d)
yr = np.arange(-b, b, d)
Xr, Yr = np.meshgrid(xr, yr)

S = np.sqrt(h ** 2 * (1 - Xr ** 2 / a ** 2 - Yr ** 2 / b ** 2))

ax = plt.axes(projection='3d')
ax.plot_surface(Xr, Yr, S)
plt.show()

# Task 6 Multiple integrals (with given nodes): volume of an aerofoil

a = open("Aerofoil.txt", "r")
a_read = a.readlines()
a.close()
aerofoil = []
for item in a_read:
    r = item.rstrip().split(",")
    row = []
    for j in range(len(r)):
        row.append(float(r[j]))
    aerofoil.append(row)


def G_top(y_num):
    x_list, y_list = [], []
    for i in range(100):
        x_list.append(aerofoil[100*y_num+i][0])
        y_list.append(aerofoil[100*y_num+i][2])
    return trapz(x_list, y_list)


def G_bottom(y_num):
    x_list, y_list = [], []
    for i in range(100):
        x_list.append(aerofoil[100 * y_num + i][0])
        y_list.append(aerofoil[100 * y_num + i][3])
    return trapz(x_list, y_list)


def G_dif(y_num):
    return G_top(y_num) - G_bottom(y_num)


ya = []
G_l = []
for i in range(15):
    ya.append(aerofoil[100*i][1])
    G_l.append(G_dif(i))

V_aerofoil = trapz(ya, G_l)
print("Aerofoil Volume = ", V_aerofoil)


Xa = np.zeros((100, 15))
Ya = np.zeros((100, 15))
A_top = np.zeros((100, 15))
A_bottom = np.zeros((100, 15))


for i in range(100):
    for j in range(15):
        Xa[i][j] = aerofoil[100 * j + i][0]
        Ya[i][j] = aerofoil[100 * j][1]
        A_top[i][j] = aerofoil[100 * j + i][2]
        A_bottom[i][j] = aerofoil[100 * j + i][3]


ax = plt.axes(projection='3d')
ax.plot_surface(Xa, Ya, A_top)
ax.plot_surface(Xa, Ya, A_bottom)
plt.axis('equal')
plt.show()


