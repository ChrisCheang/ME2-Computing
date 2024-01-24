import math
from math import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


# Task A: Heat conduction in a bar


fastforward = 10


def solve_heat(a, b, Ta, Tb, T0, dx, dt, t_end):
    # Initialisation and initial conditions
    alpha = 1.172 * (10 ** -5)
    nx = int((b - a) / dx) + 1
    nt = int(t_end / dt) + 1
    T = np.ndarray((nt, nx))
    T[0] = T0
    T[:,0] = Ta
    T[:,-1] = Tb

    # Calculation
    for t in range(1, nt):
        for x in range(1, nx - 1):
            T[t][x] = (alpha*dt/(dx**2)) * (T[t-1][x+1]+T[t-1][x-1]) + ((1 - 2*alpha*dt/dx**2) * T[t-1][x])

    return T


'''
a = 0
b = 0.5
dx = 0.01
dt = 1
t_end = 3600
T0 = [10 for i in range(int((b-a)/dx) + 1)]

X = np.arange(a, b+dx, dx)

A1 = solve_heat(0, 0.5, 50, 50, T0, dx, dt, 3600)




for i in range(int((1/fastforward) * t_end / dt)):
    plt.plot(X, A1[fastforward * i])
    plt.pause(0.01)
    plt.clf()

plt.plot(X, A1[-1])
plt.show()


A3 = solve_heat(0, 0.5, 50, 70, T0, dx, dt, 3600)

for i in range(int((1/fastforward) * t_end / dt)):
    plt.plot(X, A3[fastforward * i])
    plt.pause(0.01)
    plt.clf()

plt.plot(X, A3[-1])
plt.show()

print("Q1 ans = ", A3[-1,25])
'''


# Task B: Heat conduction in a bar with a heat source


def solve_heat_p2(a, b, Tw, Ts, T0, dx, dt, t_end):
    # Initialisation and initial conditions
    alpha = 1.172 * (10 ** -5)
    k = 40
    h = 500
    nx = int((b - a) / dx) + 1
    nt = int(t_end / dt) + 1
    T = np.ndarray((nt, nx))
    T[0] = T0
    middle = int((nx+1)/2) - 1 # setting middle temp, the +1 is only needed for even subdivisions of t.
    T[:,middle] = Ts

    # Calculation
    for t in range(1, nt):
        kxh = k/dx + h
        T[t][0] = ((k/dx) * T[t-1][1] + h * Tw)/kxh
        T[t][-1] = ((k/dx) * T[t-1][-2] + h * Tw)/kxh
        for x in [i for i in range(1, nx - 1) if i != middle]:
            T[t][x] = (alpha*dt/(dx**2)) * (T[t-1][x+1]+T[t-1][x-1]) + ((1 - 2*alpha*dt/dx**2) * T[t-1][x])

    return T


'''
a = 0
b = 0.5
dx = 0.01
dt = 1
t_end = 1200
T0 = [10 for i in range(int((b-a)/dx) + 1)]

X = np.arange(a, b+dx, dx)

B2 = solve_heat_p2(a, b, 5, 100, T0, dx, dt, t_end)


for i in range(int((1/fastforward) * t_end / dt) + 1):
    plt.plot(X, B2[fastforward * i])
    plt.pause(0.01)
    plt.clf()

plt.show()

plt.plot(X, B2[0])
plt.plot(X, B2[500])
plt.plot(X, B2[-1])

plt.show()

print("Q2 Target (time index is ans) = ", B2[232, 34])
'''


# Task C: Stability of the finite difference numerical method


def Courant(alpha, dt, dx):
    x = alpha * dt / dx ** 2
    if x < 0.5:
        return True
    else:
        return False


'''
# change case number here
case_no = 1

a = 0
b = 0.5
dx = [0.01, 0.05, 0.001, 0.001][case_no - 1]
dt = [1, 1, 1, 0.04][case_no - 1]
t_end = 3600
alpha = 1.172 * (10 ** -5)
T0 = [10 for i in range(int((b-a)/dx) + 1)]

X = np.arange(a, b+dx, dx)

C1i = solve_heat(0, 0.5, 50, 50, T0, dx, dt, 3600)

courant = Courant(alpha, dt, dx)
print(courant)

plt.plot(X, C1i[-1])
plt.show()
'''


# Task D: PDE with multiple spatial dimensions: baking a freaking potato in the oven

dx = 0.01
dy = 0.01
dt = 1
t_end = 1200

T_air = 25
T_wall = 180
T_potato = -15


def po_ta_toes_boil_it_mash_it_stick_it_in_a_stew(xb, yb, xap, xbp, yap, ybp):
    def a(x,y):
        if xap < x < xbp and ybp > y > yap:
            return 1.3 * 10 ** (-7)
        else:
            return 1.9 * 10 ** (-5)

    nx = int(xb / dx) + 1
    ny = int(yb / dx) + 1
    nt = int(t_end / dt) + 1
    T = np.ndarray((nt, nx, ny))

    T[0,:,:] = T_air
    T[0, int(xap / dx):int(xbp / dx), int(yap / dy):int(ybp / dy)] = T_potato
    T[:,:,0], T[:,:,-1] = T_wall, T_wall
    T[:,0,:], T[:,-1,:] = T_wall, T_wall

    for t in range(1, nt):
        for x in range(1, nx - 1):
            for y in range(1, ny - 1):
                adtdx2 = a(x,y) * dt / (dx ** 2)
                adtdy2 = a(x,y) * dt / (dy ** 2)
                oneminus = 1 - 2*adtdx2 - 2*adtdy2
                T[t,x,y] = adtdx2 * (T[t-1,x+1,y] + T[t-1,x-1,y]) + adtdy2 * (T[t-1,x,y+1] + T[t-1,x,y-1]) + oneminus * T[t-1,x,y]

    return T


xb = 0.4
yb = 0.4
xap = 0.15
xbp = 0.25
yap = 0.15
ybp = 0.25

potato = po_ta_toes_boil_it_mash_it_stick_it_in_a_stew(xb, yb, xap, xbp, yap, ybp)

norm = matplotlib.colors.Normalize(vmin=T_potato, vmax=T_wall)

for i in range(int((1/fastforward) * t_end / dt) + 1):
    plt.imshow(potato[fastforward*i], interpolation='bilinear', norm=norm)
    plt.colorbar()
    plt.pause(0.01)
    plt.clf()

plt.imshow(potato[-1], interpolation='bilinear', norm=norm)
plt.colorbar()
plt.show()
