import math
from math import *
import numpy as np
import matplotlib.pyplot as plt


# Task A: Direct methods


def func_a2(x):
    # returns (f(x), g(x), p(x)) in form y'' + f(x)y' + g(x)y = p(x))
    return 2*x, 2, cos(3*x)


def testfunc(x):
    # returns (f(x), g(x), p(x)) in form y'' + f(x)y' + g(x)y = p(x))
    return 0, 4, 0


def myobebc_dirichlet(f, a, b, ya, yb, N):
    h = (b - a) / N
    x = np.linspace(a, b, N + 1)
    A = np.zeros((N+1, N+1))   # coefficient matrix
    p = np.zeros((N+1, 1))   # N+1 X 1 matrix with the solutions and p's
    A[0][0], A[-1][-1] = 1, 1
    p[0], p[-1] = ya, yb
    for i in range(1, N):
        ai = 1 / h ** 2 - f(x[i])[0]
        bi = f(x[i])[1] - 2 / h**2
        ci = 1 / h**2 + f(x[i])[0]/(2*h)
        abci = [ai, bi, ci]
        A[i][i-1:i+2] = abci[0:3]
        p[i] = f(x[i])[2]
    y = np.matmul(np.linalg.inv(A), p)
    return x, y.flatten()


'''
test = myobebc_dirichlet(testfunc, 0, pi/4, -2, 10, 50)
A_2 = myobebc_dirichlet(func_a2, 0, pi, 1.5, 0, 100)

plt.plot(test[0], test[1])
plt.show()
'''


# Task B: Types of boundary conditions


def func_2b(x):
    # returns (f(x), g(x), p(x)) in form y'' + f(x)y' + g(x)y = p(x))
    return x, 1, 5*x


def testfunc2(x):
    # returns (f(x), g(x), p(x)) in form y'' + f(x)y' + g(x)y = p(x))
    return 0, 2, 0


def myobebc(f, a, b, bca, bcb, N, R):
    h = (b - a) / N
    x = np.linspace(a, b, N + 1)
    A = np.zeros((N+1, N+1))   # coefficient matrix
    p = np.zeros((N+1, 1))   # N+1 X 1 matrix with the solutions and p's

    # r0, r1 specify bc type at x = a, r2, r3 specific bc type at x = b
    # Dirichlet: r_even = 0, r_odd = 1
    # Neumann: r_even = 1, r_odd = 0
    # Robin (Mixed): r_even, r_odd != 0
    A[0][0] = R[1] - R[0] / h      # modifications to coefficient matrix for general boundary conditions
    A[0][1] = R[0] / h
    A[-1][-2] = - R[2] / h
    A[-1][-1] = R[3] + R[2] / h

    p[0], p[-1] = bca, bcb
    for i in range(1, N):
        ai = 1 / (h ** 2) - f(x[i])[0] / (2 * h)
        bi = f(x[i])[1] - 2 / (h ** 2)
        ci = 1 / (h ** 2) + f(x[i])[0] / (2 * h)
        #abci = [ai, bi, ci]
        #A[i][i-1:i+2] = abci[0:3]
        A[i][i-1] = ai
        A[i][i] = bi
        A[i][i+1] = ci
        p[i] = f(x[i])[2]
    y = np.matmul(np.linalg.inv(A), p)
    return x, y.flatten()


'''
test2 = myobebc(testfunc2, 0, pi/4, -2, 10, 50, [0, 1, 0, 1])
B_2 = myobebc(func_a2, 0, 2, 0, 5, 1000, [1, 0, 0, 1])
B_3 = myobebc(func_a2, 0, 2, 5, 0, 1000, [0, 1, 1, 0])

plt.plot(B_2[0], B_2[1])
plt.show()

plt.plot(B_3[0], B_3[1])
plt.show()
'''


# Task C: Heat transfer in a nuclear fuel rod

h = 6 * 10 ** 4
k = 16.75
R = 15
w = 3
T_w = 473


def temp_ode(x):
    # returns (f(x), g(x), p(x)) in form y'' + f(x)y' + g(x)y = p(x))
    return k/x, 0, -10 ** 8 * e ** (-x/R) / x


#'''
C_1 = myobebc(temp_ode, R, R+w, - 6.32 * 10 ** 5 / k, h * T_w / k, 1000, [0, 1, 1, h/k])
#print(C_1)

plt.plot(C_1[0], C_1[1])
plt.xlabel("radius (mm)")
plt.ylabel("temp (K)")
plt.show()
#'''
