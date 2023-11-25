import math
import numpy as np
import matplotlib.pyplot as plt
import sympy as smp
from sympy import *

# from Session_3_Numerical_Integration import lagrange

x = symbols("x", real=True)


# Task A: Simpson Integration adn Adaptive Simpson Integration


def Simpson(fx, a, b, dx):
    f = smp.sympify(fx)
    xn = np.arange(a, b + dx, dx)
    I = 0
    for i in range(len(xn)):
        yi = f.subs(x, xn[i]).evalf()
        if i == 0 or i == len(xn) - 1:
            I += yi
        elif (i % 2) == 0:  # even
            I += 2 * yi
        else:
            I += 4 * yi
    h = (b - a) / (len(xn) - 1)
    return I * h / 3


def SimpsonAdaptive(fx, a, b, error):
    dx = (b - a) / 2
    eps = Simpson(fx, a, b, (b - a) / 4) - Simpson(fx, a, b, dx)  # first error estimation
    while abs(eps) > error:
        Ihalf = Simpson(fx, a, b, dx / 2)
        Ih = Simpson(fx, a, b, dx)
        eps = (Ihalf - Ih) / 15
        dx = dx / 2
    return Simpson(fx, a, b, dx)


'''
test2 = SimpsonAdaptive(1/(1+x**2), 0, 1, 10**(-6))
print("tolerance = 10^-6: ", test2)

actual = Simpson(1/(1+x**2), 0, 1, 0.001)
print("actual = ", Simpson(1/(1+x**2), 0, 1, 0.001))
print("error = ", abs(actual - test2))
'''

Q1Ans = Simpson(4**(x-1), 0, 2, 0.01)
print("Q1 Answer = ", Q1Ans)

Q2Ans = SimpsonAdaptive(4**(x-1), 0, 2, 0.0001)
print("Q2 Answer = 17 nodes")

# Task B: K-th order derivative


def Derivative(x, xn, yn, k):
    xi = np.where(xn == x)[0][0]

    def Solver(x, xn, yn, k, forwards):
        h = xn[1] - xn[0]
        xi = np.where(xn == x)[0][0]
        if k == 1:
            return (yn[xi + forwards] - yn[xi + forwards - 1]) / h
        else:
            a = Solver(x + h * forwards, xn, yn, k - 1, forwards=forwards)
            b = Solver(x + h * forwards - h, xn, yn, k - 1, forwards=forwards)
            return (a - b) / h

    if xi < k:
        return Solver(x, xn, yn, k, forwards=True)
    else:  # default = backwards: more realistic for this situation as the future is not known
        return Solver(x, xn, yn, k, forwards=False)


def DerivativeDebug(x, xn, yn, k):
    xi = np.where(xn == x)[0][0]

    def Solverd(x, xn, yn, k, forwards):
        h = xn[1] - xn[0]
        xi = np.where(xn == x)[0][0]
        if k == 1:
            return (yn[xi + forwards] - yn[xi + forwards - 1]) / h
        else:
            a = Solverd(x + h * forwards, xn, yn, k - 1, forwards=forwards)
            b = Solverd(x + h * forwards - h, xn, yn, k - 1, forwards=forwards)
            return a

    if xi < k:
        return Solverd(x, xn, yn, k, forwards=True)
    else:  # default = backwards: more realistic for this situation as the future is not known
        return Solverd(x, xn, yn, k, forwards=False)


# test3 = Derivative(5, np.array([0,1,2,3,4,5]), np.array([0,1,2,3,4,5]), 1, forwards=False)
# print(test3)


# Task C: Smoothing derivatives with polynomial interpolation. Launch of a rocket.
# This is a weird rocket, no downrange velocity, super high altitude, space parachute?


r = open("Rocket.txt", "r")
r_read = r.readlines()
r_alt = []
for item in r_read:
    r_alt.append(float(item.rstrip()))
r.close()

times = np.arange(0, 1400, 100)
#print(times)
#print(r_alt)
# plt.plot(times,r_alt)
# plt.show()

r_vel = np.zeros(len(times))
r_accel = np.zeros(len(times))

for i in range(len(times)):
    r_vel[i] = Derivative(times[i], times, r_alt, 1)
    r_accel[i] = Derivative(times[i], times, r_alt, 2)

#print(r_vel)
#print(r_accel)


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


times_interp = np.linspace(0,1300,140)
r_alt_interp = np.zeros(140)
for i in range(len(times_interp)):
    r_alt_interp[i] = Splines(times_interp[i], times, r_alt, 0, 0)

#plt.plot(times_interp, r_alt_interp)
#plt.show()

r_vel_interp = np.zeros(len(times_interp))
r_accel_interp = np.zeros(len(times_interp))

for i in range(len(times_interp)):
    r_vel_interp[i] = Derivative(times_interp[i], times_interp, r_alt_interp, 1)
    #r_accel_interp[i] = Derivative(times_interp[i], times_interp, r_alt_interp, 2)

#print(times_interp)
#print(r_alt_interp)
#print(r_vel_interp)
#print(r_accel_interp)

# rounded times:
times_rounded = []
for i in range(len(times_interp)):
    times_rounded.append(int(times_interp[i]))
time_index = times_rounded.index(374)

Q3Ans = r_vel_interp[time_index]
print("Q3 Answer = ", Q3Ans)
