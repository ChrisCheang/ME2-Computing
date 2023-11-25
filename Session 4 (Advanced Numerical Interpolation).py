import math
import numpy as np
import matplotlib.pyplot as plt
import sympy as smp
from sympy import *

#from Session_3_Numerical_Integration import lagrange

x = symbols("x", real=True)

# Task A: Simpson Integration adn Adaptive Simpson Integration


def Simpson(fx, a, b, dx):
    f = smp.sympify(fx)
    xn = np.arange(a, b+dx, dx)
    I = 0
    for i in range(len(xn)):
        yi = f.subs(x, xn[i]).evalf()
        if i == 0 or i == len(xn)-1:
            I += yi
        elif (i % 2) == 0:  #even
            I += 2 * yi
        else:
            I += 4 * yi
    h = (b-a)/(len(xn)-1)
    return I * h/3


def SimpsonAdaptive(fx, a, b, error):
    dx = (b-a)/2
    eps = Simpson(fx, a, b, (b-a)/4) - Simpson(fx, a, b, dx) #first error estimation
    while eps > error:
        Ihalf = Simpson(fx, a, b, dx/2)
        Ih = Simpson(fx, a, b, dx)
        eps = (Ihalf-Ih)/15
        dx = dx/2
    return Simpson(fx, a, b, dx)


'''
test2 = SimpsonAdaptive(1/(1+x**2), 0, 1, 10**(-6))
print("tolerance = 10^-6: ", test2)

actual = Simpson(1/(1+x**2), 0, 1, 0.001)
print("actual = ", Simpson(1/(1+x**2), 0, 1, 0.001))
print("error = ", abs(actual - test2))
'''

# Task B: K-th order derivative


def Derivative(x, xn, yn, k):
    xi = np.where(xn == x)[0][0]
    
    def Solver(x, xn, yn, k, forwards):
        h = xn[1] - xn[0]
        xi = np.where(xn == x)[0][0]
        if k == 1:
            return (yn[xi + forwards] - yn[xi + forwards - 1]) / h
        else:
            a = Solver(x + h*forwards, xn, yn, k - 1, forwards=forwards)
            b = Solver(x + h*forwards - h, xn, yn, k - 1, forwards=forwards)
            return (a - b) / h
            
    if xi < k:
        return Solver(x, xn, yn, k, forwards=True)
    else:  # default = backwards: more realistic for this situation as the future is not known
        return Solver(x, xn, yn, k, forwards=False)


#test3 = Derivative(5, np.array([0,1,2,3,4,5]), np.array([0,1,2,3,4,5]), 1, forwards=False)
#print(test3)


# Task C: Smoothing derivatives with polynomial interpolation. Launch of a rocket.
# This is a weird rocket, no downrange velocity, super high altitude, space parachute?


r = open("Rocket.txt", "r")
r_read = r.readlines()
r_alt = []
for item in r_read:
    r_alt.append(float(item.rstrip()))
r.close()

times = np.arange(0, 1400, 100)
print(times)
print(r_alt)
#plt.plot(times,r_alt)
#plt.show()

r_vel = np.zeros(len(times))
r_accel = np.zeros(len(times))

for i in range(len(times)):
    r_vel[i] = Derivative(times[i], times, r_alt, 1)
    r_accel[i] = Derivative(times[i], times, r_alt, 2)




print(r_vel)
print(r_accel)
