import math
import numpy as np
import matplotlib.pyplot as plt
import sympy as smp
from sympy import *

y = symbols("y", real=True)
t = symbols("t", real=True)


# Task A: Explicit methods: Forward Euler and RK4


def FwEuler(f, t0, y0, h, t_end, time_nodes_output=False):
    f = smp.sympify(f)
    t_array, y_array = np.array([t0]), np.array([y0])
    tn, yn = t0, y0
    while tn < t_end - h:
        yn += h * f.subs({t: t_array[-1], y: y_array[-1]}).evalf()
        tn += h
        t_array = np.append(t_array, [tn])
        y_array = np.append(y_array, [yn])
    if not time_nodes_output:
        return y_array
    else:
        return t_array


def ODERK4(f, t0, y0, h, t_end, time_nodes_output=False):
    f = smp.sympify(f)
    t_array, y_array = np.array([t0]), np.array([y0])
    tn, yn = t0, y0
    while tn < t_end - h:
        k1 = h * f.subs({t: t_array[-1], y: y_array[-1]}).evalf()
        k2 = h * f.subs({t: t_array[-1] + 0.5 * h, y: y_array[-1] + 0.5 * k1}).evalf()
        k3 = h * f.subs({t: t_array[-1] + 0.5 * h, y: y_array[-1] + 0.5 * k2}).evalf()
        k4 = h * f.subs({t: t_array[-1] + h, y: y_array[-1] + k3})
        yn += (k1 + 2 * k2 + 2 * k3 + k4) / 6
        tn += h
        t_array = np.append(t_array, [tn])
        y_array = np.append(y_array, [yn])
    if not time_nodes_output:
        return y_array
    else:
        return t_array


# Task B: Implicit methods: Backward Euler


def BwEuler(f, t0, y0, h, t_end, time_nodes_output=False):
    f = smp.sympify(f)
    t_array, y_array = np.array([t0]), np.array([y0])
    tn, yn = t0, y0
    yn1 = symbols("yn1", real=True)
    while tn < t_end - h:
        yn = solve(yn1 - y_array[-1] - h * f.subs({t: t_array[-1] + h, y: yn1}).evalf(), yn1)[0]
        tn += h
        t_array = np.append(t_array, [tn])
        y_array = np.append(y_array, [yn])
    if not time_nodes_output:
        return y_array
    else:
        return t_array


test_ode = -2 * y * t - 2 * t ** 3

h = 0.2
xeuler = FwEuler(test_ode, 0, -3, h, 2, time_nodes_output=True)
yeuler = FwEuler(test_ode, 0, -3, h, 2, time_nodes_output=False)

xrk4 = ODERK4(test_ode, 0, -3, h, 2, time_nodes_output=True)
yrk4 = ODERK4(test_ode, 0, -3, h, 2, time_nodes_output=False)

xbeuler = BwEuler(test_ode, 0, -3, h, 2, time_nodes_output=True)
ybeuler = BwEuler(test_ode, 0, -3, h, 2, time_nodes_output=False)

xactual = np.arange(0, 2.01, 0.01)
yactual = []
for i in range(len(xactual)):
    yactual.append(1 - xactual[i] ** 2 - 4 * math.e ** (-xactual[i] ** 2))

'''
plt.plot(xeuler, yeuler, color='b')
plt.plot(xrk4, yrk4, color='r')
plt.plot(xbeuler, ybeuler, color='g')
plt.plot(xactual, yactual, color='black')
plt.show()
'''

# Task C: System of ODEs, with explicit method

# Task C1: Covid-19 model (SIR model)


def FwEulerN(Fn, t0, Y0, h, t_end, time_nodes_output=False):
    F = []
    for var in range(len(Fn)):
        F.append(smp.sympify(Fn[var]))
    t_array = np.array([t0])
    Y_array = np.expand_dims(Y0, axis=1)
    tn = t0
    Yn = np.zeros((len(Fn), 1))
    while tn < t_end - h:
        for var in range(len(Fn)):
            f_eval = F[var].subs(t, t_array[-1])
            for i in range(len(Fn)):
                f_eval = f_eval.subs(y[i], Y_array[i][-1])
            Yn[var][0] = Y_array[var][-1] + h * f_eval.evalf()
        tn += h
        t_array = np.append(t_array, [tn])
        Y_array = np.append(Y_array, Yn, axis=1)
    if not time_nodes_output:
        return Y_array
    else:
        return t_array


def FwEulerN_SIR(t0, Y0, h, t_end, time_nodes_output=False):   # runs a bit faster
    t_array = np.array([t0])
    Y_array = np.expand_dims(Y0, axis=1)
    tn = t0

    while tn < t_end - h:
        Yn = np.zeros((3,1))

        def SIR(var):
            S = Y_array[0][-1]
            I = Y_array[1][-1]
            R = Y_array[2][-1]
            if var == 0:  # S
                return - a * S * I
            elif var == 1:  # I
                return a * S * I - b * I
            else:  # R
                return b * I
        for var in range(3):
            Yn[var][0] = Y_array[var][-1] + h * SIR(var)

        tn += h
        t_array = np.append(t_array, [tn])
        Y_array = np.append(Y_array, Yn, axis=1)
    if not time_nodes_output:
        return Y_array
    else:
        return t_array

'''
y = symbols('y1:%d' % 4)

a = 0.001
b = 0.05
SIR = [-a*y[0]*y[1], a*y[0]*y[1] - b*y[1], b*y[1]]
Y0 = np.array([500, 10, 0])

model_SIR = FwEulerN(SIR, 0, Y0, 0.05, 200)
times_SIR = FwEulerN(SIR, 0, Y0, 0.05, 200, time_nodes_output=True)

# faster version with FwEulerN_SIR, not sure why its faster
model = FwEulerN_SIR(0, Y0, 0.05, 200)
times = FwEulerN_SIR(0, Y0, 0.05, 200, time_nodes_output=True)



print(model_SIR)
plt.plot(times_SIR,model_SIR[0],color='yellow')
plt.plot(times_SIR,model_SIR[1],color='red')
plt.plot(times_SIR,model_SIR[2],color='green')
plt.show()
'''


# Task C2: Financial model of the house market in London (Lotka-Volterra)

'''
y = symbols('y1:%d' % 3)

House = [0.3*y[0]*y[1] - 0.8*y[0], 1.1*y[1] - y[0]*y[1]]
Y0_H = np.array([0.8, 7])

model_house = FwEulerN(House, 0, Y0_H, 0.005, 1)
times_house = FwEulerN(House, 0, Y0_H, 0.005, 1, time_nodes_output=True)


plt.plot(times_house,model_house[0],color='yellow', label='Price (100K pounds)')
plt.plot(times_house,model_house[1],color='red', label='No. sold (K houses)')
plt.xlabel("Time / 40 months")
plt.legend()
plt.show()

plt.plot(model_house[0],model_house[1], label='price vs no.sold')
plt.xlabel("Price (100K pounds)")
plt.ylabel("No. sold (K houses)")
plt.legend()
plt.show()
'''

# Task D: Higher order ODEs

# Task D1: Damped non-linear motion of a pendulum


#'''
y = symbols('y0:%d' % 2)  # theta, Dtheta  (note: D2theta is not a variable!)

g = 9.81
m = 0.5
c = 0.05
L = 1
Pendulum = [y[1], - c*y[1]/m - g*sin(y[0])/L]
Y0_P = np.array([pi/4, 0])


model_P = FwEulerN(Pendulum, 0, Y0_P, 0.005, 15)
times_P = FwEulerN(Pendulum, 0, Y0_P, 0.005, 15, time_nodes_output=True)

print(model_P)

#plt.plot(times_P, model_P[0], color='green', label='Displacement')
#plt.xlabel("time (seconds)")
#plt.ylabel("theta (rad)")
#plt.legend()
#plt.show()

#plt.plot(times_P, model_P[1], color='green', label='Angular Velocity')
#plt.xlabel("time (seconds)")
#plt.ylabel("Dtheta (rad/s)")
#plt.legend()
#plt.show()


plt.plot(times_P, model_P[0], color='green', label='Displacement')
plt.plot(times_P, model_P[1], color='yellow', label='Angular Velocity')
plt.xlabel("time (seconds)")
plt.legend()
plt.show()
#'''
