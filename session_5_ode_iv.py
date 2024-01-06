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


def FwEulerN(Fn, t0, Y0, h, t_end, time_nodes_output=False):
    t_array = np.array([t0])
    Y_array = np.expand_dims(Y0, axis=1)
    tn, Yn = t0, Y_array
    while tn < t_end - h:
        for var in range(len(Fn)):
            f = smp.sympify(Fn[var])
            f_eval = f.subs(t, t_array[-1])
            for i in range(len(Fn)):
                f_eval = f_eval.subs(y[i], Y_array[i][-1])
            Yn[var][0] += h * f_eval.evalf()
            #return Yn[var][0]
        tn += h
        t_array = np.append(t_array, [tn])
        Y_array = np.append(Y_array, Yn, axis=1)
    if not time_nodes_output:
        return Y_array
    else:
        return t_array


y = symbols('y1:%d' % 4)


a = 0.001
b = 0.05

SIR = [-a*y[0]*y[1], a*y[0]*y[1] - b*y[1], b*y[1]]
Y0 = np.array([500, 10, 0])

model = FwEulerN(SIR, 0, Y0, 0.2, 200)
times = FwEulerN(SIR, 0, Y0, 0.2, 200, time_nodes_output=True)

print(model)



#'''
plt.plot(times,model[0],color='yellow')
plt.plot(times,model[1],color='red')
plt.plot(times,model[2],color='green')
plt.show()
#'''



# Y = np.array([[1,1,1,1,1]])
# Y = np.transpose(Y)
# Y = np.append(Y, np.array([[2],[2],[2],[2],[2]]), axis=1)

# print(Y)

#print(np.empty((3, 1)))
