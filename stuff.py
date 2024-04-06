import numpy as np
from math import *
import time
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sympy as smp

# Tutorial Sheet 11 1d,e,f)
'''

# important note: all the numbers have to be floats, or else each step will be rounded.
A = np.array([[-4.,2.,0.],[1.,-4.,1.],[0.,2.,-4.]])
b = np.array([-2.,-44.,238.])

u = np.array([0.,0.,0.])
uerr = np.zeros(3)
uanal = np.array([-2,-5,-62])

tolerance = False
iterations = 0
relaxation = 1.1

while tolerance == False:
    ulast = [u[i] for i in range(3)]
    u[0] = (b[0]-A[0,1]*u[1]-A[0,2]*u[2])/A[0,0]  #0.5*(1+u[1])#
    u[0] = relaxation*u[0] + (1-relaxation)*ulast[0]    # additional step for relaxation (repeat for the other ones)
    u[1] = (b[1]-A[1,0]*u[0]-A[1,2]*u[2])/A[1,1]  #0.25*(44+u[0]+u[2])#
    u[1] = relaxation*u[1] + (1-relaxation)*ulast[1]
    u[2] = (b[2]-A[2,0]*u[0]-A[2,1]*u[1])/A[2,2]  #-0.5*(119-u[1])#
    u[2] = relaxation*u[2] + (1-relaxation)*ulast[2]
    uerr = [100*abs((u[i] - ulast[i])/ulast[i]) for i in range(3)]  #100*abs((u[i] - ulast[i])/ulast[i]) 
    error_count = 0
    for i in range(3):
        if uerr[i] < 1:
            error_count += 1
    
    if error_count == 3:
        tolerance = True

    iterations += 1

    print(u, uerr)

print("no. of iterations = ", iterations)

'''

# Tutorial Sheet 12 2b)

'''


# explicit method

h = 0.25
k = 0.03125
r = k/h**2
t_end = 2

nx = int(1/h+1)
xs = [i*h for i in range(nx)]
nt = int(t_end/k + 1)

u = np.ndarray((nx,nt))
u[:,0] = np.ones(nx)
u[-1,0] = 0


for j in range(nt-1):
    u[0,j+1] = r*(u[1,j]+u[1,j])+(1-2*r)*(u[0,j])
    for i in range(1,nx-1):
        u[i,j+1] = r*(u[i+1,j]+u[i-1,j])+(1-2*r)*(u[i,j])
    u[nx-1,j+1] = 0



# analytical solution

def T_analytical(x,t,n=2,a=1,l=1):
    sum = 0
    for i in range(1,n):
        sum += ((-1)**(i-1)/(2*i-1)) * e**(-a*(2*i-1)**2*pi**2*t/(4*l**2)) * cos((2*i-1)*(pi*x/(2*l)))
    return 4*sum/pi

#print("12f) = ", T_analytical(0.5,1))

# Crank nicolson method

h = 0.25  # method only works with 0.25
k = 0.03125
rv = k/h**2

t_end = 2

nx = int(1/h+1)
xs = [i*h for i in range(nx)]
xsa = np.arange(0,1.05,0.05)
nt = int(t_end/k + 1)

u_cn = np.ndarray((nx,nt))
u_cn[:,0] = np.ones(nx)
u_cn[-1,0] = 0

smp.init_printing()
r = smp.symbols('r')
A = smp.Matrix([[2+2*r, -2*r,0,0], 
                [-r,2+2*r,-r,0],
                [0,-r,2+2*r,-r],
                [0,0,-r,2+2*r]])
Ainv = A.inv()

for j in range(nt-1):
    b = smp.Matrix([[(2-2*r)*u_cn[0,j] + 2*r*u_cn[1,j]],
                    [r*u_cn[0,j] + (2-2*r)*u_cn[1,j] + r*u_cn[2,j]],
                    [r*u_cn[1,j] + (2-2*r)*u_cn[2,j] + r*u_cn[3,j]],
                    [r*u_cn[2,j] + (2-2*r)*u_cn[3,j]]])
    un = Ainv*b
    u_cn[:,j+1] = [un.evalf(subs={r:rv})[i] for i in range(4)] + [0]
    u_cn[nx-1,j+1] = 0

    plt.title(f"time = {j*k}")
    plt.axis([0, 1, 0, 3])
    plt.plot(xs,u[:,j+1])   # explicit
    plt.plot(xs,u_cn[:,j+1])    # Crank-Nicolson
    plt.plot(xsa,[T_analytical(xsa[i],j*k) for i in range(21)])   # analytical
    plt.legend(['Explicit', 'Crank-Nicolson', 'Analytical'])
    plt.pause(0.05)
    plt.clf()

'''

# Tutorial Sheet 13 1a) - explicit 1D wave equation with dirichlet bc
'''

h = 0.2
k = h
r = k/h**2
t_end = 2

nx = int(1/h+1)
nt = int(t_end/k + 1)
xs = [i*h for i in range(nx)]

u = np.ndarray((nx,nt))

# Boundary conditions and initial condition in displacement and velocity
f = [sin(pi*h*i) for i in range(nx)]
u[:,0] = f
g = np.zeros(nx)
u[0,:], u[-1,:] = 0, 0

# j=1 calculations using initial condition in velocity 
u[1:nx-1,1] = [0.5*(u[i-1,0]+u[i+1,0]) + k*g[i] for i in range(1,nx-1)]

# rest of the domain
for j in range(2,nt):
    u[1:nx-1,j] = [u[i-1,j-1]+u[i+1,j-1]-u[i,j-2] for i in range(1,nx-1)]

    plt.title(f"time = {j*k}")
    plt.axis([0, 1, -1, 1])
    plt.plot(xs,u[:,j])   # explicit
    plt.pause(0.005)
    plt.clf()

'''

# Tutorial Sheet 13 1a) - explicit 1D wave equation with neumann bc

def u_analytical(x,t):
    return cos(pi*x)*cos(pi*t)

h = 0.1
k = h
r = k/h**2
t_end = 4

nx = int(1/h+1) + 2
nt = int(t_end/k + 1)
xs = [(i-1)*h for i in range(nx)]
xsa = np.arange(0,1.05,0.05)

u = np.ndarray((nx,nt))

# Boundary conditions and initial condition in displacement and velocity

f = [cos(pi*h*(i-1)) for i in range(nx)]
u[:,0] = f
g =  np.zeros(nx)

# j=1 calculations using initial condition in velocity 
# beware the difference between code index (0 to nx) and discretisation index (-1 to nx-1)

u[1,1] = f[2] + k*g[1] # x = 0 boundary for velocity i.c.
u[2:nx-2,1] = [0.5*(f[i-1]+f[i+1]) + k*g[i] for i in range(2,nx-2)]
u[-2,1] = f[-3] + k*g[-2] # x = 1 boundary for velocity i.c.
u[0,1] = u[2,1] # left ghost point
u[-1,1] = u[-3,1] # right ghost point

# rest of the domain
for j in range(2,nt):
    u[1,j] = 2*u[2,j-1] - u[1,j-2]  # x = 0 boundary with neumann b.c. considered (u0 = u2)
    u[2:nx-2,j] = [u[i-1,j-1] + u[i+1,j-1] - u[i,j-2] for i in range(2,nx-2)]
    u[-2,j] = 2*u[-3,j-1] - u[-2,j-2]  # x = 1 boundary with neumann b.c. considered (u-3 = u-1)
    u[0,j] = u[2,j] # left ghost point
    u[-1,j] = u[-3,j] # right ghost point

    plt.title(f"time = {j*k}")
    plt.axis([0, 1, -1, 1])
    plt.plot(xs[1:nx-1],u[1:nx-1,j])   # explicit
    plt.plot(xsa,[u_analytical(xsa[i],j*k) for i in range(21)])   # analytical
    plt.legend(['Explicit', 'Analytical'])
    plt.pause(h)
    plt.clf()

# error calculation

u_error = np.zeros((nx,nt))
for j in range(0,nt):
    u_error[1:nx-1,j] = [abs(u_analytical(xs[i],j*k) - u[i,j])/u_analytical(xs[i],j*k) for i in range(1,nx-1)]

