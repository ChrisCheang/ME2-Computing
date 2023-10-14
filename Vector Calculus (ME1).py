import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from sympy.vector import CoordSys3D
from matplotlib.colors import LightSource
import math
import sympy as smp
from sympy import *


t = symbols("t", real=True)
s = symbols("s", real=True)
x = symbols("x", real=True)
y = symbols("y", real=True)
z = symbols("z", real=True)


class Point:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def cr(self, d):
        if d == "x":
            return self.x
        elif d == "y":
            return self.y
        elif d == "z":
            return self.z
        else:
            return "no"

    def print(self):
        print(self.x,self.y,self.z)


class Vector:
    def __init__(self, x, y, z):
        self.i = x
        self.j = y
        self.k = z

    def print(self):
        print(f"({self.i})i + ({self.j})j + ({self.k})k")

    def draw(self, position=(0,0,0), start=0, end=1, num=50, color="black"):
        line = Curve(position[0]+ t * self.i,position[1] + t * self.j,position[2]+ t * self.k)
        line.draw(start=start, end=end, num=num, color=color)

    def cross(self, b):
        n1 = self.j*b.k - self.k*b.j
        n2 = self.k*b.i - self.i*b.k
        n3 = self.i*b.j - self.j*b.i
        return Vector(n1,n2,n3)



class Curve:
    def __init__(self, x, y, z):
        self.x = smp.sympify(x)
        self.y = smp.sympify(y)
        self.z = smp.sympify(z)

    def print(self):
        print(f"x(t) = {self.x}")
        print(f"y(t) = {self.y}")
        print(f"z(t) = {self.z}")

    def evaluate(self, a):
        xt = self.x.subs(t, a).evalf(3, chop=True)
        yt = self.y.subs(t, a).evalf(3, chop=True)
        zt = self.z.subs(t, a).evalf(3, chop=True)
        return Point(xt, yt, zt)

    def tangent_vec(self, a):
        dx = diff(self.x, t).subs(t, a).evalf(3, chop=True)
        dy = diff(self.y, t).subs(t, a).evalf(3, chop=True)
        dz = diff(self.z, t).subs(t, a).evalf(3, chop=True)
        mag = math.sqrt(dx**2 + dy**2 + dz**2)
        return Curve(self.evaluate(t, a).cr("x") + t * dx / mag, self.evaluate(t, a).cr("y") + t * dy / mag, self.evaluate(t, a).cr("z") + t * dz / mag)

    def draw(self, start=0, end=1, num=50, color='black'):
        tlist = np.linspace(start, end, num=num)
        x = np.zeros(50)
        y = np.zeros(50)
        z = np.zeros(50)
        for i in range(num):
            x[i] = self.evaluate(tlist[i]).cr("x")
            y[i] = self.evaluate(tlist[i]).cr("y")
            z[i] = self.evaluate(tlist[i]).cr("z")
        ax.plot(x, y, z, label=f"({self.x})i + ({self.y})j + ({self.z})k", color=color)
        ax.legend()

    def scalar_line_integral(self, scalar_function, start, end, expression=True):
        dx = diff(self.x, t)
        dy = diff(self.y, t)
        dz = diff(self.z, t)
        phi = smp.sympify(scalar_function)
        exp = phi * sqrt(dx**2 + dy**2 + dz**2)
        if expression:
            print(f"line integral = {integrate(exp, (t, start, end))}")
        else:
            print(f"line integral = {integrate(exp, (t, start, end)).evalf()}")


class Scalar_field:
    def __init__(self, expression):
        self.f = smp.sympify(expression)

    def print(self):
        print(f"phi = {self.f}")

    def print_grad(self):
        print(f"grad_phi_x = {diff(self.f, x)}")
        print(f"grad_phi_y = {diff(self.f, y)}")
        print(f"grad_phi_z = {diff(self.f, z)}")

    def grad(self):
        xv = diff(self.f, x)
        yv = diff(self.f, y)
        zv = diff(self.f, z)
        return Vector_field(xv,yv,zv)

    def laplacian(self):
        xl = diff(self.f, x, x)
        yl = diff(self.f, y, y)
        zl = diff(self.f, z, z)
        return Vector_field(xl,yl,zl)

    def surface(self, domain=(-5,5), range=(-5,5), size=0.25):
        if diff(self.f, z) != 0:
            return print("only 2d scalar functions allowed!")
        xs = np.arange(domain[0], domain[1], size)
        ys = np.arange(range[0], range[1], size)
        xs, ys = np.meshgrid(xs, ys)
        zs = lambdify([x, y], self.f)(xs,ys)
        ax.plot_surface(xs, ys, zs, cmap=cm.coolwarm, linewidth=0, antialiased=False)


class Vector_field:
    def __init__(self, x, y, z):
        self.x = smp.sympify(x)
        self.y = smp.sympify(y)
        self.z = smp.sympify(z)

    def print(self):
        print(f"field = ({self.x})i + ({self.y})j + ({self.z})k")

    def vis(self, cube_size=1, resolution=0.5, length=0.3, color='blue'):
        a, b, c = np.meshgrid(np.arange(- cube_size,  cube_size + resolution, resolution),
                              np.arange(- cube_size,  cube_size + resolution, resolution),
                              np.arange(- cube_size,  cube_size + resolution, resolution))
        u = lambdify([x, y, z], self.x)(a,b,c)
        v = lambdify([x, y, z], self.y)(a,b,c)
        w = lambdify([x, y, z], self.z)(a,b,c)
        ax.quiver(a, b, c, u, v, w, length=length, normalize=True, color=color)

    def div(self):
        pdx = diff(self.x, x)
        pdy = diff(self.y, y)
        pdz = diff(self.z, z)
        return Scalar_field(pdx+pdy+pdz)

    def curl(self):
        n1 = diff(self.z, y) - diff(self.y, z)
        n2 = diff(self.x, z) - diff(self.z, x)
        n3 = diff(self.y, x) - diff(self.x, y)
        return Vector_field(n1,n2,n3)


class Surface:
    def __init__(self, x, y, z):
        self.x = smp.sympify(x)
        self.y = smp.sympify(y)
        self.z = smp.sympify(z)

    def print(self):
        print(f"field = ({self.x})i + ({self.y})j + ({self.z})k")

    def evaluate(self, a, b):
        xt = lambdify([s, t], self.x)(a, b)
        yt = lambdify([s, t], self.y)(a, b)
        zt = lambdify([s, t], self.z)(a, b)
        return Point(xt, yt, zt)

    def draw(self, s_range=(-1,1), t_range=(-1,1), size=0.05, color="black"):
        sv = np.arange(s_range[0], s_range[1]+size, size)
        tv = np.arange(t_range[0], t_range[1]+size, size)
        sv, tv = np.meshgrid(sv, tv)
        xs = lambdify([s, t], self.x)(sv, tv)
        ys = lambdify([s, t], self.y)(sv, tv)
        zs = lambdify([s, t], self.z)(sv, tv)
        ax.plot_surface(xs, ys, zs, cmap=cm.coolwarm, linewidth=0, antialiased=False)

    def normal(self, a, b):
        u1 = lambdify([s, t], diff(self.x, s))(a, b)
        u2 = lambdify([s, t], diff(self.y, s))(a, b)
        u3 = lambdify([s, t], diff(self.z, s))(a, b)
        v1= lambdify([s, t], diff(self.x, t))(a, b)
        v2 = lambdify([s, t], diff(self.y, t))(a, b)
        v3 = lambdify([s, t], diff(self.z, t))(a, b)
        n1 = u2*v3 - u3*v2
        n2 = u3*v1 - u1*v3
        n3 = u1*v2 - u2*v1
        return Curve(self.evaluate(a,b).cr("x") + t * n1, self.evaluate(a,b).cr("y") + t * n2, self.evaluate(a,b).cr("z") + t * n3)





ax = plt.figure(figsize=(12,9)).add_subplot(projection='3d')

#ax.set_xlim3d(-1, 1)
#ax.set_ylim3d(-1, 1)
#ax.set_zlim3d(-1, 1)

b = Curve(t,t**2,t**3)
c = Scalar_field(2*y*sin(x))
d = Vector_field(y,z*exp(x),5*x)
e = Surface(1+cos(s),2*t,sin(s))
f = Vector(3,4,3)
g = Vector(3,7,5)



d.vis()
d.curl().vis(color='red')
d.print()
d.curl().print()



plt.show()
