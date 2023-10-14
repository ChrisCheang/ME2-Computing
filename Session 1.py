import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import plotly.graph_objects as go

# Task 1
a = np.arange(-5,-1.5,0.5)
b = np.arange(-1.95,3,0.05)
c = np.arange(3,5.5,0.5)
x = np.concatenate((a,b,c))

f, g = [], []
for i in range(x.size):
    f.append(math.sin(x[i]))
    g.append(math.sin((x[i])**2+np.pi))

#fig, ax = plt.subplots()

#ax.scatter(x, y_1, color='red', marker='d')
#ax.scatter(x, y_2, color='purple', marker='o')

#plt.show()

#Task 2

dx = 0.1#*np.pi
x = np.arange(-2*np.pi,2*np.pi+dx,dx)
y = np.arange(-np.pi,2*np.pi+dx,dx)
Nx = len(x)
Ny = len(y)
Xg, Yg = np.meshgrid(x, y)

F = np.sin(Xg)*np.cos(Yg)   # element-wise operation using meshgrid
G = np.cos(Xg)*np.sin(Yg)
S = F + G
P = F * G

print(P[2,10])

#Task 3

# Comment away the others to show a particular graph
#ax = plt.axes(projection='3d')   # comment this away for contour plots
#ax.plot_surface(Xg,Yg,S)
#ax.plot_surface(Xg,Yg,P)
#plt.contour(Xg,Yg,S)
#plt.contour(Xg,Yg,P)
#plt.show()



R = F * np.e**(-0.5*0)
#ax.plot_surface(Xg,Yg,R)

R = F * np.e**(-0.5*5)
#ax.plot_surface(Xg,Yg,R)

#plt.show()

t = np.arange(0,10.05,0.05)
R_spec = np.sin(np.pi)*np.cos(-np.pi/2)*np.e**(-0.5*t)
#plt.scatter(t,R_spec)
#plt.show()

#Task 4

dl = 0.5
l = np.arange(-5,5+dl,dl)
Xg, Yg = np.meshgrid(l, l)
f_1i, f_1j = Xg, Yg
f_2i, f_2j = Yg, -Xg

#plt.quiver(Xg,Yg,f_1i,f_1j)
#plt.streamplot(Xg,Yg,f_1i,f_1j)
#plt.quiver(Xg,Yg,f_2i,f_2j)
#plt.streamplot(Xg,Yg,f_2i,f_2j)
#plt.show()

#f_1 is irrotational with positive divergence everywhere, f_2 is rotational with no divergence everywhere

#Task 5

c = open("Computing.txt", "r")
m = open("Maths.txt", "r")
c_read = c.readlines()
m_read = m.readlines()
c_marks = []
m_marks = []
for item in c_read:
    c_marks.append(int(float(item.rstrip())))
for item in m_read:
    m_marks.append(int(float(item.rstrip())))

c_values, c_counts = np.unique(c_marks, return_counts=True)
m_values, m_counts = np.unique(m_marks, return_counts=True)

#plt.bar(c_values, c_counts)
#plt.bar(m_values, m_counts)

#plt.scatter(c_marks,m_marks)
#plt.xlabel("Computing Marks")
#plt.ylabel("Maths Marks")

#plt.show()

d = 0.1 #*np.pi
x = np.arange(-2,2+d,d)
y = np.arange(-3,3+d,d)
z = np.arange(-np.pi,2*np.pi+d,d)

Xg, Yg, Zg = np.meshgrid(x,y,z)
F = Xg**2 + Yg**2 + Zg**2 - 5*np.sin(Zg)**2

fig = go.Figure(data=go.Isosurface(
    x=Xg.flatten(),
    y=Yg.flatten(),
    z=Zg.flatten(),
    value=F.flatten(),
    isomin=1,
    isomax=10,
    surface_count=5,
    opacity=0.6,
    caps=dict(x_show=False, y_show=False)
    ))
fig.show()


