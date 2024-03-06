# In this section I am importing all the libraries I will need
import math
from math import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import keyboard  # just so erronous plots can be stopped more QUICKLy


# In this section I am setting the domain of solution and the discretised grid

xb = 1
yb = 1 # xb and yb give the coordinates of the upper right corner of the bounding rectangle, with the origin being the other

hx = 0.02 # step sizes and end time
hy = 0.02
k = 0.02
t_end = 15

nx = int(xb / hx) + 1
ny = int(yb / hy) + 1
nt = int(t_end / k) + 1

c = 0.5

# In this section I am defining arrays I would need (if needed)

U = np.ndarray((nt, nx, ny))

x = np.arange(0,xb+hx,hx)  #meshgrid for surface
y = np.arange(0,yb+hy,hy)

Xg, Yg = np.meshgrid(x, y)

# In this section I am setting the boundary conditions/initial values

# first version: square box bounding values

U[0,:,:] = 0   # no initial displacement everywhere - boundary value in time
U[1,:,:] = 0   # need another boundary condition for accel (should improve later!)
U[:,:,0], U[:,:,-1] = 0, 0 # no displacement at the boundaries always
U[:,0,:], U[:,-1,:] = 0, 0

# larger central perterbation version
#for x in [i for i in range(1, nx-1) if sqrt((i-nx/2)**2) < 1]:
#    for y in [i for i in range(1, ny-1) if sqrt((x-nx/2)**2+(i-ny/2)**2) < 1]:
#        U[:,x,y] = [sin(0.1*i) for i in range(nt)]
#

# point central perterbation version
U[:,int(nx/2),int(ny/2)] = [2*sin(0.09*i) for i in range(nt)]   # oscillating point at the centre


# In this section I am implementing the numerical methode

for t in range(2, nt):
        for x in range(1, nx-1):#[i for i in range(1, nx - 1) if i != int(nx/2)]:
            for y in range(1, ny-1):#[i for i in range(1, ny - 1) if i != int(ny/2)]:
                if x != int(nx/2) or y != int(ny/2):# central area version: sqrt((x-nx/2)**2+(y-ny/2)**2) >= 1:
                    uxx = (1/hx**2) * (U[t-1,x+1,y] - 2*U[t-1,x,y] + U[t-1,x-1,y])
                    uyy = (1/hy**2) * (U[t-1,x,y+1] - 2*U[t-1,x,y] + U[t-1,x,y-1])
                    U[t,x,y] = 2*U[t-1,x,y] - U[t-2,x,y] + k**2*c**2*uxx + k**2*c**2*uyy
                #U[t,x,y] = (t*c/(hx*hy)) * (U[t-1,x-1,y] + U[t-1,x+1,y] + U[t-1,x,y-1] + U[t-1,x,y+1] - 4*U[t-1,x,y]) + 2*U[t-1,x,y] - U[t-2,x,y]


# In this section I am showing the results

norm = matplotlib.colors.Normalize(vmin=-1, vmax=1)

'''
for i in range(int(t_end / k) + 1):
    plt.imshow(U[i], interpolation='bilinear', norm=norm)
    plt.colorbar()
    plt.pause(0.000001)
    plt.clf()
    if keyboard.is_pressed('e'):  # if key 'q' is pressed 
            print('Abort!')
            break  # finishing the loop
'''

#testing surface plot



#'''
ax = plt.axes(projection='3d')   # comment this away for contour plots
for i in range(int(t_end / k) + 1):   
    ax.plot_surface(Xg,Yg,U[i,:,:])
    ax.plot_surface(Xg,Yg,U[i,:,:])
    ax.set(xlim=(0, xb), ylim=(0, yb), zlim=(-3,3))
    plt.pause(0.000001)
    plt.cla()
    if keyboard.is_pressed('e'):  # if key 'q' is pressed 
            print('Abort!')
            break  # finishing the loop
#'''
    



# In this section I am celebrating
print('CW done: I deserve a good mark')

