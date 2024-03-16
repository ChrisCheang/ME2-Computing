# In this section I am importing all the libraries I will need
import math
from math import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import keyboard  # just so erronous plots can be stopped more QUICKLy


# Resources
# https://uwaterloo.ca/computational-mathematics/sites/ca.computational-mathematics/files/uploads/files/linqi_shao_report_pdf.pdf#page16
# http://spiff.rit.edu/classes/phys283/lectures/two_d/two_d.html
# https://www.comsol.com/blogs/how-do-chladni-plates-make-it-possible-to-visualize-sound/
# https://thelig.ht/chladni/


class Plate:
    def __init__(self, xb, yb, h, k, t_end, freq, c):
        self.xb = xb
        self.yb = yb
        self.h = h
        self.k = k
        self.freq = freq
        self.t_end = t_end
        self.c = c

        self.nx = int(xb / h) + 1
        self.ny = int(yb / h) + 1
        self.nt = int(t_end / k) + 1

    def mesh(self):
        x = np.arange(0,self.xb+self.h,self.h)  #meshgrid for surface
        y = np.arange(0,self.yb+self.h,self.h)

        Xg, Yg = np.meshgrid(x, y)
        return (Xg, Yg)
    
    def solve_matrix(self):
        # In this section I am defining arrays I would need (if needed)

        U = np.zeros((self.nt, self.nx, self.ny))

        # In this section I am setting the boundary conditions/initial values

        # first version: square box bounding values

        U[0,:,:] = 0   # no initial displacement everywhere - boundary value in time
        U[1,:,:] = 0   # need another boundary condition for accel (should improve later!)
        #U[:,:,0], U[:,:,-1] = 0, 0 # no displacement at the boundaries always
        #U[:,0,:], U[:,-1,:] = 0, 0

        # larger central perterbation version
        #for x in [i for i in range(1, nx-1) if sqrt((i-nx/2)**2) < 1]:
        #    for y in [i for i in range(1, ny-1) if sqrt((x-nx/2)**2+(i-ny/2)**2) < 1]:
        #        U[:,x,y] = [sin(0.1*i) for i in range(nt)]
        #

        # point central perterbation version
        U[:,int(self.nx/2),int(self.ny/2)] = [1*sin(2*pi*self.freq*i) for i in range(self.nt)]# + [0 for i in range(nt-20)]   # oscillating point at the centre


        # In this section I am implementing the numerical method

        for t in range(2, self.nt):
            for x in range(1, self.nx-1):#[i for i in range(1, nx - 1) if i != int(nx/2)]:
                for y in range(1, self.ny-1):#[i for i in range(1, ny - 1) if i != int(ny/2)]:
                    if x != int(self.nx/2) or y != int(self.ny/2):# central area version: sqrt((x-nx/2)**2+(y-ny/2)**2) >= 1:
                        uxx = (1/self.h**2) * (U[t-1,x+1,y] - 2*U[t-1,x,y] + U[t-1,x-1,y])
                        uyy = (1/self.h**2) * (U[t-1,x,y+1] - 2*U[t-1,x,y] + U[t-1,x,y-1])
                        U[t,x,y] = 2*U[t-1,x,y] - U[t-2,x,y] + self.k**2*self.c**2*uxx + self.k**2*self.c**2*uyy
                        #U[t,x,y] = (t*c/(hx*hy)) * (U[t-1,x-1,y] + U[t-1,x+1,y] + U[t-1,x,y-1] + U[t-1,x,y+1] - 4*U[t-1,x,y]) + 2*U[t-1,x,y] - U[t-2,x,y]
            U[t,0,:] = U[t,1,:]
            U[t,-1,:] = U[t,-2,:]
            U[t,:,0] = U[t,:,1]
            U[t,:,-1] = U[t,:,-2]
        
        return U
    
    def visualise(self):
        # In this section I am showing the results

        U = self.solve_matrix()
        Xg = self.mesh()[0]
        Yg = self.mesh()[1]

        norm = matplotlib.colors.Normalize(vmin=-2, vmax=2)

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


        #'''
        fig = plt.figure()
        ax1 = fig.add_subplot(221, projection='3d')
        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')
        ax1.set_zlabel('U')

        ax2 = fig.add_subplot(212)
        ax2.set_xlabel('X')
        ax2.set_ylabel('Y')

        ax3 = fig.add_subplot(223)
        ax3.set_xlabel('X')
        ax3.set_ylabel('Y')

        ax1.set_box_aspect([1, 1, 1])  # set aspect ratio of 3D plot
        ax2.set_aspect('equal')  # set aspect ratio of contour plot
        ax3.set_aspect('equal')  # set aspect ratio of colour plot

        colour = 'coolwarm'

        # Plot color key (snapshot chosen in the middle of time range so we get the right max and mins - find more refined version)
        fig.colorbar(ax2.contour(Xg, Yg, U[int(self.nt/2)], cmap=colour), ax=ax2)

        for i in range(int(self.t_end / self.k) + 1):
            ax1.clear()
            ax2.clear()
            ax3.clear()
            ax1.plot_surface(Xg, Yg, U[i], cmap=colour, vmin=-2, vmax=2)
            ax1.set_zlim([-2, 2])  # set the z-axis limits
            ax2.contour(Xg, Yg, U[i], cmap=colour, vmin=-2, vmax=2)
            ax3.imshow(U[i], interpolation='bilinear', norm=norm, extent=[0, self.xb, 0, self.yb], origin='lower', cmap=colour)
            #ax3.colorbars()

            plt.pause(0.000001)

            if keyboard.is_pressed('e'):  # if key 'e' is pressed 
                plt.close()
                print('Abort!')
                break  # finishing the loop
            
        plt.show()

        #'''

test_plate = Plate(1,1,0.02,0.02,15,0.03/pi,0.5)
test_plate.visualise()





