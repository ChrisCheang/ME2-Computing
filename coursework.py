
# In this section I am importing all the libraries I will need

import math
from math import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

import keyboard  # used for 


class Plate:
    def __init__(self, L, h, k,amp,m,n,p,pv,rplt,q):
        self.L = L
        self.h = h
        self.k = k
        
        self.freq = p*(1/(L)) * sqrt((m**2 + n**2))
        
        self.t_end = q/self.freq
        
        self.amp = amp
        self.m = m
        self.n = n
        
        self.p = p
        
        self.pv = pv
        self.rplt = rplt
        
        self.nx = int(L / h) + 1
        self.ny = self.nx
        self.nt = int(self.t_end / k) + 1


    def mesh(self):
        
        #In this section I am setting the discretized mesh grid
        
        x = np.arange(-self.L/2,self.L/2+self.h,self.h) 
        y = np.arange(-self.L/2,self.L/2+self.h,self.h)

        Xg, Yg = np.meshgrid(x, y)
        return (Xg, Yg)
    
    
    def compare(self):
        
        #This section is the analytical solution for the resonance node pattern, for visual comparison
        
        Xg, Yg = self.mesh()
        Ua = np.cos(2*self.m*np.pi*Xg/self.L)*np.cos(2*self.n*np.pi*Yg/self.L) + np.cos(2*self.n*np.pi*Xg/self.L)*np.cos(2*self.m*np.pi*Yg/self.L)
        return Ua
    
    
    def solve_matrix(self):
        
        # In this section I am defining arrays I would need (if needed)

        U = np.zeros((self.nt, self.nx, self.ny))

        # In this section I am setting the boundary conditions/initial values

       
        U[0,:,:] = 0   
    
        U[1,:,:] = 0   
 
    
        U[:,int(self.nx/2),int(self.ny/2)] = [self.amp*sin((k)*self.freq*i) for i in range(self.nt)]
        
        
        # In this section I am implementing the numerical method
        
        for t in range(1, self.nt-1):
            U[t,0,:] = U[t,1,:]
            U[t,-1,:] = U[t,-2,:]
            U[t,:,0] = U[t,:,1]
            U[t,:,-1] = U[t,:,-2]
            for x in range(1, self.nx-1):
                for y in range(1, self.ny-1):
                    if x != int(self.nx/2) or y != int(self.ny/2):
                        uxx = (1/self.h**2) * (U[t,x+1,y] - 2*U[t,x,y] + U[t,x-1,y])
                        uyy = (1/self.h**2) * (U[t,x,y+1] - 2*U[t,x,y] + U[t,x,y-1])
                        U[t+1,x,y] = 2*U[t,x,y] - U[t-1,x,y] + self.k**2*uxx + self.k**2*uyy
        return U
    
    
    def max_amplitude(self):
        U = self.solve_matrix()
        max_amp = np.max(np.abs(U))
        return max_amp
    
    
    def visualise(self):
        
        # In this section I am showing the results        

        U = self.solve_matrix()
        Xg = self.mesh()[0]
        Yg = self.mesh()[1]

        ################

        Ua = self.compare()
        
        #this is the analytical solution for the resonance pattern, for comparison

        ################


        fig = plt.figure()
        ax1 = fig.add_subplot(221, projection='3d')
        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')
        ax1.set_zlabel('U')


        ax2 = fig.add_subplot(222)
        ax2.set_xlabel('X')
        ax2.set_ylabel('Y')
        

        ax3 = fig.add_subplot(223)
        ax3.set_xlabel('X')
        ax3.set_ylabel('Y')
        
        ax4 = fig.add_subplot(224)
        ax4.set_xlabel('X')
        ax4.set_ylabel('Y')
        
        ax1.set_box_aspect([1, 1, 1])  
        ax2.set_aspect('equal')  
        ax3.set_aspect('equal')  
        ax4.set_aspect('equal')
        
        colour = 'coolwarm'

 
    
        norm = matplotlib.colors.Normalize(vmin=-self.rplt, vmax=self.rplt)
        
        fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=colour), ax=ax3, orientation='vertical', label='Displacement (m)')
      

        
        ax4.contour(Xg, Yg, Ua, levels=[0], cmap=colour, vmin=-self.rplt, vmax=self.rplt)
        ax4.set_title('Analytical node pattern solution')
        
        for i in range(int(self.pv*(self.t_end / self.k)) ,int(self.t_end / self.k) + 1):
            ax1.clear()
            ax2.clear()
            ax3.clear()
            ax1.plot_surface(Xg, Yg, U[i], cmap=colour, vmin=-self.rplt, vmax=self.rplt)
            ax1.set_zlim([-2, 2])  # set the z-axis limits
            
            #contour plot of 0 displacement
    
            ax2.contour(Xg, Yg, U[i], levels=[0], cmap=colour, vmin=-self.rplt, vmax=self.rplt)
            ax2.set_title('Zero displacement')
            
            ax3.imshow(U[i], interpolation='bilinear', norm=norm, extent=[-self.L/2, self.L/2, -self.L/2, self.L/2], origin='lower', cmap=colour)

            plt.pause(0.000001)

            if keyboard.is_pressed('e'):  # if key 'e' is pressed - press to stop erroneous plots
                plt.close()
                print('Abort!')
                break  # finishing the loop
        
        
        plt.show()



############### In this section I am setting the parameters for the simulation ##################

L = 1
#domain width

h = 0.02
#mesh step h

k = 0.01
#timestep k

amp = 1
#this is the amplitude of the driven oscillation

p = 6.70
#proportionality factor of the frequency-mode relationship

q = 100
#proportionality factor between k/w and t_end

rplt = 3
#this is the displacement axis plot range

pv = 0.8
#this is how far through the timespan the visualisation is started
#i.e 0.8 means only show the last 20% of the timespan



#Mode number input

m = 1

n = 2

#'''

test_plate = Plate(L,h,k,amp,m,n,p,pv,rplt,q)

test_plate.visualise()

print('p =',p)
print('wn =', test_plate.freq,'rad/s')
print('nt =',test_plate.nt)
print('t_end =',test_plate.t_end)


#'''


'''

#This is a section of code we used to find the value of p which led to maximum amplitude oscillations of the correct mode
#This could have been automated further but since it was a one-off calculation we manually set the bounds and step size
#The resonance pattern shape must also be observed to ensure the correct mode is being excited

dp = 0.05
pmin = 6.0
pmax = 7.0
p = np.arange(pmin,pmax,dp)
A = []

for ip in range(len(p)):
    test_plate = Plate(L,h,k,amp,m,n,p[ip],pv,rplt,q)
    amax = test_plate.max_amplitude()
    A.append(amax)
    
    print('factor =',p[ip])
    print('w =', test_plate.freq,'rad/s')
    print('Max amplitide =',amax)
    print('')

plt.plot(p,A)

plt.xlabel('p')
plt.ylabel('Maximum amplitude')
plt.show()

'''
