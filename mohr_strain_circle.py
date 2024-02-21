import math
from math import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


class Rosette:
    def __init__(self, theta1, theta2, theta3, str1, str2, str3):
        self.thetaa = theta1
        self.thetab = theta2
        self.thetac = theta3
        self.str1 = str1
        self.str2 = str2
        self.str3 = str3

    def solve(self):
        M = np.zeros((3,3))
        M[0,:] = [0.5*(1+cos(2*self.thetaa)),0.5*(1-cos(2*self.thetaa)),0.5*sin(2*self.thetaa)]
        M[1,:] = [0.5*(1+cos(2*self.thetab)),0.5*(1-cos(2*self.thetab)),0.5*sin(2*self.thetab)]
        M[2,:] = [0.5*(1+cos(2*self.thetac)),0.5*(1-cos(2*self.thetac)),0.5*sin(2*self.thetac)]

        s = np.zeros((3,1))
        s[:,0] = [self.str1, self.str2, self.str3]

        sol = np.matmul(np.linalg.inv(M),s)
        return [sol[i,0] for i in range(3)]
    
    def draw_mohr_circle(self):
        plt.axvline(x=self.str1)
        plt.axvline(x=self.str2)
        plt.axvline(x=self.str3)
        #plt.set_aspect('equal', 'box')
        plt.show()



rosette789 = Rosette(0,pi/3,135*pi/180,300,130,-50)


rosette789.draw_mohr_circle()
print(rosette789.solve())
