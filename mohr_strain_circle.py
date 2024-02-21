import math
from math import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from sympy import Point, Segment


def volt_to_strain_gauge(mVolt):
    return mVolt*1000*4/(2.1*2.5)


class Rosette:

    def __init__(self, theta1, theta2, theta3, str1, str2, str3):
        if (str1 < str2 and str3 < str2) or (str1 > str2 and str3 > str2):
            print("Make sure strain 2 is the interim strain")
            return
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
    
    def solve_principal(self):
        e1 = 0.5*(self.solve()[0] + self.solve()[1]) + 0.5*sqrt((self.solve()[0] - self.solve()[1])**2 + self.solve()[2]**2)
        e2 = 0.5*(self.solve()[0] + self.solve()[1]) - 0.5*sqrt((self.solve()[0] - self.solve()[1])**2 + self.solve()[2]**2)

        return [e1,e2]
    
    def draw_mohr_circle(self):
        plt.axvline(x=self.str1, linestyle="--")
        plt.axvline(x=self.str2, linestyle="--")
        plt.axvline(x=self.str3, linestyle="--")
        
        theta1 = self.thetaa-self.thetab
        theta2 = self.thetac-self.thetab

        switcha = abs(cos(theta1))/cos(theta1)
        switchc = abs(cos(theta2))/cos(theta2)

        print(theta1)
        print(theta2)


        pointA = Point(self.str1, switcha*(self.str2 - self.str1)/tan(theta1))
        pointB = Point(self.str2, 0)
        pointC = Point(self.str3, switchc*(self.str2 - self.str3)/tan(theta2))

        segmentA = Segment(pointA,pointB)
        segmentB = Segment(pointB,pointC)

        perpBisecA = segmentA.perpendicular_bisector() 
        perpBisecB = segmentB.perpendicular_bisector() 

        O = perpBisecA.intersection(perpBisecB)[0]

        yoffset = -O.y

        plt.plot([self.str2,self.str1],[yoffset, yoffset + pointA.y], color = "black", linestyle="--")
        plt.plot([self.str2,self.str3],[yoffset, yoffset + pointC.y], color = "black", linestyle="--")
        
        Amid = [0.5*(self.str2 + self.str1),0.5*(2 * yoffset + pointA.y)]
        Bmid = [0.5*(self.str2 + self.str3),0.5*(2 * yoffset + pointC.y)]

        plt.plot([Amid[0],O.x],[Amid[1],0], color = "green", linestyle="--")
        plt.plot([Bmid[0],O.x],[Bmid[1],0], color = "green", linestyle="--")

        theta = np.linspace( 0 , 2 * np.pi , 150 )
 
        radius = O.distance(pointA.midpoint(pointB))
 
        a = radius * np.cos( theta ) + O.x
        b = radius * np.sin( theta )

        plt.plot(a,b)

        plt.axis('square')
        plt.show()


milliVolts = [-0.05919,-0.04729,0.35239]

strains = [volt_to_strain_gauge(milliVolts[i]) for i in range(3)]

angles = [0 + 20*pi/180,pi/4 + 20*pi/180,pi/2 + 20*pi/180]


rosettetest = Rosette(0,4*pi/3,2*pi/3,108,90,64)   
rosettetest2 = Rosette(0,pi/3,135*pi/180,300,130,-50)   
rosette789 = Rosette(angles[0],angles[1],angles[2],strains[2],strains[1],strains[0])


rosettetest.draw_mohr_circle()
#rosette789.draw_mohr_circle()

#print(rosette789.solve())
