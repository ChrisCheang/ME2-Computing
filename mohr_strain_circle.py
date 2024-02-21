import math
from math import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from sympy import Point, Segment


class Rosette:
    def __init__(self, theta1, theta2, theta3, str1, str2, str3):
        if (str1 < str2 and str3 < str2) or (str1 > str2 and str3 > str2):
            print("Make sure strain 2 is the interim strain")
            return
        if theta1 > theta2 or theta2 > theta3 or theta1 > theta3:
            print("Make sure angles are in increasing order")
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
    
    def draw_mohr_circle(self):
        plt.axvline(x=self.str1, linestyle="--")
        plt.axvline(x=self.str2, linestyle="--")
        plt.axvline(x=self.str3, linestyle="--")
        
        theta1 = self.thetaa-self.thetab
        theta2 = self.thetac-self.thetab

        pointA = Point(self.str1,(self.str2 - self.str1)/tan(theta1))
        pointB = Point(self.str2, 0)
        pointC = Point(self.str3, (self.str2 - self.str3)/tan(theta2))

        segmentA = Segment(pointA,pointB)
        segmentB = Segment(pointB,pointC)

        perpBisecA = segmentA.perpendicular_bisector() 
        perpBisecB = segmentB.perpendicular_bisector() 

        O = perpBisecA.intersection(perpBisecB)

        yoffset = -O[0].y

        plt.plot([self.str2,self.str1],[yoffset, yoffset + (self.str2 - self.str1)/tan(theta1)], color = "black", linestyle="--")
        plt.plot([self.str2,self.str3],[yoffset, yoffset + (self.str2 - self.str3)/tan(theta2)], color = "black", linestyle="--")

        Amid = [0.5*(self.str2 + self.str1),0.5*(yoffset + yoffset + (self.str2 - self.str1)/tan(theta1))]
        Bmid = [0.5*(self.str2 + self.str3),0.5*(yoffset + yoffset + (self.str2 - self.str3)/tan(theta2))]

        plt.plot([Amid[0],O[0].x],[Amid[1],0])
        plt.plot([Bmid[0],O[0].x],[Bmid[1],0])

        plt.axis('square')
        plt.show()


rosette789 = Rosette(0,pi/3,135*pi/180,300,130,-50)


rosette789.draw_mohr_circle()

