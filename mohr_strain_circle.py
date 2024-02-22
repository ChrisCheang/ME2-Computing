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
        plt.xlabel("en (microstrain)")
        plt.ylabel("0.5es")

        plt.axvline(x=0, color="black", linestyle="-", linewidth=0.8)

        plt.axvline(x=self.str1, linestyle="--", linewidth=0.8)
        plt.axvline(x=self.str2, linestyle="--", linewidth=0.8)
        plt.axvline(x=self.str3, linestyle="--", linewidth=0.8)

        #plt.annotate('gauge 7 = -47.9', xy=(self.str1, 100), xytext=(self.str1 + 100, 200), arrowprops=dict(arrowstyle="->", facecolor='black'))
        #plt.annotate('gauge 8 = -34.0', xy=(self.str2, 100), xytext=(self.str2 + 100, 100), arrowprops=dict(arrowstyle="->", facecolor='black'))
        #plt.annotate('gauge 9 = 265.8', xy=(self.str3, 100), xytext=(self.str3 - 100, 200), arrowprops=dict(arrowstyle="->", facecolor='black'))
        
        #plt.gca().get_yaxis().set_visible(False)
        #plt.axis('square')
        #plt.show()

        theta1 = self.thetaa-self.thetab
        theta2 = self.thetac-self.thetab

        switcha = -abs(cos(theta1))/cos(theta1)       # this bit is a quick fix that doesn't really work all the time - fix later
        switchc = -abs(cos(theta2))/cos(theta2)


        pointA = Point(self.str1, switcha*(self.str2 - self.str1)/tan(theta1))
        pointB = Point(self.str2, 0)
        pointC = Point(self.str3, switchc*(self.str2 - self.str3)/tan(theta2))

        segmentA = Segment(pointA,pointB)
        segmentB = Segment(pointB,pointC)

        perpBisecA = segmentA.perpendicular_bisector() 
        perpBisecB = segmentB.perpendicular_bisector() 

        O = perpBisecA.intersection(perpBisecB)[0]

        yoffset = -O.y

        #plt.plot([self.str2,self.str1],[yoffset, yoffset + pointA.y], color = "black", linestyle="--", linewidth=0.8)
        #plt.plot([self.str2,self.str3],[yoffset, yoffset + pointC.y], color = "black", linestyle="--", linewidth=0.8)

        theta = np.linspace( 0 , pi/4 , 20 )

        a = 10 * np.cos( theta + pi/2) + self.str2
        b = 10 * np.sin( theta + pi/2) - O.y
        
        #plt.plot(a,b, color="blue", linewidth=0.5)

        a = 20 * np.cos( -theta + pi/2) + self.str2
        b = 20 * np.sin( -theta + pi/2) - O.y
        
        #plt.plot(a,b, color="blue", linewidth=0.5)

        #plt.annotate('45 degrees anticlockwise from gauge 8 (intermediate) to 7', xy=(self.str2 - 5, - O.y + 10), xytext=(-50, 0), arrowprops=dict(arrowstyle="->", facecolor='black'))
        #plt.annotate('45 degrees clockwise from gauge 8 (intermediate) to 9', xy=(self.str2 + 10, - O.y + 15), xytext=(50, -50), arrowprops=dict(arrowstyle="->", facecolor='black'))

        #plt.gca().get_yaxis().set_visible(False)
        #plt.axis('square')
        #plt.show()

        Amid = [0.5*(self.str2 + self.str1),0.5*(2 * yoffset + pointA.y)]
        Bmid = [0.5*(self.str2 + self.str3),0.5*(2 * yoffset + pointC.y)]

        #plt.plot([Amid[0],O.x],[Amid[1],0], color = "green", linestyle="--", linewidth=0.8)
        #plt.plot([Bmid[0],O.x],[Bmid[1],0], color = "green", linestyle="--", linewidth=0.8)

        theta = np.linspace( 0 , 2 * np.pi , 150 )
 
        radius = O.distance(pointA)
 
        a = radius * np.cos( theta ) + O.x
        b = radius * np.sin( theta )

        plt.axhline(y=0, color="black", linestyle="-", linewidth=0.8)

        plt.annotate('O', xy=(O.x, 0), xytext=(O.x - 50, 10), arrowprops=dict(arrowstyle="->", facecolor='black'))

        plt.plot(a,b, linewidth=0.5)

        #plt.axis('square')
        #plt.show()

        plt.axis('square')
        
        #plt.plot([pointC.x],[pointC.y - O.y], 'ro')

        thetaOC = atan2(pointC.y - O.y,pointC.x - O.x)

        plt.plot([radius * np.cos( thetaOC - 40*pi/180) + O.x], [radius * np.sin( thetaOC - 40*pi/180)], 'ro')

    
        plt.annotate('O', xy=(O.x, 0), xytext=(O.x - 50, 10), arrowprops=dict(arrowstyle="->", facecolor='black'))


        

        theta = np.linspace( thetaOC - 40*pi/180 , thetaOC , 200 )

        a = 100 * np.cos( theta) + O.x
        b = 100 * np.sin( theta)


        #plt.plot(a,b, color="black", linewidth=0.5)


        #plt.plot([O.x,pointC.x],[0,pointC.y - O.y], color="black", linewidth=0.5)
        #plt.plot([O.x,radius*np.cos(thetaOC - 40*pi/180) + O.x],[0,radius*np.sin(thetaOC - 40*pi/180)], color="black", linewidth=0.5)


        #plt.annotate('Clockwise by 2 x 20 deg to get to x direction from gauge 9', xy=(O.x+80,50), xytext=(0, -100), arrowprops=dict(arrowstyle="->", facecolor='black'))


        #plt.axvline(x=radius * np.cos( thetaOC - 40*pi/180) + O.x, color="green", linestyle="--", linewidth=0.8)

        #plt.annotate('Clockwise by 2 x 20 deg to get to x direction from gauge 9', xy=(O.x+80,50), xytext=(0, -100), arrowprops=dict(arrowstyle="->", facecolor='black'))








        




milliVolts = [-0.05919,-0.04729,0.35239]

#strains = [volt_to_strain_gauge(milliVolts[i]) for i in range(3)]
strains = [-47.93386667,-33.99817143,265.758019]

offset = 20*pi/180
angles = [0 + offset,pi/4 + offset,pi/2 + offset]


rosettetest = Rosette(0,4*pi/3,2*pi/3,108,90,64)   
rosettetest2 = Rosette(0,pi/3,135*pi/180,300,130,-50)   
rosettetest3 = Rosette(0,pi/3,135*pi/180,300,130,-50)   
rosette789 = Rosette(angles[0],angles[1],angles[2],strains[0],strains[1],strains[2])


#rosettetest.draw_mohr_circle()
rosette789.draw_mohr_circle()
plt.show()

print(rosette789.solve())
print(rosette789.solve_principal())
