import math
from math import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


# Task A: Complex numbers and phasors


tl = np.linspace(0, 2*np.pi, 201)
t = np.linspace(0, 2, 201)

ai1 = [10*cos(pi*x) for x in tl]
ai2 = [5*cos(pi*x + pi/2) for x in tl]

aii = [5*cos(pi*x + pi/2) + 5 for x in tl]

aiii1 = [10*cos(pi*x) for x in t]
aiii2 = [10*cos(2*pi*x) for x in t]

aiv1 = [10*cos(pi*x) for x in t]
aiv2 = [10*cos(2*pi*x - pi/4) for x in t]

#plt.plot(t,aiv1)
#plt.plot(t,aiv2)


class Complex:
    def __init__(self, re, im):
        self.re = re
        self.im = im

    def mod(self):
        return sqrt(self.re**2+self.im**2)

    def phase(self):
        return atan(self.im/self.re)

    def add(self, b):
        return Complex(self.re+b.re,self.im+b.im)

    def multiply(self, b):
        mag = self.mod() * b.mod()
        phase = self.phase() + b.phase()
        return Complex(mag*cos(phase), mag*sin(phase))

    def divide(self, b):
        mag = self.mod() / b.mod()
        phase = self.phase() - b.phase()
        return Complex(mag * cos(phase), mag * sin(phase))

    def plot(self):
        plt.plot([0,self.re],[0,self.im])


av1 = Complex(10*cos(0), 10*sin(0))
av2 = Complex(10*cos(-pi/4), 10*sin(-pi/4))

#plt.plot([0,av1.re],[0,av1.im])
#plt.plot([0,av2.re],[0,av2.im])

avi = [aiv1[i] + aiv2[i] for i in range(len(t))]

#av1.add(av2).plot()
#plt.plot(t, avi)
#plt.show()

index = 20
#print("Q1 Ans = at time ", t[index], "s, sum is ", avi[index])


# Task B: Complex functions: analogue filters and bode plots


class Phasor:
    def __init__(self,magnitude,phase):
        self.mod = magnitude
        self.phase = phase

    def add(self, b):
        mag = self.mod*cos(self.phase) + b.mod*cos(b.phase)
        phase = self.mod*sin(self.phase) + b.mod*sin(b.phase)
        return Phasor(mag,phase)

    def multiply(self, b):
        return Phasor(self.mod*b.mod, self.phase+b.phase)

    def divide(self, b):
        return Phasor(self.mod/b.mod, self.phase-b.phase)


def Bi(w):
    j01w = Complex(0, 0.1*w)
    denom = Complex(1,0).add(j01w)
    return Complex(1,0).divide(denom)


wlog = np.linspace(-3,5,301)  #log x
w = [10**x for x in wlog]


biH = [Bi(x) for x in w]
biHa = [20*log10(biH[i].mod()) for i in range(len(w))]
biHp = [biH[i].phase()*180/pi for i in range(len(w))]


plt.plot(wlog, biHa)
plt.plot(wlog,biHp)
plt.legend(['amplitude', 'phase'])
plt.show()


def Bii(w):
    R1 = 1000
    R2 = 2000
    C1 = 0.001
    C2 = 0.002

    jwr2c2plus1 = Complex(0, R2*C2*w).add(Complex(1, 0))
    jwr1c1plus1 = Complex(0, R1*C1*w).add(Complex(1, 0))
    jrcfrac = jwr2c2plus1.divide(jwr1c1plus1).multiply(Complex(R1/R2, 0))

    return Complex(1, 0).divide(Complex(1, 0).add(jrcfrac))


biiH = [Bii(x) for x in w]

biiHa = [20*log10(biiH[i].mod()) for i in range(len(w))]
biiHp = [biiH[i].phase()*180/pi for i in range(len(w))]

plt.plot(wlog, biiHa)
plt.plot(wlog, biiHp)
plt.legend(['amplitude', 'phase'])
plt.show()

print(Bii(0.0001).mod())
