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
        return atan2(self.im,self.re)

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
print("Q1 Ans = at time ", t[index], "s, sum is ", avi[index])


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
    #j01w = Complex(0, 0.1*w)
    #denom = Complex(1,0).add(j01w)
    #return Complex(1,0).divide(denom)

    return 1/(1+1j*0.1*w)



wlog = np.linspace(-3,5,301)  #log x
w = [10**x for x in wlog]


biH = [Bi(x) for x in w]
biHa = [20*log10(np.abs(biH[i])) for i in range(len(w))]
biHp = [np.angle(biH[i])*180/pi for i in range(len(w))]


#plt.plot(wlog, biHa)
#plt.plot(wlog,biHp)
#plt.legend(['amplitude', 'phase'])
#plt.show()


def Bii(w):
    R1 = 1000
    R2 = 2000
    C1 = 0.001
    C2 = 0.002

    jwr2c2plus1 = 1+(R2*C2*1j*w)   #Complex(0, R2*C2*w).add(Complex(1, 0))
    jwr1c1plus1 = 1+(R1*C1*1j*w)   #Complex(0, R1*C1*w).add(Complex(1, 0))
    jrcfrac = (R1/R2)*jwr2c2plus1/jwr1c1plus1      # jwr2c2plus1.divide(jwr1c1plus1).multiply(Complex(R1/R2, 0))

    return 1/(1+jrcfrac)

    #return Complex(1, 0).divide(Complex(1, 0).add(jrcfrac))


wlog = np.linspace(-3,1,301)  #log x
w = [10**x for x in wlog]

biiH = [Bii(x) for x in w]

biiHa = [20*log10(np.abs(biiH[i])) for i in range(len(w))]
biiHp = [np.angle(biiH[i])*180/pi for i in range(len(w))]

#plt.plot(wlog, biiHa)
#plt.plot(wlog, biiHp)
#plt.legend(['amplitude', 'phase'])
#plt.show()

w = 0.1
print("Q2 Ans: ", 20*log10(np.abs(Bii(0.1))))


def HLP(w):
    j001w = Complex(0,0.01*w)
    return Complex(1,0).divide(Complex(1,0).add(j001w))


def HHP(w):
    j40w = Complex(0,40*w)
    return j40w.divide(Complex(1,0).add(j40w))


def HPLP(w):
    mod = HHP(w).mod()*HLP(w).mod()
    phase = HHP(w).phase()+HLP(w).phase()
    return Complex(mod*cos(phase), mod*sin(phase))


wlog = np.linspace(-5,5,301)  #log x
w = [10**x for x in wlog]

biv = [HPLP(x) for x in w]

biva = [20*log10(biv[i].mod()) for i in range(len(w))]
bivp = [biv[i].phase()*180/pi for i in range(len(w))]

#plt.plot(wlog, biva)
#plt.plot(wlog, bivp)
#plt.legend(['amplitude','phase'])
#plt.show()


# Task C: Fourier Series


def ci(t, T, N):
    y = 1/2
    for n in range(1,N+1):
        y -= sin(2*n*pi*t/T) / (pi * n)
    return y


T = 2
N = 50

t = np.linspace(0,2*T,101)
yi = [ci(i, T, N) for i in t]

#plt.plot(t,yi)
#plt.show()


def cii(t, T, N):
    y = 0
    for n in [n for n in range(1, N+1) if n % 2 == 1]:
        y += 4 * sin(2*n*pi*t/T) / (n * pi)
    return y


T = 5
N = 20     # note: changing this sometimes results in a very smooth wave, but this is because of insufficient sampling i.e. Nyquist criterion not reached.

t = np.linspace(0,2*T,201)
yi = [cii(i, T, N) for i in t]

#plt.plot(t,yi)
#plt.show()

print("Q3 ans: ", cii(0.1,T,N))

# Task D: Discrete Fourier Transform


def DFT(yn):
    N = len(yn)
    FTk = np.zeros(N).astype(complex)
    for k in range(N):
        sum = 0
        for n in range(N):
            sum += yn[n]*math.e**(-2*pi*1j*k*n/N)
        FTk[k] = sum
    return FTk   # note only the front half spectrum is valid.


def DFTInv(FTk):
    N = len(FTk)
    yn = np.zeros(N).astype(complex)
    for n in range(N):
        sum = 0
        for k in range(N):
            sum += FTk[k]*e**(2*pi*1j*k*n/N)
        yn[n] = sum/N
    return yn

dt = 6*pi/19
N = 20
df = 1/(N*dt)

t = np.linspace(0,6*pi,20)
yn = [e**(-(x-5)**2/4) for x in t]
ynft = np.fft.fft(yn)#
ynftinv = DFTInv(ynft)   #np.fft.ifft(ynft)#

index = 5
print("Q4 Ans = at f ", t[index]/(N*dt), "Hz, sum is ", ynft[index].real)

#plt.plot(t,yn)
#plt.plot(t,ynft.real)
#plt.plot(t,ynftinv.real)
#plt.show()




# Task E: Signal Processing


r = open("Vibration.txt", "r")
r_read = r.readlines()
vibrations = []
for item in r_read:
    vibrations.append(float(item.rstrip()))
r.close()

dt = 0.01
N = len(vibrations)
df = 1/(N*dt)

t = np.arange(0, dt*N, dt)
ftVibrations = DFT(vibrations)

resonantFreqIndex = np.where(ftVibrations == ftVibrations.max())[0][0]
resonantFreq = resonantFreqIndex * df   # ! actually check if this method is correct! - the time x axis should translate directly to the frequencies.

#print("Resonant Freq = ", resonantFreq)

#plt.plot(t, fftvibrations)
#plt.show()


r = open("Noisy.txt", "r")
r_read = r.readlines()
noise = []
for item in r_read:
    noise.append(float(item.rstrip()))
r.close()

dt = 1/20
N = len(noise)
df = 1/(N*dt)

t = np.arange(0, dt*N, dt)
noiseFT = DFT(noise)
noiseFTCleaned = [0 if (0.5 < t[i] < 19.5) else noiseFT[i] for i in range(N)]
noiseCleaned = DFTInv(noiseFTCleaned)

index = 100
print("Q5 Ans = at t ", t[index], "s, signal is ", noiseCleaned[index].real)

plt.plot(t, noise)
#plt.plot(t, noiseFT)
plt.plot(t, noiseCleaned)
plt.show()

