from __future__ import division
from pylab import *
from visual import *

from scipy.integrate import odeint
from scipy.constants import G

#Timespan
years = float(raw_input("How many years should I run? Enter a positive number (A few years are best for inner planets, up to 164ish best for outer planets)"))
rat = int(raw_input("Should I run fast or slow? Enter a positive integer.  Small is slow, big is fast."))
timespan = 3600 * 24 * 365 * years
num = int(raw_input("How many data points should I calculate? 10000 is a good baseline. More makes it take longer to calculate."))
t = linspace(0, timespan, num)

#Sun
#6.96e8 * 20
sun = sphere(pos = vector(-7e8, 0, 0), radius = 6.96e8 * 20 * (1/20),
             color = color.yellow, material = materials.emissive,
             make_trail = true)
local_light(pos = sun.pos, color = color.yellow)

q_C = [sun.pos[0], 0, sun.pos[1], -13808.0476,
       sun.pos[2], 0] 


def fC(q_i, t_i):
    m = 1.899e27
    M = 2.0e30
    k = G * M * m
    
    x = q_i[0]
    y = q_i[2]
    z = q_i[4]
    r2 = x**2 + y**2 + z**2
    
    return [q_i[1], -k/M * x/(r2**(3/2)), q_i[3], -k/M * y/(r2**(3/2)), q_i[5],
            -k/M * z/(r2**(3/2))]

C = odeint(fC, q_C, t)

#Earth
earth = sphere(pos = vector(1.47e11, 0, 0), radius = 6.37e6 * 1000,
               material = materials.earth, make_trail = True)

q_E = [earth.pos[0], 0, earth.pos[1], 29366,
       earth.pos[2], 0] 


def fE(q_i, t_i):
    m = 6.0e24
    M = 2.0e30
    k = G * M * m
    
    x = q_i[0]
    y = q_i[2]
    z = q_i[4]
    r2 = x**2 + y**2 + z**2
    
    return [q_i[1], -k/m * x/(r2**(3/2)), q_i[3], -k/m * y/(r2**(3/2)), q_i[5],
            -k/m * z/(r2**(3/2))]

E = odeint(fE, q_E, t)

#Mercury
mercury = sphere(pos = vector(4.6e10, 0, 0), radius = 2439700 * 1000,
               color = (1, .7, .2) , make_trail = True)

q_Me = [mercury.pos[0], 0, mercury.pos[1], 59127.22,
       mercury.pos[2], 7265] 


def fMe(q_i, t_i):
    m = 3.3e23
    M = 2.0e30
    k = G * M * m
    
    x = q_i[0]
    y = q_i[2]
    z = q_i[4]
    r2 = x**2 + y**2 + z**2
    
    return [q_i[1], -k/m * x/(r2**(3/2)), q_i[3], -k/m * y/(r2**(3/2)), q_i[5],
            -k/m * z/(r2**(3/2))]

Me = odeint(fMe, q_Me, t)

#Venus
venus = sphere(pos = vector(1.07e11, 0, 0), radius = 6051800 * 1000,
               color = color.yellow , make_trail = True)

q_V = [venus.pos[0], 0, venus.pos[1], 35472.14,
       venus.pos[2], 2103] 


def fV(q_i, t_i):
    m = 4.87e24
    M = 2.0e30
    k = G * M * m
    
    x = q_i[0]
    y = q_i[2]
    z = q_i[4]
    r2 = x**2 + y**2 + z**2
    
    return [q_i[1], -k/m * x/(r2**(3/2)), q_i[3], -k/m * y/(r2**(3/2)), q_i[5],
            -k/m * z/(r2**(3/2))]

V = odeint(fV, q_V, t)

#Mars
mars = sphere(pos = vector(2.07e11, 0, 0), radius = 3396200 * 1000,
               color = color.red , make_trail = True)

q_Ma = [mars.pos[0], 0, mars.pos[1], 26529.25,
       mars.pos[2], 856.88] 


def fMa(q_i, t_i):
    m = 6.42e23
    M = 2.0e30
    k = G * M * m
    
    x = q_i[0]
    y = q_i[2]
    z = q_i[4]
    r2 = x**2 + y**2 + z**2
    
    return [q_i[1], -k/m * x/(r2**(3/2)), q_i[3], -k/m * y/(r2**(3/2)), q_i[5],
            -k/m * z/(r2**(3/2))]

Ma = odeint(fMa, q_Ma, t)

#Jupiter
jupiter = sphere(pos = vector(7.41e11, 0, 0), radius = 71492000 * 1000,
               color = color.white , make_trail = True)

q_J = [jupiter.pos[0], 0, jupiter.pos[1], 13740.77,
       jupiter.pos[2], 313] 


def fJ(q_i, t_i):
    m = 1.899e27
    M = 2.0e30
    k = G * M * m
    
    x = q_i[0]
    y = q_i[2]
    z = q_i[4]
    r2 = x**2 + y**2 + z**2
    
    return [q_i[1], -k/m * x*(r2**(2)), q_i[3], -k/m * y*(r2**(2)), q_i[5],
            -k/m * z*(r2**(2))]

J = odeint(fJ, q_J, t)

#Saturn
saturn = sphere(pos = vector(1.35e12, 0, 0), radius = 60000000 * 1000,
               color = color.white , make_trail = True)

q_S = [saturn.pos[0], 0, saturn.pos[1], 10214.84,
       saturn.pos[2], 442.4] 


def fS(q_i, t_i):
    m = 5.7e26
    M = 2.0e30
    k = G * M * m
    
    x = q_i[0]
    y = q_i[2]
    z = q_i[4]
    r2 = x**2 + y**2 + z**2
    
    return [q_i[1], -k/m * x/(r2**(3/2)), q_i[3], -k/m * y/(r2**(3/2)), q_i[5],
            -k/m * z/(r2**(3/2))]

S = odeint(fS, q_S, t)

#Uranus
uranus = sphere(pos = vector(2.7e12, 0, 0), radius = 25559000 * 1000,
               color = color.cyan , make_trail = True)

q_U = [uranus.pos[0], 0, uranus.pos[1], 7234.22,
       uranus.pos[2], 97.22] 


def fU(q_i, t_i):
    m = 8.68e25
    M = 2.0e30
    k = G * M * m
    
    x = q_i[0]
    y = q_i[2]
    z = q_i[4]
    r2 = x**2 + y**2 + z**2
    
    return [q_i[1], -k/m * x/(r2**(3/2)), q_i[3], -k/m * y/(r2**(3/2)), q_i[5],
            -k/m * z/(r2**(3/2))]

U = odeint(fU, q_U, t)

#Neptune
neptune = sphere(pos = vector(4.46e12, 0, 0), radius = 24764000 * 1000,
               color = color.blue , make_trail = True)

q_N = [neptune.pos[0], 0, neptune.pos[1], 5493.29,
       neptune.pos[2], 169.75] 


def fN(q_i, t_i):
    m = 1.0243e26
    M = 2.0e30
    k = G * M * m
    
    x = q_i[0]
    y = q_i[2]
    z = q_i[4]
    r2 = x**2 + y**2 + z**2
    
    return [q_i[1], -k/m * x/(r2**(3/2)), q_i[3], -k/m * y/(r2**(3/2)), q_i[5],
            -k/m * z/(r2**(3/2))]

N = odeint(fN, q_N, t)

#Position updating
i = 0
while i < num:
    rate(rat)

    sun.pos = [C[i][0], C[i][2], C[i][4]]

    earth.pos = [E[i][0], E[i][2], E[i][4]]
    earth.rotate(angle = (2*pi)/(3600*24), axis = earth.axis, origin = earth.pos)
    
    mercury.pos = [Me[i][0], Me[i][2], Me[i][4]]

    venus.pos = [V[i][0], V[i][2], V[i][4]]

    mars.pos = [Ma[i][0], Ma[i][2], Ma[i][4]]

    jupiter.pos = [J[i][0], J[i][2], J[i][4]]

    saturn.pos = [S[i][0], S[i][2], S[i][4]]

    uranus.pos = [U[i][0], U[i][2], U[i][4]]

    neptune.pos = [N[i][0], N[i][2], N[i][4]]
    
    i = i + 1
    


