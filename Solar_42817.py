# -*- coding: utf-8 -*-
"""
Created on Tue Mar 07 20:12:52 2017

@author: Daniel C. Keylon

This python file simulates the orbits of the planetary bodies and the sun
using newtonian 2-body equations.  This file depends on python 2.7, vpython 
classic, scipy.integrate, scipy.constants, pylab, and the future version of
division.  If you try to run this file without any of these, the simulation
probably won't run correctly.

This project was initially completed for a college astrophysics course.  The 
cohesion and readability of the project reflected my poor programming skills.
This can be seen in solar_old.py, which is the original version.  This new
version reorganizes the project for readability and simplicity.  I've removed
copious redundancy through the implementation of a basic Class structure.

In the future, I plan on reducing run time so that
more datapoints can be used.  
"""

from __future__ import division
from pylab import *
from visual import *

from scipy.integrate import odeint
from scipy.constants import G

class Two_Body_Orbit:
    #The class creates objects that automatically calculates orbits based on
    #known parameters of the body's orbit
    
    def __init__(self, M_body, M_attractor, body_state_vector, time_vector):
        #
        self._M1 = M_body
        self._M2 = M_attractor
        self._init_vector = body_state_vector
        self._time = time_vector
        
        self.calc_orbit()
        
        
    def newton_acc(self, state_vector, time_vector):
        #Calculates the acceleration state vector using Newton's gravity law
        # F = G (M1 * M2)/r^2 * r^hat
        k = G * self._M1 * self._M2
        
        x = state_vector[0]
        y = state_vector[2]
        z = state_vector[4]
        r2 = x**2 + y**2 + z**2
        
        self._acc_vector = [state_vector[1], -k/self._M1 * x/(r2**(3/2)), 
                            state_vector[3], -k/self._M1 * y/(r2**(3/2)), 
                            state_vector[5], -k/self._M1 * z/(r2**(3/2))]
        
        return self._acc_vector
    
    def calc_orbit(self):   
        self._orbit = odeint(self.newton_acc, self._init_vector, self._time)
        
    def get_orbit(self):
        return self._orbit

#number of years simulated
years = 200
#rate at which the simulation progresses
rat = 5000
#total number of seconds
timespan = 3600 * 24 * 365 * years
#number of datapoints to calculate
num = 400000
t = linspace(0, timespan, num)

#masses of solar bodies in kilograms
masses = {"sun": 2.0e30, "mer": 3.3e23, "ven": 4.87e24, 
          "ear": 6.0e24, "mar": 6.42e23, "jup": 1.899e27,
          "sat": 5.7e26, "ura": 8.68e25, "nep": 1.0242e26}

#equatorial radii of solar bodies in meters
radii = {"sun": 6.96e8, "mer": 2.44e6, "ven": 6.05e6, 
          "ear": 6.37e6, "mar": 3.40e6, "jup": 7.15e7,
          "sat": 6e7, "ura": 2.56e7, "nep": 2.48e7}

#initial body positions in 3D space
init_pos = {"sun": vector(-7e8, 0, 0), "mer": vector(4.6e10, 0, 0), 
            "ven": vector(1.07e11, 0, 0), "ear": vector(1.47e11, 0, 0), 
            "mar": vector(2.07e11, 0, 0), "jup": vector(7.41e11, 0, 0),
            "sat": vector(1.35e12, 0, 0), "ura": vector(2.7e12, 0, 0), 
            "nep": vector(4.46e12, 0, 0)}

#initial velocities
#possibly an unuseful dictionary, but good for reference if nothing else
init_vel = {"sun": [0, -13808.0476, 0], "mer": [0, 59127.22, 7265], 
            "ven": [0, 35472.14, 2103], "ear": [0, 29366, 0], 
            "mar": [0, 26529.25, 856.88], "jup": [0, 13740.77, 313],
            "sat": [0, 10214.84, 442.4], "ura": [0, 7234.22, 97.22], 
            "nep": [0, 5493.29, 169.75]}

#Scene Generation: generates the visual planets

#scale factor makes small planets visible
#in the scale of the solar system, all planets are small
planet_scale = 1000

sun = sphere(pos = init_pos["sun"], radius = radii["sun"] * 20,
             color = color.yellow, material = materials.emissive,
             make_trail = true)
local_light(pos = sun.pos, color = color.yellow)

mercury = sphere(pos = init_pos["mer"], radius = radii["mer"] * planet_scale,
               color = (1, .7, .2) , make_trail = True)

venus = sphere(pos = init_pos["ven"], radius = radii["ven"] * planet_scale,
               color = color.yellow , make_trail = True)

earth = sphere(pos = init_pos["ear"], radius = radii["ear"] * planet_scale,
               material = materials.earth, make_trail = True)

mars = sphere(pos = init_pos["mar"], radius = radii["mar"] * planet_scale,
               color = color.red , make_trail = True)

jupiter = sphere(pos = init_pos["jup"], radius = radii["jup"] * planet_scale,
               color = color.white , make_trail = True)

saturn = sphere(pos = init_pos["sat"], radius = radii["sat"] * planet_scale,
               color = color.white , make_trail = True)

uranus = sphere(pos = init_pos["ura"], radius = radii["ura"] * planet_scale,
               color = color.cyan , make_trail = True)

neptune = sphere(pos = init_pos["nep"], radius = radii["nep"] * planet_scale,
               color = color.blue , make_trail = True)

#initial state vectors of the form [x, v_x, y, v_y, z, z_y]
q_sun = [sun.pos[0], 0, sun.pos[1], -13808.0476,
       sun.pos[2], 0] 

q_mer = [mercury.pos[0], 0, mercury.pos[1], 59127.22,
       mercury.pos[2], 7265] 

q_ven = [venus.pos[0], 0, venus.pos[1], 35472.14,
       venus.pos[2], 2103] 

q_ear = [earth.pos[0], 0, earth.pos[1], 29366,
       earth.pos[2], 0] 

q_mar = [mars.pos[0], 0, mars.pos[1], 26529.25,
       mars.pos[2], 856.88] 

q_jup = [jupiter.pos[0], 0, jupiter.pos[1], 13740.77,
       jupiter.pos[2], 313]

q_sat = [saturn.pos[0], 0, saturn.pos[1], 10214.84,
       saturn.pos[2], 442.4] 

q_ura = [uranus.pos[0], 0, uranus.pos[1], 7234.22,
       uranus.pos[2], 97.22] 

q_nep = [neptune.pos[0], 0, neptune.pos[1], 5493.29,
       neptune.pos[2], 169.75]

#orbit calculations using class declaration
sun_orbit = Two_Body_Orbit(masses["sun"],masses["jup"], q_sun, t).get_orbit()

mer_orbit = Two_Body_Orbit(masses["mer"],masses["sun"], q_mer, t).get_orbit()

ven_orbit = Two_Body_Orbit(masses["ven"],masses["sun"], q_ven, t).get_orbit()

ear_orbit = Two_Body_Orbit(masses["ear"],masses["sun"], q_ear, t).get_orbit()

mar_orbit = Two_Body_Orbit(masses["mar"],masses["sun"], q_mar, t).get_orbit()

jup_orbit = Two_Body_Orbit(masses["jup"],masses["sun"], q_jup, t).get_orbit()

sat_orbit = Two_Body_Orbit(masses["sat"],masses["sun"], q_sat, t).get_orbit()

ura_orbit = Two_Body_Orbit(masses["ura"],masses["sun"], q_ura, t).get_orbit()

nep_orbit = Two_Body_Orbit(masses["nep"],masses["sun"], q_nep, t).get_orbit()

#Position updating
for i in range(0, num):
    rate(rat)

    sun.pos = [sun_orbit[i][0], sun_orbit[i][2], sun_orbit[i][4]]
    
    mercury.pos = [mer_orbit[i][0], mer_orbit[i][2], mer_orbit[i][4]]

    venus.pos = [ven_orbit[i][0], ven_orbit[i][2], ven_orbit[i][4]]
    
    earth.pos = [ear_orbit[i][0], ear_orbit[i][2], ear_orbit[i][4]]

    mars.pos = [mar_orbit[i][0], mar_orbit[i][2], mar_orbit[i][4]]

    jupiter.pos = [jup_orbit[i][0], jup_orbit[i][2], jup_orbit[i][4]]

    saturn.pos = [sat_orbit[i][0], sat_orbit[i][2], sat_orbit[i][4]]

    uranus.pos = [ura_orbit[i][0], ura_orbit[i][2], ura_orbit[i][4]]

    neptune.pos = [nep_orbit[i][0], nep_orbit[i][2], nep_orbit[i][4]]