#import all necessary libraries
import rebound
import math

#initialize the simulation and setup the basic parameters
sim = rebound.Simulation()
sim.G = 1.0
sim.integrator = "ias15"    #Method of calculating
sim.units = ('yr', 'AU', 'Msun')

#add the particles and their parameters to the simulation
sim.add(m=1.0, x=-5.0, vy=3.0,r=0.005) #Star 0
sim.add(m=1.0, x=5.0, vy=-3.0,r=0.005) #Star 1
sim.add(m=1.0, y=5.0, vx=0.001, vy=-0.1,r=0.005) #Star 2

#abbreviation for convenience
sp = sim.particles

#move the reference frame to the center of mass of the system
sim.move_to_com()

#functions for calculating the distance between particles
def gendist1(particle):
    distance1 = abs(particle[0].x - particle[1].x)
    distance2 = abs(particle[0].y - particle[1].y)
    return math.sqrt(distance1**2 + distance2**2)
def gendist2(particle):
    distance1 = abs(particle[0].x - particle[2].x)
    distance2 = abs(particle[0].y - particle[2].y)
    return math.sqrt(distance1**2 + distance2**2)
def gendist3(particle):
    distance1 = abs(particle[1].x - particle[2].x)
    distance2 = abs(particle[1].y - particle[2].y)
    return math.sqrt(distance1**2 + distance2**2)

#function for calculating the general velocity of a particle
def genvel(particle):
    xvel = particle.vx
    yvel = particle.vy
    return math.sqrt(xvel**2 + yvel**2)

#set "max" time for the simulation, but will be extended indefinitely until a collision or escape occurs
max_time = 100

#loop to run the simulation, and check for collisions or ejections at each step
while sim.t < max_time:
    sim.integrate(sim.t + 0.01)
    max_time = sim.t + 0.01
    #collision checks
    if gendist1(sp) <= (sp[0].r) + (sp[1].r):
        print(sim.t)
        sim.stop()
        break
    if gendist2(sp) <= (sp[0].r) + (sp[2].r):
        print(sim.t)
        sim.stop()
        break
    if gendist3(sp) <= (sp[1].r) + (sp[2].r):
        print(sim.t)
        sim.stop()
        break
    #ejection checks
    if genvel(sp[0]) >= math.sqrt(2 * sim.G * 2*(sp[1].m) / ((gendist1(sp) + gendist2(sp))/2)):
        print(sim.t)
        sim.stop()
        break
    if genvel(sp[1]) >= math.sqrt(2 * sim.G * 2*(sp[0].m) / ((gendist1(sp) + gendist3(sp))/2)):
        print(sim.t)
        sim.stop()
        break
    if genvel(sp[2]) >= math.sqrt(2 * sim.G * 2*(sp[0].m) / ((gendist2(sp) + gendist3(sp))/2)):
        print(sim.t)
        sim.stop()
        break

rebound.OrbitPlot(sim, color=True)