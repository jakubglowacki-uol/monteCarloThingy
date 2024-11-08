from SimulatorUnit import SimulatorUnit 
from typing import List
import matplotlib.pyplot as plt
from rdfpy import rdf
import time

epsilon = 128000
r0 = 342  

def current_milli_time():
    return round(time.time() * 1000)

def initSim():
    
    print("Initializing simulation")
    steps = 6
    
    #create simulator cube
    newSim = SimulatorUnit([], 24, 0)
    #populate cube with molecules
    newSim.populateArbitrary(1, 1, 2, newSim.L/2)
     #create figure
    fig = plt.figure(1)
    fig2 = plt.figure(2)
    fig3 = plt.figure(3)
    ax = fig2.add_subplot(projection='3d')
    ab = fig.add_subplot(projection='3d')
    #ac = fig3.add_subplot()
    #ac.set_title("RDF")
    ax.set_title("Argon after "+ str(steps)+" steps, at " +str(newSim.T) + "K, with side length " + str(newSim.L) + ", and " + str(len(newSim.molecules)) + " molecules")
    ab.set_title("Argon before "+ str(steps)+" steps, at " +str(newSim.T) + "K, with side length " + str(newSim.L) + ", and " + str(len(newSim.molecules)) + " molecules")
    xb = [molecule.x for molecule in newSim.molecules]
    yb = [molecule.y for molecule in newSim.molecules]
    zb = [molecule.z for molecule in newSim.molecules]
    
    accepted = 0
    notaccepted = 0
    time1 = current_milli_time()
    for i in range(steps):
        stepAccepted = newSim.MonteCarloStep()
        if stepAccepted:
            accepted += 1
        else:
            notaccepted += 1
    time2 = current_milli_time()
    print("Time taken: " + str(time2-time1) + "ms")
    print("Time per step: " + str((time2-time1)/steps) + "ms")
    xs = [molecule.x for molecule in newSim.molecules]
    ys = [molecule.y for molecule in newSim.molecules]
    zs = [molecule.z for molecule in newSim.molecules]
    
    ax.scatter(xs, ys, zs)
    ab.scatter(xb, yb, zb)
    ax.set_xlim(0,newSim.L)
    ab.set_xlim(0,newSim.L)
    ax.set_ylim(0,newSim.L)
    ab.set_ylim(0,newSim.L)
    ax.set_zlim(0,newSim.L)
    ab.set_zlim(0,newSim.L)
    
    print("Acceptance Rate: " + str(accepted/(accepted+notaccepted)*100) +"%")
    print("Average Potential Energy: " + str(newSim.AvgPotentialEnergy(newSim.molecules)/epsilon) + " kJ/mol")
    print("Standard Deviation Of Potential Energies: " + str(newSim.stdDevPotentialEnergy(newSim.molecules)/epsilon) + " kJ/mol")
    #12.386jk-1mol-1
    print("Specific Heat Capacity: " + str(SpecHeatCap(newSim.AvgPotentialEnergy(newSim.molecules)/epsilon, newSim.stdDevPotentialEnergy(newSim.molecules)/epsilon, 300)))
    #g_r, r = rdf(newSim.getMoleculeCoors(), dr=0.05, parallel=False)
    #ac.plot(r, g_r)
    plt.show()

def SpecHeatCap(avg: float, stdev:float, T:float):
        kB=1.38064852e-23
        V = 3
        beta2 = 1/(kB*(T**2))
        kinetic=(V/2)*kB
        Cv=((beta2)*(stdev) - ((avg)**2)) + kinetic
        return Cv

initSim()