from SimulatorUnit import SimulatorUnit 
from typing import List
import matplotlib.pyplot as plt

epsilon = 128
r0 = 342

def initSim():
    print("Initializing simulation")
    steps = 100

    #create simulator cube
    newSim = SimulatorUnit([], 100, 0.5)
    #populate cube with molecules
    newSim.populateSolid(1, 1, 256, newSim.L/2)
     #create figure
    fig = plt.figure(1)
    fig2 = plt.figure(2)
    fig3 = plt.figure(3)
    ax = fig2.add_subplot(projection='3d')
    ab = fig.add_subplot(projection='3d')
    ac = fig3.add_subplot()
    ax.set_title("Argon after "+ str(steps)+" steps, at " +str(newSim.T) + "K, with side length " + str(newSim.L) + ", and " + str(len(newSim.molecules)) + " molecules")
    ab.set_title("Argon before "+ str(steps)+" steps, at " +str(newSim.T) + "K, with side length " + str(newSim.L) + ", and " + str(len(newSim.molecules)) + " molecules")
    ac.set_title("Potential Energy Over Steps")
    xc = [i for i in range(steps)]
    yc = [potEnergies]
    xb = [molecule.x for molecule in newSim.molecules]
    yb = [molecule.y for molecule in newSim.molecules]
    zb = [molecule.z for molecule in newSim.molecules]
    
    accepted = 0
    notaccepted = 0
    for i in range(steps):
        stepAccepted = newSim.MonteCarloStep()
        if stepAccepted:
            accepted += 1
        else:
            notaccepted += 1
    xs = [molecule.x for molecule in newSim.molecules]
    ys = [molecule.y for molecule in newSim.molecules]
    zs = [molecule.z for molecule in newSim.molecules]
    
    ax.scatter(xs, ys, zs)
    ab.scatter(xb, yb, zb)
    ac.plot(xc, yc)
    ax.set_xlim(0,newSim.L)
    ab.set_xlim(0,newSim.L)
    ax.set_ylim(0,newSim.L)
    ab.set_ylim(0,newSim.L)
    ax.set_zlim(0,newSim.L)
    ab.set_zlim(0,newSim.L)
    
    #ax.text(0,10,0,"Acceptance Rate: " + str(accepted/(accepted+notaccepted)*100) +"%", fontsize=12)
    #ax.text(-10,25,0,"Average Potential Energy: " + str(newSim.AvgPotentialEnergy(newSim.molecules)*epsilon) + " kJ/mol", fontsize=12)
    #ax.text(0,40,0,"Standard Deviation Of Potential Energies: " + str(newSim.stdDevPotentialEnergy(newSim.molecules)*epsilon) + " kJ/mol", fontsize=12)
    plt.show()

initSim()