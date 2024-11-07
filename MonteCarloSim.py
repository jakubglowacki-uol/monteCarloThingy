import matplotlib.pyplot as plt
from SimulatorUnit import SimulatorUnit 
from typing import List



def initSim():
    print("Initializing simulation")
    #create simulator cube
    newSim = SimulatorUnit([], 1000)
    #populate cube with molecules
    newSim.populateRandom(128, 342, 20)

    for n,molecule in enumerate(newSim.molecules):
        print ("Molecule " + str(n) + " has coords: " + str(molecule.x) + " " + str(molecule.y) + " " + str(molecule.z))

    print("molecules 0 to 4 have LJ potential: " + str(newSim.molecules[0].LJs[4]) + ", " + str(newSim.molecules[4].LJs[0]))

initSim()