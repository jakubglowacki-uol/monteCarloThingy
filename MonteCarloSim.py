from SimulatorUnit import SimulatorUnit 
from typing import List



def initSim():
    print("Initializing simulation")
    #create simulator cube
    newSim = SimulatorUnit([], 250, 295.15)
    #populate cube with molecules
    newSim.populateSolid(1, 1, 256)

    for n,molecule in enumerate(newSim.molecules):
        print ("Molecule " + str(n) + " has coords: " + str(molecule.x) + " " + str(molecule.y) + " " + str(molecule.z))

    accepted = 0
    notaccepted = 0
    for i in range(1000):
        stepAccepted = newSim.MonteCarloStep()
        if stepAccepted:
            accepted += 1
        else:
            notaccepted += 1
    print("Acceptance rate: " + str(accepted/(accepted+notaccepted)))
    for n,molecule in enumerate(newSim.molecules):
         print ("Molecule " + str(n) + " has coords: " + str(molecule.x) + " " + str(molecule.y) + " " + str(molecule.z))

    # print("distance between molecule 1 and 5: " + str(newSim.distance(newSim.molecules[1], newSim.molecules[5])))

initSim()