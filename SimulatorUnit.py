from typing import List
from Molecule import Molecule
import numpy as np
import random

class SimulatorUnit:
    def __init__(self, molecules: List[Molecule], L: float):
        self.molecules = molecules
        self.L = L


    def MonteCarloStep(self):
        #select random molecule
        moleculeIndex = random.randint(0, len(self.molecules)-1)
        #store current state of molecules
        moleculeStateCandidate = self.molecules
        #get current potential energy
        oldPotentialEnergy = self.PotentialEnergy()
        #randomly move molecule
        ##TODO: Decide movement range
        moleculeStateCandidate[moleculeIndex].x += random.uniform(-1,1)
        moleculeStateCandidate[moleculeIndex].y += random.uniform(-1,1)
        moleculeStateCandidate[moleculeIndex].z += random.uniform(-1,1)
        #get new potential energy
        self.recalculateMolecule(moleculeIndex)
        newPotentialEnergy = self.PotentialEnergy()
        #if new potential energy is lower, accept move
        if newPotentialEnergy < oldPotentialEnergy:
            self.molecules = moleculeStateCandidate
        #if new potential energy is higher, accept move with probability e^(-beta*(new-old))
        else:
            pass
        pass

    ##recalculate molecule relationships after a move
    def recalculateMolecule(self, moleculeIndex: int):
        for i in range(len(self.molecules)):
            if i != moleculeIndex:
                dist=self.distance(self.molecules[i], self.molecules[moleculeIndex])
                LJcalc = self.LJ(dist, self.molecules[i].epsilon, self.molecules[i].r0)
                self.molecules[moleculeIndex].LJs[i] = LJcalc
                self.molecules[i].LJs[moleculeIndex] = LJcalc
        


    #Calculating potential energy of system
    def PotentialEnergy(self, moleculeList: List[Molecule]):
        potentialEnergy = 0
        for i in range(len(moleculeList)):
            for j in range(i+1, len(moleculeList)):
                potentialEnergy += moleculeList[i].LJs[j]
        return potentialEnergy

    #calculating LJ potential between two molecules based on distance
    def LJ(self, r: float, epsilon:float=1, r0:float=1):
        attraction = self.AttractionComponent(r, r0)
        return 4*epsilon*(self.RepulsionComponent(attraction)-attraction)

    def AttractionComponent(self, r: float, r0:float):
        return (r0/r)**6

    def RepulsionComponent(self, attractionComponent:float):
        return attractionComponent**2

    def distance(self, molecule1: Molecule, molecule2: Molecule): ##update to use minimum image convention
        return np.linalg.norm(np.array([molecule1.x, molecule1.y, molecule1.z])-np.array([molecule2.x, molecule2.y, molecule2.z]))

    #Add a bunch of randomly placed molecules
    def populateRandom(self, epsilon, r0, N):
        for i in range(N):
            x = random.uniform(0, self.L)
            y = random.uniform(0, self.L)
            z = random.uniform(0, self.L)
            self.molecules.append(Molecule(x, y, z, r0, epsilon, [None]*N))

        #iterate through every molecule relationship (just combinations, not permutations, so no duplicates)
        for i in range(N):
            for j in range(i+1, N): #for each relationship, calculate LJ potential and store it in both molecules' lists, with index corresponding to the other molecule
                dist=self.distance(self.molecules[i], self.molecules[j])
                LJcalc = self.LJ(dist, self.molecules[i].epsilon, self.molecules[i].r0)
                self.molecules[j].LJs[i] = LJcalc
                self.molecules[i].LJs[j] = LJcalc

    ##maths on this is wrong, placeholder
    def populateSolid(self, epsilon, r0, N):
        for i in range(N):
            x = i*self.L/N
            y = 0
            z = 0
            self.molecules.append(Molecule(x, y, z, r0, epsilon, []))

