from typing import List
from Molecule import Molecule
import numpy as np
import random
import math
import statistics
import rdfpy


class SimulatorUnit:
    
    def __init__(self, molecules: List[Molecule], L: float, T: float):
        self.molecules = molecules
        self.L = L
        self.T=T
        self.kB = 1.38064852e-23/128000 #boltzmann constant in epsilon units
        self.potEnergies = []

    def getMoleculeCoors(self):
        return [[molecule.x, molecule.y, molecule.z] for molecule in self.molecules]
    
    def getMolecules(self):
        return self.molecules

    def MonteCarloStep(self):
        accepted = False
        #select random molecule
        moleculeIndex = random.randint(0, len(self.molecules)-1)
        #store current state of molecules
        moleculeStateCandidate = self.molecules
        #get current potential energy
        oldPotentialEnergy = self.PotentialEnergy(moleculeStateCandidate)
        self.potEnergies.append(oldPotentialEnergy)
        #randomly move molecule
        movementSize = self.T/25
        moleculeStateCandidate[moleculeIndex].x += random.uniform(-movementSize,movementSize)
        if (moleculeStateCandidate[moleculeIndex].x > self.L):
            moleculeStateCandidate[moleculeIndex].x = moleculeStateCandidate[moleculeIndex].x-self.L
        if (moleculeStateCandidate[moleculeIndex].x < 0):
            moleculeStateCandidate[moleculeIndex].x = self.L + moleculeStateCandidate[moleculeIndex].x
        moleculeStateCandidate[moleculeIndex].y += random.uniform(-movementSize,movementSize)
        if (moleculeStateCandidate[moleculeIndex].y > self.L):
            moleculeStateCandidate[moleculeIndex].y = moleculeStateCandidate[moleculeIndex].y-self.L
        if (moleculeStateCandidate[moleculeIndex].y < 0):
            moleculeStateCandidate[moleculeIndex].y = self.L + moleculeStateCandidate[moleculeIndex].y
        moleculeStateCandidate[moleculeIndex].z += random.uniform(-movementSize,movementSize)
        if (moleculeStateCandidate[moleculeIndex].z > self.L):
            moleculeStateCandidate[moleculeIndex].z = moleculeStateCandidate[moleculeIndex].z-self.L
        if (moleculeStateCandidate[moleculeIndex].z < 0):
            moleculeStateCandidate[moleculeIndex].z = self.L + moleculeStateCandidate[moleculeIndex].z
        #get new potential energy
        self.recalculateMolecule(moleculeIndex, moleculeStateCandidate)
        newPotentialEnergy = self.PotentialEnergy(moleculeStateCandidate)
        #if new potential energy is lower, accept move
        if newPotentialEnergy <= oldPotentialEnergy:
            self.molecules = moleculeStateCandidate
            accepted = True
        #if new potential energy is higher, accept move with probability e^(-beta*(new-old))
        else:
            
            probablity = min(1, math.exp(-((newPotentialEnergy-oldPotentialEnergy)/self.kB*self.T)))
            print(-((newPotentialEnergy-oldPotentialEnergy)/self.kB*self.T))
            if random.uniform(0,1) < probablity:
                self.molecules = moleculeStateCandidate
                accepted = True
        return accepted
        

    ##recalculate molecule relationships after a move
    def recalculateMolecule(self, moleculeIndex: int, moleculeList: List[Molecule]):
        for i in range(len(moleculeList)):
            if i != moleculeIndex:
                dist=self.distance(moleculeList[moleculeIndex], moleculeList[i])
                LJcalc = self.LJ(dist)
                moleculeList[moleculeIndex].LJs[i] = LJcalc
                moleculeList[i].LJs[moleculeIndex] = LJcalc
            else:
                moleculeList[moleculeIndex].LJs[i] = 0

    def AvgPotentialEnergy(self, moleculeList: List[Molecule]):
        return self.PotentialEnergy(moleculeList)/len(moleculeList)
    
    def stdDevPotentialEnergy(self, moleculeList: List[Molecule]):
        PEs =[]
        for molecule in moleculeList:
            moleculePE = 0
            for LJval in molecule.LJs:
                if LJval != None:
                    moleculePE += LJval
            PEs.append(moleculePE)
        return statistics.stdev(PEs)
    
    #Calculating potential energy of system
    def PotentialEnergy(self, moleculeList: List[Molecule]):
        potentialEnergy = 0
        for i in range(len(moleculeList)):
            for j in range(i+1, len(moleculeList)):
                potentialEnergy += moleculeList[i].LJs[j]
        return potentialEnergy

    #calculating LJ potential between two molecules based on distance
    def LJ(self, r: float):
        return 4*(self.RepulsionComponent(r)-self.AttractionComponent(r))

    def AttractionComponent(self, r: float):
        return (1/r)**6

    def RepulsionComponent(self, r:float):
        return (1/r)**12

    def distance(self, molecule1: Molecule, molecule2: Molecule): #minimum image convention
        x1 = molecule1.x
        y1 = molecule1.y
        z1 = molecule1.z
        x2 = molecule2.x
        y2 = molecule2.y
        z2 = molecule2.z
        L = self.L
        delta_x=0
        delta_y=0
        delta_z=0
        # Periodic boundary conditions for x
        if (-0.5)*L<(x1-x2)<0.5*L:
            delta_x=x1-x2
        if (x1-x2)<-0.5*L:
            delta_x=x1-x2+L
        if (x1-x2)>0.5*L:
            delta_x=x1-x2-L
        # Periodic boundary conditions for y
        if (-0.5)*L<(y1-y2)<0.5*L:
            delta_y=y1-y2
        if (y1-y2)<-0.5*L:
            delta_y=y1-y2+L
        if (y1-y2)>0.5*L:
            delta_y=y1-y2-L
        # Periodic boundary conditions for z
        if (-0.5)*L<(z1-z2)<0.5*L:
            delta_z=z1-z2
        if (z1-z2)<-0.5*L:
            delta_z=z1-z2+L
        if (z1-z2)>0.5*L:
            delta_z=z1-z2-L
        dis=math.sqrt((delta_x)**2 + (delta_y)**2 + (delta_z)**2)
        return dis

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
                LJcalc = self.LJ(dist)
                self.molecules[j].LJs[i] = LJcalc
                self.molecules[i].LJs[j] = LJcalc

    ##maths on this is wrong, placeholder
    def populateSolid(self, epsilon, r0, N, offset=0):
        r_sig = 0.77
        a_sig = r_sig*2
        iterNo = math.floor(N/4)
        cube_root = math.ceil(iterNo**(1./3))
        for x in range(cube_root):
            for y in range(cube_root):
                for z in range(cube_root):
                    x1 = 0 + (a_sig*x) + offset
                    y1 = r_sig + (a_sig*y)+ offset
                    z1 = r_sig + (a_sig*z)+ offset

                    x2 = r_sig + (a_sig*x)+ offset
                    y2 = 0 + (a_sig*y)+ offset
                    z2 = r_sig+ (a_sig*z)+ offset

                    x3 = 0+ (a_sig*x)+ offset
                    y3 = 0+ (a_sig*y)+ offset
                    z3 = 0+ (a_sig*z)+ offset

                    x4 = r_sig+ (a_sig*x)+ offset
                    y4 = r_sig+ (a_sig*y)+ offset
                    z4 = 0+ (a_sig*z)+ offset

                    self.molecules.append(Molecule(x1, y1, z1, r0, epsilon, [None]*N))
                    self.molecules.append(Molecule(x2, y2, z2, r0, epsilon, [None]*N))
                    self.molecules.append(Molecule(x3, y3, z3, r0, epsilon, [None]*N))
                    self.molecules.append(Molecule(x4, y4, z4, r0, epsilon, [None]*N))
        
        #iterate through every molecule relationship (just combinations, not permutations, so no duplicates)
        for i in range(N):
            for j in range(i+1, N): #for each relationship, calculate LJ potential and store it in both molecules' lists, with index corresponding to the other molecule
                dist=self.distance(self.molecules[i], self.molecules[j])
                LJcalc = self.LJ(dist)
                self.molecules[j].LJs[i] = LJcalc
                self.molecules[i].LJs[j] = LJcalc

