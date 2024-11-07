from typing import List
from Molecule import Molecule
import numpy as np
import random
import math



class SimulatorUnit:
    
    def __init__(self, molecules: List[Molecule], L: float, T: float):
        self.molecules = molecules
        self.L = L
        self.T=T
        kB = 1.380649*10**-23
        self.beta = 1/(kB*T)

    def MonteCarloStep(self):
        accepted = False
        #select random molecule
        moleculeIndex = random.randint(0, len(self.molecules)-1)
        #store current state of molecules
        moleculeStateCandidate = self.molecules
        #get current potential energy
        oldPotentialEnergy = self.PotentialEnergy(moleculeStateCandidate)
        #randomly move molecule
        movementSize = 100
        moleculeStateCandidate[moleculeIndex].x += random.uniform(-movementSize,movementSize)
        if (moleculeStateCandidate[moleculeIndex].x > self.L):
            moleculeStateCandidate[moleculeIndex].x = moleculeStateCandidate[moleculeIndex].x-self.L
        if (moleculeStateCandidate[moleculeIndex].x < 0):
            moleculeStateCandidate[moleculeIndex].x = moleculeStateCandidate[moleculeIndex].x+self.L
        moleculeStateCandidate[moleculeIndex].y += random.uniform(-movementSize,movementSize)
        if (moleculeStateCandidate[moleculeIndex].y > self.L):
            moleculeStateCandidate[moleculeIndex].y = moleculeStateCandidate[moleculeIndex].y-self.L
        if (moleculeStateCandidate[moleculeIndex].y < 0):
            moleculeStateCandidate[moleculeIndex].y = moleculeStateCandidate[moleculeIndex].y+self.L
        moleculeStateCandidate[moleculeIndex].z += random.uniform(-movementSize,movementSize)
        if (moleculeStateCandidate[moleculeIndex].z > self.L):
            moleculeStateCandidate[moleculeIndex].z = moleculeStateCandidate[moleculeIndex].z-self.L
        if (moleculeStateCandidate[moleculeIndex].z < 0):
            moleculeStateCandidate[moleculeIndex].z = moleculeStateCandidate[moleculeIndex].z+self.L
        #get new potential energy
        self.recalculateMolecule(moleculeIndex)
        newPotentialEnergy = self.PotentialEnergy(moleculeStateCandidate)
        #if new potential energy is lower, accept move
        if newPotentialEnergy < oldPotentialEnergy:
            self.molecules = moleculeStateCandidate
            accepted = True
        #if new potential energy is higher, accept move with probability e^(-beta*(new-old))
        else:
            probablity = min(1, math.exp(-(self.beta)*(newPotentialEnergy-oldPotentialEnergy)))
            if random.uniform(0,1) < probablity:
                self.molecules = moleculeStateCandidate
                accepted = True
        return accepted
        

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
                LJcalc = self.LJ(dist, self.molecules[i].epsilon, self.molecules[i].r0)
                self.molecules[j].LJs[i] = LJcalc
                self.molecules[i].LJs[j] = LJcalc

    ##maths on this is wrong, placeholder
    def populateSolid(self, epsilon, r0, N: int):
        rSigma = 0.77
        aSigma = rSigma*2
        for i in range(math.floor(N/4)):
            x1 = 0 + (aSigma*i)
            y1 = rSigma + (aSigma*i)
            z1 = rSigma + (aSigma*i)
            x2 = rSigma + (aSigma*i)
            y2 = 0 + (aSigma*i)
            z2 = rSigma+ (aSigma*i)
            x3 = 0+ (aSigma*i)
            y3 = 0+ (aSigma*i)
            z3 = 0+ (aSigma*i)
            x4 = rSigma+ (aSigma*i)
            y4 = rSigma+ (aSigma*i)
            z4 = 0+ (aSigma*i)
            self.molecules.append(Molecule(x1, y1, z1, r0, epsilon, [None]*N))
            self.molecules.append(Molecule(x2, y2, z2, r0, epsilon, [None]*N))
            self.molecules.append(Molecule(x3, y3, z3, r0, epsilon, [None]*N))
            self.molecules.append(Molecule(x4, y4, z4, r0, epsilon, [None]*N))
        
        #iterate through every molecule relationship (just combinations, not permutations, so no duplicates)
        for i in range(N):
            for j in range(i+1, N): #for each relationship, calculate LJ potential and store it in both molecules' lists, with index corresponding to the other molecule
                dist=self.distance(self.molecules[i], self.molecules[j])
                LJcalc = self.LJ(dist, self.molecules[i].epsilon, self.molecules[i].r0)
                self.molecules[j].LJs[i] = LJcalc
                self.molecules[i].LJs[j] = LJcalc

