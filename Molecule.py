from typing import List
import numpy as np

class Molecule:
    def __init__(self, x:float, y:float, z:float, r0:float, epsilon:float, LJs:List[float]):
        self.x = x
        self.y = y
        self.z = z
        self.r0 = r0
        self.epsilon = epsilon
        self.LJs = LJs
    
    