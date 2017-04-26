from enum import Enum

class BoundaryConditionType(Enum):
    """Type of Boundary Condition for Wavefunction or Hamiltonian"""
    Default = 0
    SquareWellWavefunction = 1
    SquareWellHamiltonian = 2
    PeriodicWavefunction = 3
    PeriodicHamiltonian = 4


class BoundaryConditionSetter():
    """Sets points on boundaries of Wavefunction or Hamiltonian"""

    def __init__(self, bcType:BoundaryConditionType):
        self.BoundaryConditionType = bcType

    def UpdateLineBoundaryConditions(self, line):
        if self.BoundaryConditionType is BoundaryConditionType.SquareWellWavefunction:
            line.pos[0, 1] = 0.
            line.pos[-1, 1] = 0.
        elif self.BoundaryConditionType is BoundaryConditionType.PeriodicWavefunction:
            #do nothing, since periodic second derivative is fine
            #Note that for a range of [a:b] inclusive, a and b will be physically distinct points,
            # separated (like all other points) by deltaX distance
            #Since they are distinct, we do not need to set: line.pos[-1, 1] = line.pos[0, 1]
            pass

class GridDimensions():
    """Dimensions of simulation grid"""
    def __init__(self, boxExtent: float, numPoints: int):
        self.BoxMin = -boxExtent
        self.BoxMax = boxExtent
        self.NumPoints = numPoints
        self.PointExtent = 0.5 * (self.BoxMax - self.BoxMin) / (numPoints - 1)
        self.BoxWidth = self.BoxMax - self.BoxMin