import math

import numpy as np
from vispy import scene

from SimulationState import GridDimensions, BoundaryConditionSetter, BoundaryConditionType
from SimulationLines import QuantumSimulationLines, EditableSimulationLine

class SquareWellUpdater():
    """Updates one-dimensional wavefunctions in a finite-sized periodic or square well domain."""
    def __init__(self, simulationLines: QuantumSimulationLines, gridExtents: GridDimensions):
        self.simulationLines = simulationLines
        self.real = simulationLines.Real
        self.imag = simulationLines.Imaginary
        self.hamil = simulationLines.Hamiltonian
        self.mod = simulationLines.ModulusSquared
        self.gridExtents = gridExtents
        self.deltaPosition = gridExtents.PointExtent
        self.deltaPositionSquared = self.deltaPosition * self.deltaPosition
        self.hBarSquared = pow(1.0, 2)
        self.mass = 1
        self.deltaTValue = 0.1
        self.simulationStepsPerFrame = 30
        self.boundaryConditionSetter = BoundaryConditionSetter(BoundaryConditionType.SquareWellWavefunction)

        self.InitializeEditLineVisual(self.real)
        self.InitializeEditLineVisual(self.imag)
        self.InitializeEditLineVisual(self.hamil)
        self.InitializeLine(self.mod)

        for index in range(0, len(self.real.pos)):
            fraction = index / len(self.real.pos)
            self.real.pos[index, 1] = math.sin(math.pi * 2 * fraction) * 20
            self.imag.pos[index, 1] = -self.real.pos[index, 1]

        self.real.update_after_simulation()
        self.imag.update_after_simulation()
        self.set_norm_squared()
        self.updateModAfterSimulation()

    def InitializeEditLineVisual(self, lineVisual: EditableSimulationLine):
        self.InitializeLine(lineVisual)
        lineVisual.update_after_simulation()

    def InitializeLine(self, line: scene.visuals.Line):
        line.pos[:, 0] = np.linspace(self.gridExtents.BoxMin, self.gridExtents.BoxMax, self.gridExtents.NumPoints)

    def ContinueSimulation(self):
        """Continue the simulation"""
        self.real.beforeSimulation = False
        self.imag.beforeSimulation = False
        self.set_norm_squared()

    def PauseSimulation(self):
        """Pause the simulation"""
        self.real.beforeSimulation = True
        self.imag.beforeSimulation = True
        self.set_norm_squared()

    def set_norm_squared(self):
        """Determine and set the rough integral of modulus squared of the wavefunction"""
        realNormSquare = np.sum(np.square(self.real.pos[:,1]))
        imagNormSquare = np.sum(np.square(self.imag.pos[:, 1]))
        sumOfSquares = realNormSquare + imagNormSquare
        roughIntegral = sumOfSquares / self.gridExtents.BoxWidth
        self.integralPsiSquared = roughIntegral

    def update(self, event):
        #So we are discretizing R at integer steps of delta T, and I at half-Integer steps
        #R(t + deltaT / 2) = R(t - deltaT / 2) + deltaT*H*I(t)
        #I(t + deltaT / 2) = I(t - deltaT / 2) - deltaT*H*R(t)
        #at any given integer time step deltaT * i, we will have R at deltaT * i, I at deltaT * (i + 1/2). Call these R_i, I_i
        #this includes initial condition i == 0
        #so  R_i = R_(i-1) + deltaT * H * I_(i-1)
        #and I_i = I_(i-1) - deltaT * H * R_(i)
        #they can be calculated in that order
        #P(x,t)=R^2(x,t)+I(x,t+.5Δt)I(x,t−.5Δt) at integer time steps, but we won't need to find this exactly.

        for i in range(0, self.simulationStepsPerFrame):
            self.updateRealAndImaginaryPosValuesFromHamiltonianStep()

        self.real.update_after_simulation()
        self.imag.update_after_simulation()
        self.updateModAfterSimulation()

    def updateModAfterSimulation(self):
        """Update modulus line after simulation"""
        if self.integralPsiSquared > 0:
            denominator = np.sqrt(self.integralPsiSquared)
            for index, previousModValue in enumerate(self.mod.pos[:, 1]):
                self.mod.pos[index, 1] = self.modulusSquared(self.real.pos[index, 1], self.imag.pos[index, 1]) / denominator
                #this will apply rough rescaling of modulus squared so that mod squared and wavefunction real, imag parts are in same proportion as if we'd actually normalized the whole thing
            self.mod.set_data(pos = self.mod.pos)

    def modulusSquared(self, realValue, imagValue):
        """Modulus squared from real, imaginary value"""
        return realValue * realValue + imagValue * imagValue;

    def updateRealAndImaginaryPosValuesFromHamiltonianStep(self):
        """Perform evolution of real, imaginary parts of wavefuction"""
        secondDerivativeImaginaryPart = self.second_derivative_periodic(self.imag.pos)
        potentialFraction = 0.0001#rescaling fraction of Potential for stability, since we've not assumed realistic values for mass, hBar
        self.real.pos[:, 1] += (self.deltaTValue) * (- secondDerivativeImaginaryPart * (0.5 * self.hBarSquared / self.mass) + potentialFraction * self.hamil.pos[:, 1] * self.imag.pos[:, 1])
        self.boundaryConditionSetter.UpdateLineBoundaryConditions(self.real)
        secondDerivativeRealPart = self.second_derivative_periodic(self.real.pos)
        self.imag.pos[:, 1] -= (self.deltaTValue) * (- secondDerivativeRealPart * (0.5 * self.hBarSquared / self.mass) + potentialFraction * self.hamil.pos[:, 1] * self.real.pos[:, 1])
        self.boundaryConditionSetter.UpdateLineBoundaryConditions(self.imag)

    def second_derivative_periodic(self, pos):
        """Second derivative, periodic"""
        posLen = len(pos)
        posVal = pos[:,1]
        deriv = np.array([0.] * posLen)
        for index, positionValue in enumerate(posVal):
            if index == 0:
                #deriv[index] = (posVal[index + 2] + positionValue - posVal[index + 1] * 2.) / (self.DeltaPositionSquared);
                deriv[index] = (posVal[index + 1] + posVal[index - 1 + posLen] - positionValue * 2.) / (self.deltaPositionSquared);
            elif index == posLen - 1:
                #deriv[index] = (posVal[index - 2] + positionValue - posVal[index - 1] * 2.) / (self.DeltaPositionSquared);
                deriv[index] = (posVal[index + 1 - posLen] + posVal[index - 1] - positionValue * 2.) / (self.deltaPositionSquared);
            else:
                deriv[index] = (posVal[index + 1] + posVal[index - 1] - positionValue * 2.) / (self.deltaPositionSquared);
        return deriv

    def second_derivative_clamped(self, pos):
        """Second derivative, clamped at domain edges"""
        posLen = len(pos)
        posVal = pos[:,1]
        deriv = np.array([0.] * posLen)
        for index, positionValue in enumerate(posVal):
            if index == 0:
                deriv[index] = (posVal[index + 2] + positionValue - posVal[index + 1] * 2.) / (self.deltaPositionSquared);
            elif index == posLen - 1:
                deriv[index] = (posVal[index - 2] + positionValue - posVal[index - 1] * 2.) / (self.deltaPositionSquared);
            else:
                deriv[index] = (posVal[index + 1] + posVal[index - 1] - positionValue * 2.) / (self.deltaPositionSquared);
        return deriv