from vispy import app, scene, color as coloring

from SimulationLines import QuantumSimulationLines
from FiniteDomainSimulation import SquareWellUpdater
from SimulationState import GridDimensions, BoundaryConditionType


class Canvas(scene.SceneCanvas):
    """ A simple test canvas for visualizing the 1-D evolution """
    def __init__(self):
        """Initialize canvas and simulation"""
        scene.SceneCanvas.__init__(self, keys='interactive', size=(1000, 800))

        #set cameras and views
        self.view = self.central_widget.add_view()
        self.view.camera = scene.PanZoomCamera(rect=(-100, -100, 200, 200), aspect=1.0)
        # disable mouse move and wheel events for the camera
        self.view.camera._viewbox.events.mouse_move.disconnect(self.view.camera.viewbox_mouse_event)
        self.view.camera._viewbox.events.mouse_wheel.disconnect(self.view.camera.viewbox_mouse_event)
        scene.visuals.GridLines(parent=self.view.scene)

        # create new editable lines with appropriate extents
        gridExtents = GridDimensions(boxExtent=50, numPoints=51)
        simulationLines = QuantumSimulationLines(gridExtents)
        self.stateHandler = SimpleStateHandler(simulationLines, gridExtents)
        simulationLines.AddToView(self.view)

        self.instructionsLabel = scene.widgets.Console(text_color='white', font_size=10, bgcolor=coloring.Color('black', 0.5))
        self.instructionsLabel.parent = self.central_widget

        self.update_label_text()

        self.timer = app.Timer('auto', connect=self.on_timer, start=True)

        self.show()

    def update_label_text(self):
        """Update instructional label text"""
        self.instructionsLabel.write("Click and drag mouse button to draw:")
        self.instructionsLabel.write("[Left]: Real part (blue). [Right]: Imaginary part (red). [Middle]: Potential (green).")
        if self.stateHandler.running:
            self.instructionsLabel.write("[Space]: Stop simulation.")
        else:
            self.instructionsLabel.write("[Space]: Play simulation.")
        if self.stateHandler.squareWellUpdater.boundaryConditionSetter.BoundaryConditionType is BoundaryConditionType.SquareWellWavefunction:
            self.instructionsLabel.write("[P]: Switch to periodic boundary conditions. * Currently: square well.")
        else:
            self.instructionsLabel.write("[W]: Switch to square well boundary conditions. * Currently: periodic.")

    def on_timer(self, event):
        """Update simulation and canvas on timer tick"""
        self.stateHandler.update(event)
        self.update()

    def on_draw(self, event):
        """Draws canvas"""
        self.instructionsLabel.size = (self.size[0], 65)
        self.instructionsLabel.pos = (0, self.size[1] - self.instructionsLabel.size[1])
        scene.SceneCanvas.on_draw(self, event)

    def on_mouse_move(self, event):
        """Handle mouse move"""
        self.stateHandler.on_mouse_move(event)

    def on_key_press(self, event):
        """Handle key press"""
        self.stateHandler.on_key_press(event)
        self.update_label_text()





class SimpleStateHandler():
    """Simple state handler for (un)pausing and switching boundary conditions."""
    def __init__(self, simulationLines : QuantumSimulationLines, gridExtents: GridDimensions):
        self.squareWellUpdater = SquareWellUpdater(simulationLines, gridExtents)
        self.running = False
        self.SetSquareWellBoundaryConditions()

    def update(self, event):
        """Update, if running"""
        if self.running:
            self.squareWellUpdater.update(event)

    def on_key_press(self, event):
        """Handle key press for (un)pausing or switching boundary conditions"""
        if self.running:
            if event.key.name is 'Space':
                self.running = False
        else:
            if event.key.name is 'Space':
                self.running = True

        if event.key.name is 'W':
            self.SetSquareWellBoundaryConditions()
            print("Square well")
        elif event.key.name is 'P':
            self.SetPeriodicBoundaryConditions()
            print("Periodic")

    def SetPeriodicBoundaryConditions(self):
        """Set periodic boundary conditions"""
        self.squareWellUpdater.real.boundaryConditionSetter.BoundaryConditionType = BoundaryConditionType.PeriodicWavefunction;
        self.squareWellUpdater.imag.boundaryConditionSetter.BoundaryConditionType = BoundaryConditionType.PeriodicWavefunction;
        self.squareWellUpdater.hamil.boundaryConditionSetter.BoundaryConditionType = BoundaryConditionType.PeriodicHamiltonian;
        self.squareWellUpdater.boundaryConditionSetter.BoundaryConditionType = BoundaryConditionType.PeriodicWavefunction;

    def SetSquareWellBoundaryConditions(self):
        """"Set square well boundary conditions"""
        self.squareWellUpdater.real.boundaryConditionSetter.BoundaryConditionType = BoundaryConditionType.SquareWellWavefunction;
        self.squareWellUpdater.imag.boundaryConditionSetter.BoundaryConditionType = BoundaryConditionType.SquareWellWavefunction;
        self.squareWellUpdater.hamil.boundaryConditionSetter.BoundaryConditionType = BoundaryConditionType.SquareWellHamiltonian;
        self.squareWellUpdater.boundaryConditionSetter.BoundaryConditionType = BoundaryConditionType.SquareWellWavefunction;

    def on_mouse_move(self, event):
        """Handle mouse moved"""
        if not self.running:
            if event.button is not None:
                self.squareWellUpdater.set_norm_squared()
                self.squareWellUpdater.updateModAfterSimulation()

    @property
    def running(self):
        """If the simulation is currently running and unpaused"""
        return self.__running

    @running.setter
    def running(self, value):
        """If the simulation is currently running and unpaused"""
        if isinstance(value, bool):
            self.__running = value
            if value:
                self.squareWellUpdater.ContinueSimulation()
            else:
                self.squareWellUpdater.PauseSimulation()




if __name__ == '__main__':
    win = Canvas()
    app.run()