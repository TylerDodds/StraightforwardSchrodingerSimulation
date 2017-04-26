from vispy import scene, color as coloring
from SimulationState import GridDimensions, BoundaryConditionType, BoundaryConditionSetter
import numpy as np

class QuantumSimulationLines():
    def __init__(self, gridExtents : GridDimensions):
        self.Real = EditableSimulationLine.InitialLine(1, "blue", gridExtents)
        self.Imaginary = EditableSimulationLine.InitialLine(2, "red", gridExtents)
        self.Hamiltonian = EditableSimulationLine.InitialLine(3, "green", gridExtents)
        self.ModulusSquared = scene.visuals.Line(pos=np.zeros((gridExtents.NumPoints, 3), dtype=np.float32), color='w', width=3,
                                     antialias=True, method='gl')

    def AddToView(self, view: scene.widgets.ViewBox):
        view.add(self.Real)
        view.add(self.Imaginary)
        view.add(self.Hamiltonian)
        view.add(self.ModulusSquared)

class EditableSimulationLine(scene.visuals.Line):
    """
    Mouse editable extension to the Line visual.
    This class adds mouse picking for line points, mouse_move handling for dragging existing points.
    """

    def __init__(self, buttonToUseIndex = -1, gridExtents = GridDimensions(0, 0), edgeColor ="red",
                 boundaryConditionSetter = BoundaryConditionSetter(BoundaryConditionType.Default), *args, **kwargs):
        scene.visuals.Line.__init__(self, *args, **kwargs)

        # initialize point markers
        self.edgeColor = coloring.Color(edgeColor)
        self.markers = scene.visuals.Markers()
        self.marker_colors = np.ones((len(self.pos), 4), dtype=np.float32)
        self.markers.set_data(pos=self.pos, symbol="s", edge_color=edgeColor, size=6)
        self.selected_point = None
        self.selected_index = -1
        self.buttonToUse = buttonToUseIndex
        self.gridExtents = gridExtents
        self.boundaryConditionSetter = boundaryConditionSetter
        self.beforeSimulation = True

    @classmethod
    def InitialLine(cls, buttonToUseIndex, edgeColor, gridExtents, *args, **kwargs):
        return cls(buttonToUseIndex, gridExtents, edgeColor, BoundaryConditionSetter(BoundaryConditionType.SquareWellWavefunction),
                   pos=np.zeros((gridExtents.NumPoints, 3), dtype=np.float32), color=edgeColor, width=3, antialias=True, method='gl')

    def draw(self, transforms):
        # draw line and markers
        scene.visuals.Line.draw(self, transforms)
        self.markers.draw(transforms)

    def select_point_by_x_position_and_mouse_button(self, event):
        """ If the appropriate mouse button is pressed, get the relevant point of this line from the x-position of the mouse """
        if event.button == self.buttonToUse:
            # position in scene/document coordinates
            pos_scene = event.pos[:3]
            return self.get_closest_point_by_x(pos_scene)
        # no point found, return None
        return None, -1

    def get_closest_point_by_x(self, pos):
        """Finds the closest point by x-position, within the point grid spacing"""
        xPos = pos[0]
        index = (xPos - self.gridExtents.BoxMin) * (self.gridExtents.NumPoints - 1) / (self.gridExtents.BoxWidth)
        index = int(round(index))
        if(index >= 0 and index < self.gridExtents.NumPoints):
            return self.pos[index], index
        else:
            return None, index

    def select_point(self, event):
        """Selects a point on this line, if appropriate"""
        return self.select_point_by_x_position_and_mouse_button(event)

    def update_markers(self, selected_index=-1):
        """ update marker colors, and highlight a marker with a given color """
        self.marker_colors.fill(1)#fill with white
        # default shape and size (non-highlighted)
        shape = "o"
        size = 6
        if 0 <= selected_index < len(self.marker_colors):
            self.marker_colors[selected_index] = self.edgeColor.rgba
            # if there is a highlighted marker, change all marker shapes to a square shape and larger size
            shape = "s"
            size = 8
        self.markers.set_data(pos=self.pos, symbol=shape, edge_color=self.edgeColor,
                              size=size, face_color=self.marker_colors)

    def on_mouse_press(self, event):
        if self.beforeSimulation:
            pos_scene = event.pos[:3]

            if event.button == self.buttonToUse:
                # find closest point to mouse and select it
                self.selected_point, self.selected_index = self.select_point(event)
                self.update_markers(self.selected_index)
                self.lastMouseDownPosition = pos_scene

    def on_mouse_release(self, event):
        if self.beforeSimulation:
            #self.print_mouse_event(event, 'Mouse release')
            self.selected_point = None
            self.update_markers()
            self.lastMouseDownPosition = None

    def on_mouse_move(self, event):
        if self.beforeSimulation:
            eventPos = event.pos[:3]
            if event.button == self.buttonToUse:
                if self.lastMouseDownPosition is not None:
                    lastPoint, lastIndex = self.get_closest_point_by_x(self.lastMouseDownPosition)
                    newPoint, newIndex = self.get_closest_point_by_x(eventPos)
                    indexDifference = newIndex - lastIndex
                    if indexDifference == 0:
                        self.try_select_and_update_point(event, eventPos)
                    else:
                        step = -1 if indexDifference < 0 else 1
                        indicesToTry = range(lastIndex, newIndex, step)
                        for index in indicesToTry:
                            if index >= 0 and index < self.gridExtents.NumPoints:
                                fraction = (index - lastIndex) / (newIndex - lastIndex)
                                interpolated = self.lastMouseDownPosition + (eventPos - self.lastMouseDownPosition) * fraction
                                self.selected_point = self.pos[index]
                                self.selected_index = index
                                self.updatePoint(interpolated)
                else:
                    # find closest point to mouse and select it
                    self.try_select_and_update_point(event, eventPos)

                self.lastMouseDownPosition = eventPos

            else:
                self.lastMouseDownPosition = None

    def try_select_and_update_point(self, event, eventPos):
        """Try to select and update a point based on event position"""
        self.selected_point, self.selected_index = self.select_point(event)
        if self.selected_point is not None:
            self.updatePoint(eventPos)

    def updatePoint(self, pos_scene):
        """Update currently-selected point to new position given by scene position"""
        self.selected_point[1] = pos_scene[1]
        self.boundaryConditionSetter.UpdateLineBoundaryConditions(self)
        self.set_data(pos=self.pos)
        self.update_markers(self.selected_index)

    def update_after_simulation(self):
        self.set_data(pos=self.pos)
        self.markers.set_data(pos=self.pos, symbol='o', edge_color=self.edgeColor, size=6, face_color=self.marker_colors)