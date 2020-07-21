#!/usr/bin/env python3
import sys
import math
import numpy as np
import akantu as aka


element_type = aka._triangle_6
cohesive_type = aka._cohesive_2d_6
spatial_dimension = 2
max_steps = 1000


class InternalFieldSetter:
    def __init__(self, model, material_id):
        self.model = model
        self.material = model.getMaterial(material_id)
        self.fe_engine = model.getFEEngine()
        element_filter = self.material.getElementFilter(element_type)
        nb_quads = self.fe_engine.getNbIntegrationPoints(element_type) * \
            len(element_filter)

        self.quadrature_points = np.zeros((nb_quads, spatial_dimension))

        self.fe_engine.interpolateOnIntegrationPoints(
            model.getMesh().getNodes(), self.quadrature_points,
            spatial_dimension, element_type, aka._not_ghost,
            element_filter)

    def setField(self, field_id, func):
        _field = self.material.getInternalFieldReal(field_id, element_type)
        nb_quads = len(_field)
        for q in range(nb_quads):
            _field[q, :] = func(self.quadrature_points[q, :])


aka.parseInput("material.dat")
mesh = aka.Mesh(spatial_dimension)
mesh.read("triangle_2.msh")

model = aka.SolidMechanicsModelCohesive(mesh)

rules = {
    ("cathode", "cathode"): "cathode_cathode",
    ("ceramic", "ceramic"): "ceramic_ceramic",
    ("anode", "anode"): "anode_anode",
    ("cathode", "ceramic"): "cathode_ceramic",
    ("anode", "ceramic"): "anode_ceramic"
  }

material_selector = aka.MaterialCohesiveRulesSelector(model, rules)
model.setMaterialSelector(material_selector)

model.initFull(_analysis_method=aka._explicit_lumped_mass,
               _is_extrinsic=True)


cathode_setter = InternalFieldSetter(model, "cathode")
ceramic_setter = InternalFieldSetter(model, "ceramic")
anode_setter = InternalFieldSetter(model, "anode")

fe_engine = model.getFEEngine()
nb_element = model.getMesh().getNbElement(element_type)
nb_quads = fe_engine.getNbIntegrationPoints(element_type) * nb_element
all_quadrature_points = np.zeros((nb_quads, spatial_dimension))
aka.printBacktrace(True)
fe_engine.interpolateOnIntegrationPoints(
    model.getMesh().getNodes(), all_quadrature_points,
    spatial_dimension, element_type)


time_step = model.getStableTimeStep() * 0.05
model.setTimeStep(time_step)
print(f"Time step: {time_step}")

inserter = model.getElementInserter()
inserter.setLimit(aka._y, 0.0, 0.50)
model.updateAutomaticInsertion()

position = mesh.getNodes()
velocity = model.getVelocity()
boundary = model.getBlockedDOFs()
displacement = model.getDisplacement()

nb_nodes = mesh.getNbNodes()

# boundary conditions
for n in range(nb_nodes):
    if position[n, 1] < 0.01:
        boundary[n, 1] = True

    if position[n, 0] < 0.01:
        boundary[n, 0] = True

model.setBaseName("extrinsic")
model.addDumpFieldVector("displacement")
model.addDumpField("velocity")
model.addDumpField("acceleration")
model.addDumpField("internal_force")
model.addDumpField("stress")
model.addDumpField("grad_u")
model.addDumpField("delta_T")
model.dump()

# initial conditions
a = np.array([-1., -1.])
b = np.array([1., 1.])
L = np.linalg.norm(b - a)
direction = (b - a) / L

# Main loop
for s in range(1, max_steps):
    # Here I would like to create a temperature field which linearly changes
    # from point (-1,-1) to point (1,1), the change is then gradually
    # increasing to cause some thermal cracks.

    scale = 0.1 * s
    cathode_setter.setField(
        "delta_T", lambda pos: pos.dot(direction) / L * scale)
    ceramic_setter.setField(
        "delta_T", lambda pos: pos.dot(direction) / L * scale)
    anode_setter.setField(
        "delta_T", lambda pos: pos.dot(direction) / L * scale)
    model.checkCohesiveStress()
    model.solveStep()

    if s % 10 == 0:
        model.dump(time_step * s, s)
        print(f"passing step {s}/{max_steps}")

    Ed = model.getEnergy("dissipated")

    Edt = 200 * math.sqrt(2)
    print(f"{Ed} {Edt}")
    if (Ed < Edt * 0.999) or (Ed > Edt * 1.001) or math.isnan(Ed):
        print("The dissipated energy is incorrect")
        sys.exit(1)
