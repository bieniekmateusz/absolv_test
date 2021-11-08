from simtk.openmm import unit
from simtk import openmm
from simtk.openmm import app
import numpy

from absolv.models import EquilibriumProtocol, State, System, TransferFreeEnergySchema

schema = TransferFreeEnergySchema(
    # Define the solutes in the system. There may be multiple in the case of,
    # e.g., ion pairs like Na+ + Cl-. Here we use `None` to specify that the solute 
    # will be transferred into a vacuum.
    system=System(solutes={"[H]c1c(c(c2c(c1[H])C(=C(N2[H])C3(C(C(C(C3([H])[H])([H])[H])([H])[H])([H])[H])[H])[H])[H])[H]": 1}, solvent_a={"O": 895}, solvent_b=None),
    # Define the state that the calculation will be performed at.
    state=State(temperature=298 * unit.kelvin, pressure=1.0 * unit.atmosphere),
    # Define the alchemical pathway to transform the solute along in the first
    # and second (i.e. vacuum) solvent respectively.
    alchemical_protocol_a=EquilibriumProtocol(
        lambda_sterics=[  # fmt: off
            1.00, 1.00, 1.00, 1.00, 1.00, 0.95, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40,
            0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00,
        ],
        lambda_electrostatics=[  # fmt: off
            1.00, 0.75, 0.50, 0.25, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
            0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
        ],
    ),
    alchemical_protocol_b=EquilibriumProtocol(
        lambda_sterics=[  # fmt: off
            1.00, 1.00, 1.00, 1.00, 1.00, 0.95, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40,
            0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00,
        ],
        lambda_electrostatics=[  # fmt: off
            1.00, 0.75, 0.50, 0.25, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
            0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
        ],
    ),
)

forcefield = app.ForceField('mol.xml', '../jm_gmx_tip3p.xml')

def myff(off_top, crd):
    # take the off topology and crd and return the omm system
    omm_top = off_top.to_openmm()
    # was needed but should be determined based on the number of waters?
    omm_top.setPeriodicBoxVectors(unit.Quantity(numpy.diag([34, 34, 34]), unit.angstrom))
    # import pdb ; pdb.set_trace()
    system = forcefield.createSystem(
        omm_top,
        nonbondedMethod=openmm.app.PME,
        nonbondedCutoff=1.0*unit.nanometer,
        constraints=openmm.app.HBonds
    )
    return system

from absolv.runners.equilibrium import EquilibriumRunner
EquilibriumRunner.setup(schema, myff)

EquilibriumRunner.run(schema, platform="CUDA")

free_energies = EquilibriumRunner.analyze(schema)
print(free_energies)
