# NumPy and standard library imports
import numpy as np

# PySPH base and carray imports
from pysph.base.nnps import DomainManager
from pysph.base.utils import get_particle_array_gasd as gpa
from pysph.base.kernels import Gaussian

# PySPH solver and integrator
from pysph.solver.application import Application
from pysph.solver.solver import Solver
from pysph.sph.integrator import PECIntegrator
from pysph.sph.integrator_step import GasDFluidStep

# PySPH sph imports
from pysph.sph.equation import Group
from pysph.sph.gas_dynamics.basic import ScaleSmoothingLength, UpdateSmoothingLengthFromVolume,\
	SummationDensity, IdealGasEOS, MPMAccelerations

# Numerical constants
dim = 1
