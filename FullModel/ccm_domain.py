"""
This file creates the geometrical topology of the 
Chlamy model using FEniCS.
"""

# ==================================================
# ============= A. import packages =================
from dolfin import *
import importlib.util
import numpy as np
import ufl
import os
import sys
set_log_active(False)

# current directory
cdir          = os.path.dirname(os.path.abspath(__file__))
# import ccm_constants
spec_const    = importlib.util.spec_from_file_location("ccm_constants",\
                                                       cdir + "/ccm_constants.py")
ccm_constants = importlib.util.module_from_spec(spec_const)
spec_const.loader.exec_module(ccm_constants)
# ==================================================
# ==================================================


# get outer radius of the thylakoid tubules R_tub, which
# defines how far they extend toward the chloroplast envelope
def get_R_tub(TUBGEO, R_tub_val = None):
    # TUBGEO:    identifier of the tubule geometry
    # R_tub_val: optional user input to specify the R_tub
    if TUBGEO == "extended":
        return ccm_constants.R_chlor
    elif TUBGEO == "limited":
        if R_tub_val == None:
            R_tub_o = ccm_constants.R_pyr
        else:
            R_tub_o = R_tub_val
        return R_tub_o
    else:
        print("Warning (get_R_tub): Unrecognized identifier '%s'\n"%TUBGEO)
        return ccm_constants.R_pyr



def check_inner_boundary(mesho, innbdo, expectedNum):
    count = 0
    for vt in vertices(mesho):
        if innbdo.inside(vt.point(), True):
            count = count + 1
    return (count == expectedNum)



# generating mesh with marked domains
def genMesh(L, N_points, R_tub, tubule_marker, pyr_marker, pyr_mem_marker,
            chlor_mem_marker):

    mesh = IntervalMesh(N_points, 0, L)

    class Tubule(SubDomain):
        def inside(self, x, on_boundary):
            return x[0] <= R_tub + DOLFIN_EPS

    class Pyrenoid(SubDomain):
        def inside(self, x, on_boundary):
            return x[0] <= ccm_constants.R_pyr + DOLFIN_EPS

    class Right(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[0], ccm_constants.R_chlor)

    class PyrenoidMembrane(SubDomain):
        def inside(self, x, on_boundary):
            return (x[0] > ccm_constants.R_pyr - L / N_points / 4) and (
                x[0] < ccm_constants.R_pyr + L / N_points / 4)

    domain_all = MeshFunction('size_t', mesh, mesh.topology().dim())
    domain_all.set_all(0)
    domain_tub = Tubule()
    domain_tub.mark(domain_all, tubule_marker)
    domain_pyr = Pyrenoid()
    domain_pyr.mark(domain_all, pyr_marker)

    boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    boundaries.set_all(0)

    right = Right()
    right.mark(boundaries, chlor_mem_marker)
    pyrmem = PyrenoidMembrane()
    pyrmem.mark(boundaries, pyr_mem_marker)

    if check_inner_boundary(mesh, pyrmem, 1) == 0:
        print(
            "Error in genMesh: number of points on the inner boundary is unexpected!"
        )
        exit(0)

    return mesh, domain_all, boundaries
