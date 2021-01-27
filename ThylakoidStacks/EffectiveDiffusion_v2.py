# =============================================== #
# ============ Import head packages  ============ #
from fenics import *
import pygmsh as pg
import os
import meshio
import numpy as np
import csv
from EffectiveDiffusion_mesh import createMeshFull
# =============================================== #
# =============================================== #




# =============================================== #
# =========== Set geometry parameters =========== #
N_STACK = 10                                       # Number of thylakoid stacks in the simulation domain
f_nrrw  = 0.4                                      # Fraction of thylakoid stacks with narrow gaps
N_nrrw  = int(np.ceil(N_STACK * f_nrrw))
d_t     = 5/800                                    # thickness of the thylkaoid membranes
d_h     = 10/800                                   # height of the thylakoid lumen
d_s     = 4/800                                    # small spacing between adjacent thylakoid stacks
d_l     = 1/N_STACK - 4*d_t - 2*d_h - d_s          # large spacing between adjacent thylakoid stacks

Delta_n_array = [10 / 800]                         # list of the widths of narrow gaps
diffDelta     = (50.2 - 5.6)/800                   # width difference between the narrow and wide gaps
# Delta_n can be varied hypothetically, and we assume that
# diffDelta = Delta_w - Delta_n is fixed.

# mesh_size
lmesh1  = 0.002
lmesh2  = 0.002

D_0       = 1.0                                    # free diffusion coefficient
D_m_array = pow(10,np.linspace(-3.0,-0.3,28))      # hypothetical diffusion coefficent across the membrane
# =============================================== #
# =============================================== #




# =============================================== #
# =========== Boundaries and domains ============ #
             
def DiffusionConstant(mesh, domain_all, D_0, D_1):
    V0   = FunctionSpace(mesh,'DG',0)
    Dc   = Function(V0)
    for i in np.arange(len(domain_all)):
        if domain_all[i]%2 == 0:
            Dc.vector()[i] = D_1
        else:
            Dc.vector()[i] = D_0
    return Dc
                
def bottom_boundary(x, on_boundary):
    return on_boundary and near(x[1],0.);
    
def top_boundary(x, on_boundary):
    return on_boundary and near(x[1],1.);
    
class Top(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 1.)

class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 0.)


# =============================================== #
# =============================================== #




# =============================================== #
# ============== I/O directories ================ #

maindir = "data_dt=" + "%.2E_"%d_t\
             + "dh=" + "%.2E_"%d_h\
             + "N="  + "%d_"%N_STACK\
             + "f="  + "%.2E"%f_nrrw + "/"

if not os.path.exists(maindir):
    os.makedirs(maindir)

for Delta_n in Delta_n_array:
    if not os.path.exists(maindir + 'Delta_n=%.2E'%Delta_n + '/'):
        os.makedirs(maindir + 'Delta_n=%.2E'%Delta_n + '/')

print(maindir)
# =============================================== #
# =============================================== #




# =============================================== #
# ================= Main loop  ================== #

for Delta_n in Delta_n_array:
    
    Delta_dir = maindir + 'Delta_n=%.2E'%Delta_n + '/'
    geoname   = Delta_dir + "test.geo"
    mshname   = Delta_dir + "test.msh"
    xmlname   = Delta_dir + "test.xml"
    
    # create mesh
    createMeshFull(N_STACK, N_nrrw, d_t, d_h, d_s, d_l, Delta_n, Delta_n + diffDelta, lmesh1, lmesh2, geoname, mshname, xmlname)
    mesh = Mesh(xmlname)
    
    
    # function spaces
    V    = FunctionSpace(mesh,'P',1)
    V0   = FunctionSpace(mesh,'DG',0)
    domain_all = MeshFunction("size_t", mesh, Delta_dir +"test_physical_region.xml")
    
    # bcs and boundary markers
    bcs = [DirichletBC(V, Constant(1.0), top_boundary),\
    DirichletBC(V, Constant(0.0), bottom_boundary)]
    
    boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    boundaries.set_all(0)
    bottom =  Bottom()
    bottom.mark(boundaries, 1)
    top    = Top()
    top.mark(boundaries, 2)
    ds = Measure("ds", domain=mesh, subdomain_data=boundaries)
    
    # loop over different Dm values
    Deff = np.zeros(len(D_m_array))
    for loop in np.arange(len(D_m_array)):
        
        D_m = D_m_array[loop]
        D_c = DiffusionConstant(mesh, domain_all, D_0, D_m)
        
        if loop == 0:
            dfile = File( Delta_dir + "D.pvd")
            dfile << project(D_c,V0)
        
        u    = Function(V)
        u_ = TestFunction(V)
        du = TrialFunction(V)
        
        # compute concentration profile
        L = project(D_c,V0)*inner(grad(u),grad(u_))*dx
        solve(L==0, u, bcs=bcs)
        
        # compute effective diffusion
        J1         = assemble(grad(u)[1]*ds(2))
        J2         = assemble(grad(u)[1]*ds(1))
        Deff[loop] = 0.5*(J1+J2)
        
    # save Deff
    Deff_file = open(Delta_dir + "Deff_value.txt", "w")
    for loop in np.arange(len(D_m_array)):
        Deff_file.write('%f'%D_m_array[loop] + '\t%f\n'%Deff[loop])
    Deff_file.close()
# =============================================== #
# =============================================== #



