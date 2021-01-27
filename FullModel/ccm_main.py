'''
This script is the main file for the model.
Change model parameters in ccmParaDict, this version
outputs the concentration profile across the chloroplast
with respect to radius.

For more output type options, use ccm_scan.py
----------------------------------------------------------
'''
# ==================================================
# ============= A. import packages =================
from dolfin import *
import numpy as np
from scipy import interpolate as itp
import ufl
import os
import sys
import csv
set_log_active(True)

cdir = os.path.dirname(os.path.abspath(__file__))

# + Add in-house package
import importlib.util
spec_const = importlib.util.spec_from_file_location("ccm_constants",
                                                    cdir + "/ccm_constants.py")
ccm_constants = importlib.util.module_from_spec(spec_const)
spec_const.loader.exec_module(ccm_constants)

spec_domain = importlib.util.spec_from_file_location("ccm_domain",
                                                     cdir + "/ccm_domain.py")
ccm_domain = importlib.util.module_from_spec(spec_domain)
spec_domain.loader.exec_module(ccm_domain)

spec_eqn = importlib.util.spec_from_file_location("ccm_eqn",
                                                  cdir + "/ccm_eqn.py")
ccm_eqn = importlib.util.module_from_spec(spec_eqn)
spec_eqn.loader.exec_module(ccm_eqn)

# # + parameter setting
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["quadrature_degree"] = 4
parameters["form_compiler"]["representation"] = "uflacs"
# ==================================================
# ==================================================


def main():
    return


def getFolderName(problemDict):
    """
    problemDict is a dictionary with all the relevant parameters
    returns a folder name that is specific to each problem.
    """
    namefolder = ""
    geofolder = ""

    for key in problemDict.keys():
        if "_s" in key:
            continue
        if "_e" in key:
            continue
        if "Geo" in key:
            geofolder += key
            if type(problemDict[key]) == str:
                geofolder += "_" + problemDict[key] + "_"
            else:
                geofolder += ("_%.2E_" % problemDict[key])
        else:
            namefolder += key
            if type(problemDict[key]) == str:
                namefolder += "_" + problemDict[key] + "_"
            else:
                namefolder += ("_%.2E_" % problemDict[key])

    return geofolder[:-1] + "/" + namefolder[:-1] + "/"


def makeDir(namefolder):

    if os.path.exists(namefolder):
        os.stat(namefolder)
    else:
        os.makedirs(namefolder)


def makePath(simDict, problemDict):
    """
    It creates the name of the output folder and the path of the
    output file.
    """
    lookup = "version"
    if not (lookup in simDict):
        print("Error specifying version!")
        exit()
    pwd = os.getcwd()
    namefolder = pwd + "/data_" + simDict["version"] + "/" + getFolderName(
        problemDict)
    makeDir(namefolder)

    return namefolder


def saveTOcsv(dict_name, file_name):
    fields = []
    vals = []
    wfile = csv.writer(open(file_name, "w"))
    for key in dict_name:
        fields.append(key)
        vals.append(dict_name[key])
    wfile.writerow(fields)
    wfile.writerow(vals)


def interpDeff(Delta, Permeability, D0):
    Dm_sim = []  # initialize list
    Deff_sim = []

    f = open('Deff_calib/Deff_d=10_Delta=' + str(Delta) + '.txt', 'r')
    for line in f:
        Dm_this, Deff_this = line.split()
        Dm_sim.append(float(Dm_this))
        Deff_sim.append(float(Deff_this))
    f.close()
    fDeff = itp.interp1d(Dm_sim, Deff_sim, fill_value="extrapolate")
    Dm = Permeability * 0.005 / D0
    Deff = fDeff(Dm) * D0

    return Deff


# def get_range(localization, r_start=None, r_end=None):
#     default_geometry = {
#         "diffuse": (0.0, ccm_constants.R_chlor),
#         "stroma": (ccm_constants.R_pyr, ccm_constants.R_chlor),
#         "pyrenoid": (0.0, ccm_constants.R_pyr),
#         "ring": (ccm_constants.R_pyr, ccm_constants.R_pyr + 0.1),
#         "ring2": (ccm_constants.R_pyr, ccm_constants.R_pyr + 0.1)
#     }
#     if (r_start == None) or (r_end == None):
#         r_start, r_end = default_geometry.get(localization, (None, None))
#     return r_start, r_end


def get_permeability(pyrmem_barrier, kcPyrMem=None, khPyrMem=None):
    default_permeability = {"on": (1e-7, 1e-7), "off": (1e1, 1e1)}
    if (kcPyrMem == None) or (khPyrMem == None):
        kcPyrMem, khPyrMem = default_permeability.get(pyrmem_barrier,
                                                      (None, None))
    return kcPyrMem, khPyrMem


def get_rescale_factor(dx, f_geo, r_s_0, r_e_0, r_s_n, r_e_n):
    r = Expression("x[0]", degree=1)
    fv_0 = conditional(le(r, r_e_0),
                       conditional(ge(r, r_s_0), f_geo, Constant(0.0)),
                       Constant(0.0))
    fv_n = conditional(le(r, r_e_n),
                       conditional(ge(r, r_s_n), f_geo, Constant(0.0)),
                       Constant(0.0))

    vol_0 = assemble(r**2 * fv_0 * dx)
    vol_n = assemble(r**2 * fv_n * dx)

    return vol_0 / vol_n


if __name__ == '__main__':

    TubuleMarker = 1
    PyrMarker = 2
    PyrMemMarker = 3
    ChlorMemMarker = 4
    Nmesh = 5000

    ccmSimDict = {
        "version": "v5",
        "Noutput": 2001,
        "oFileName": "conc.txt",
        "dFileName": "Deff.txt",
        "ParaDictName": "para_dict.csv",
        "SimDictName": "sim_dict.csv",
        "GeoDictName": "geo_dict.csv",
        "SaveFile": True
    }

    ccmParaDict = {
        "Ccyt": 1e-5,
        "kc": 0.0003,
        "Delta": 10,  # in nanometer
        "TubGeometry": "extended",
        "ThyGeometry": "chlamy",
        "RubGeo": "pyrenoid",
        "Rub_s": 0.0,  # default 0.0
        "Rub_e": ccm_constants.R_pyr,  # default ccm_constants.R_pyr
        "LCIB_s": ccm_constants.R_pyr,  # default ccm_constants.R_pyr
        "LCIB_e": ccm_constants.R_chlor,  # default ccm_constants.R_chlor
        "CAH3_s": 0.0,  # default 0.0
        "CAH3_e": ccm_constants.R_pyr,  # default ccm_constants.R_pyr
        "BST_s": 0.0,  # default 0.0
        "BST_e": ccm_constants.R_chlor,  # default ccm_constants.R_chlor
        "LCIB": "on",
        "vLCIB": 1e3,
        "vCAH3": 1e4,
        "TubTrans": "channel",
        "ChlorTrans": "active",
        "StromaDiffB": "off",
        "PyrMemDiffB": "on",  # off is mostly permeable
        "vTub": 1e-3,
        "vChlor": 1e-2,
        "gamma": pow(10, -1.5),
    }

    ccmGeoDict = {
        "Rubisco_Geometry": {
            "r_start": 0.0,
            "r_end": 0.0
        },
        "LCIB_Geometry": {
            "r_start": 0.0,
            "r_end": 0.0
        },
        "CAH3_Geometry": {
            "r_start": 0.0,
            "r_end": 0.0
        },
        "BST_Geometry": {
            "r_start": 0.0,
            "r_end": 0.0
        },
        "Starch_Permeability": {
            "kcPyrMem": 0.0,
            "khPyrMem": 0.0
        }
    }

    namefolder = makePath(ccmSimDict, ccmParaDict)
    print(namefolder + "\n")

    R_tub = ccm_domain.get_R_tub(ccmParaDict["TubGeometry"])

    mesh, domain_all, boundaries = ccm_domain.genMesh(ccm_constants.R_chlor,
                                                      Nmesh, R_tub,
                                                      TubuleMarker, PyrMarker,
                                                      PyrMemMarker,
                                                      ChlorMemMarker)

    # + volume/surface integration element
    dx = Measure("dx", domain=mesh, subdomain_data=domain_all)
    ds = Measure("ds", domain=mesh, subdomain_data=boundaries)
    dS = Measure("dS", domain=mesh, subdomain_data=boundaries)

    # + cytosolic CO2 and bicarb
    C_cyt = ccmParaDict["Ccyt"] / ccm_constants.Conc
    H_cyt = 10 * C_cyt

    # rescale rate constant
    k_c = ccmParaDict["kc"] / (ccm_constants.Lnth / ccm_constants.Time)
    v_tub = ccmParaDict["vTub"] / (ccm_constants.Lnth / ccm_constants.Time)
    v_chlor = ccmParaDict["vChlor"] / (ccm_constants.Lnth / ccm_constants.Time)
    D_c_slow = interpDeff(ccmParaDict["Delta"], k_c, ccm_constants.D_c)
    D_h_slow = interpDeff(ccmParaDict["Delta"], v_tub, ccm_constants.D_h_chlor)
    v_cah3 = ccmParaDict["vCAH3"]
    v_lcib = ccmParaDict["vLCIB"]

    # parse geometry input
    ccmGeoDict["Starch_Permeability"]["kcPyrMem"], ccmGeoDict[
        "Starch_Permeability"]["khPyrMem"] = get_permeability(
            ccmParaDict["PyrMemDiffB"])
    ccmGeoDict["Rubisco_Geometry"]["r_start"] = ccmParaDict["Rub_s"]
    ccmGeoDict["Rubisco_Geometry"]["r_end"] = ccmParaDict["Rub_e"]
    ccmGeoDict["LCIB_Geometry"]["r_start"] = ccmParaDict["LCIB_s"]
    ccmGeoDict["LCIB_Geometry"]["r_end"] = ccmParaDict["LCIB_e"]
    ccmGeoDict["CAH3_Geometry"]["r_start"] = ccmParaDict["CAH3_s"]
    ccmGeoDict["CAH3_Geometry"]["r_end"] = ccmParaDict["CAH3_e"]
    ccmGeoDict["BST_Geometry"]["r_start"] = ccmParaDict["BST_s"]
    ccmGeoDict["BST_Geometry"]["r_end"] = ccmParaDict["BST_e"]

    # if enzyme domain varies from default, rescale rate
    # Rubisco rate is rescaled in ccm_eqn and ccm_erg
    fv_tub, fs_tub = ccm_eqn.get_geo_factors(ccmParaDict["ThyGeometry"], R_tub)
    if (ccmParaDict["LCIB_s"] != ccm_constants.R_pyr) or (
            ccmParaDict["LCIB_e"] != ccm_constants.R_chlor):
        v_lcib = ccmParaDict["vLCIB"] * get_rescale_factor(
            dx, 1.0 - fv_tub, ccm_constants.R_pyr, ccm_constants.R_chlor,
            ccmParaDict["LCIB_s"], ccmParaDict["LCIB_e"])
    if (ccmParaDict["BST_s"] != 0.0) or (ccmParaDict["BST_e"] !=
                                         ccm_constants.R_chlor):
        v_tub = ccmParaDict["vTub"] * get_rescale_factor(
            dx, fv_tub * fs_tub, 0.0, ccm_constants.R_chlor,
            ccmParaDict["BST_s"], ccmParaDict["BST_e"])
    if (ccmParaDict["CAH3_s"] != 0.0) or (ccmParaDict["CAH3_e"] !=
                                          ccm_constants.R_pyr):
        v_cah3 = ccmParaDict["vCAH3"] * get_rescale_factor(
            dx, fv_tub, 0.0, ccm_constants.R_pyr, ccmParaDict["CAH3_s"],
            ccmParaDict["CAH3_e"])

    # + define functionspace and variable
    P1 = FiniteElement('P', mesh.ufl_cell(), 1)
    element = MixedElement(
        [P1, P1, P1, P1, P1,
         P1])  # 1 = tubule CO2, 2 = pyrenoid CO2, 3 = stroma CO2
    # 4 = tubule bicarb, 5 = pyrenoid bicarb, 6 = stroma bicarb
    V = FunctionSpace(mesh, element)
    VDG = FunctionSpace(mesh, "DG", 0)

    u = Function(V)

    # + assemble eqns
    L_tot, Jacobian = ccm_eqn.gen_weak_form_eqn(u, VDG, k_c, D_c_slow, D_h_slow,\
                      R_tub, C_cyt, H_cyt,\
                      ccmParaDict["TubTrans"], ccmParaDict["ChlorTrans"],\
                      v_tub, v_chlor, ccmParaDict["gamma"],\
                      ccmParaDict["LCIB"], ccmParaDict["vLCIB"],ccmParaDict["vCAH3"],\
                      ccmParaDict["StromaDiffB"], ccmGeoDict, ccmParaDict["ThyGeometry"], ccmParaDict["RubGeo"],\
                      dx, ds, dS, pyr_marker = PyrMarker, tubule_marker=TubuleMarker,\
                      pyr_mem_marker=PyrMemMarker, chlor_mem_marker=ChlorMemMarker)

    class ReactionDiffusionProblem(NonlinearProblem):
        def __init__(self, L, a):
            NonlinearProblem.__init__(self)
            self.L = L
            self.a = a

        def F(self, b, x):
            assemble(self.L, tensor=b)

        def J(self, A, x):
            assemble(self.a, tensor=A)
            A.ident_zeros()

    problem = ReactionDiffusionProblem(L_tot, Jacobian)

    # + set up solvers
    solver = NewtonSolver()
    solver.parameters["linear_solver"] = "lu"
    solver.parameters["convergence_criterion"] = "incremental"
    solver.parameters["relative_tolerance"] = 1e-8
    solver.parameters["maximum_iterations"] = 100
    solver.parameters["report"] = True

    # + initialization
    u.vector()[:] = 0.1

    # + solve eqns
    solver.solve(problem, u.vector())
    (c1p, c2p, c2p_, h1p, h2p, h2p_) = u.split(deepcopy=True)

    print("Complete calculating concentration profiles...\n")

    # + save concentration profile
    if ccmSimDict["SaveFile"] == True:

        ofile = open(namefolder + "/" + ccmSimDict["oFileName"], 'w')
        for x_ in np.linspace(0, ccm_constants.R_chlor, ccmSimDict["Noutput"]):
            ofile.write(
                str(x_) + ' ' + str(c1p(x_)) + ' ' + str(c2p(x_)) + ' ' +
                str(c2p_(x_)) + ' ' + str(h1p(x_)) + ' ' + str(h2p(x_)) + ' ' +
                str(h2p_(x_)) + '\n')

        ofile.close()

        ofile2 = open(namefolder + "/" + ccmSimDict["dFileName"], 'w')
        ofile2.write(str(D_c_slow) + ' ' + str(D_h_slow))
        ofile2.close()

        print("Complete saving data...\n")

    # + save simulation settings

    if ccmSimDict["SaveFile"] == True:

        saveTOcsv(ccmSimDict, namefolder + ccmSimDict["SimDictName"])
        saveTOcsv(ccmParaDict, namefolder + ccmSimDict["ParaDictName"])
        saveTOcsv(ccmGeoDict, namefolder + ccmSimDict["GeoDictName"])

        print("Complete saving simulation settings...\n")

    print("Finished!")
