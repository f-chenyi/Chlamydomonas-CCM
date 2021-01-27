'''
This script allows for a 3d parameter scan of the model.

(when I'm done with this I can make the same edits to Main.py,
which will allow for single quick runs of one set of params.)

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
# import ccm_constants
spec_const = importlib.util.spec_from_file_location("ccm_constants",
                                                    cdir + "/ccm_constants.py")
ccm_constants = importlib.util.module_from_spec(spec_const)
spec_const.loader.exec_module(ccm_constants)
# import ccm_domain
spec_domain = importlib.util.spec_from_file_location("ccm_domain",
                                                     cdir + "/ccm_domain.py")
ccm_domain = importlib.util.module_from_spec(spec_domain)
spec_domain.loader.exec_module(ccm_domain)
# import ccm_eqn
spec_eqn = importlib.util.spec_from_file_location("ccm_eqn",
                                                  cdir + "/ccm_eqn.py")
ccm_eqn = importlib.util.module_from_spec(spec_eqn)
spec_eqn.loader.exec_module(ccm_eqn)
# import ccm_main
spec_main = importlib.util.spec_from_file_location("ccm_main",
                                                   cdir + "/ccm_main.py")
ccm_main = importlib.util.module_from_spec(spec_main)
spec_main.loader.exec_module(ccm_main)
# import ccm_erg
spec_erg = importlib.util.spec_from_file_location("ccm_erg",
                                                  cdir + "/ccm_erg.py")
ccm_erg = importlib.util.module_from_spec(spec_erg)
spec_erg.loader.exec_module(ccm_erg)

# # + parameter setting
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["quadrature_degree"] = 4
parameters["form_compiler"]["representation"] = "uflacs"
# ==================================================
# ==================================================


def main():
    return


def makeDir(namefolder):

    if os.path.exists(namefolder):
        os.stat(namefolder)
    else:
        os.makedirs(namefolder)


def saveTOcsv(dict_name, file_name):
    fields = []
    vals = []
    wfile = csv.writer(open(file_name, "w"))
    for key in dict_name:
        fields.append(key)
        vals.append(dict_name[key])
    wfile.writerow(fields)
    wfile.writerow(vals)


def makeFilePath(paraScan, problemDict):
    """
    It creates the name of the output folder and the path of the
    output file.
    """
    FigIndex = 'null'
    if problemDict["ChlorTrans"] == "channel":
        FigIndex = 'passive_'
    if problemDict["ChlorTrans"] == "active":
        FigIndex = 'active_'
    for ScanName in paraScan:
        FigIndex += ScanName
        FigIndex += "_"

    pwd = os.getcwd()
    geofolder = ""
    DiffBfolder = ""

    for key in problemDict.keys():
        if "Geo" in key:
            geofolder += key
            if type(problemDict[key]) == str:
                geofolder += "_" + problemDict[key] + "_"
            else:
                geofolder += ("_%.2E_" % problemDict[key])
        if "DiffB" in key:
            DiffBfolder += key
            if type(problemDict[key]) == str:
                DiffBfolder += "_" + problemDict[key] + "_"
            else:
                DiffBfolder += ("_%.2E_" % problemDict[key])

    # make different folders for different output types
    if problemDict["output_type"] == "screen":
        namefolder = pwd + "/screen_data/" + FigIndex + "/" + geofolder + "/" + DiffBfolder + "/"
    elif problemDict["output_type"] == "flux":
        namefolder = pwd + "/flux_data/" + FigIndex + "/" + geofolder + "/" + DiffBfolder + "/"
    elif problemDict["output_type"] == "erg":
        namefolder = pwd + "/erg_data/" + FigIndex + "/" + geofolder + "/" + DiffBfolder + "/"
    elif problemDict["output_type"] == "conc":
        namefolder = pwd + "/conc_data/" + FigIndex + "/" + geofolder + "/" + DiffBfolder + "/"
    else:
        namefolder = pwd + "/para_scan_data/" + FigIndex + "/" + geofolder + "/" + DiffBfolder + "/"

    makeDir(namefolder)
    filename = namefolder

    for key in problemDict.keys():
        if key in paraScan:
            continue
        if "Geo" in key:
            continue
        if "DiffB" in key:
            continue
        if "_s" in key:
            continue
        if "_e" in key:
            continue
        if "output" in key:
            continue
        filename += key
        if type(problemDict[key]) == str:
            filename += "_" + problemDict[key] + "_"
        else:
            filename += ("_%.2E_" % problemDict[key])

    return namefolder, filename + '.txt'


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
        "oFileName": "para.txt",
        "ParaDictName": "para_dict.csv",
        "SimDictName": "sim_dict.csv",
        "GeoDictName": "geo_dict.csv",
        "SaveFile": True,
        "N_para_A": 2,
        "N_para_B": 2,
        "loop_para": "Ccyt",
        "scan_para": ['vLCIB', 'vChlor']
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
        # "screen", "flux", "erg", or "conc"
        "output_type": "erg"
    }

    # store geometric variables for passing to other functions
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

    para_A_list = pow(10, np.linspace(0, 6, ccmSimDict["N_para_A"]))
    para_B_list = pow(10, np.linspace(-6, -2, ccmSimDict["N_para_B"]))

    loop_para_val = pow(10, np.linspace(-5, -5, 1))

    N_test = len(para_A_list) * len(para_B_list)

    for para_val in loop_para_val:

        column = 3
        if ccmParaDict["output_type"] == "screen":
            column = 3
        elif ccmParaDict["output_type"] == "flux":
            column = 17
        elif ccmParaDict["output_type"] == "erg":
            column = 13
        elif ccmParaDict["output_type"] == "conc":
            column = 10

        print("Number of columns:" + str(column))
        data_test = np.zeros([N_test, column])
        count = 0
        ccmParaDict[ccmSimDict["loop_para"]] = para_val
        namefolder, txtofilename = makeFilePath(ccmSimDict["scan_para"],
                                                ccmParaDict)

        for para_A in para_A_list:
            for para_B in para_B_list:

                R_tub = ccm_domain.get_R_tub(ccmParaDict["TubGeometry"])
                #is this allowed
                ccmParaDict[ccmSimDict["scan_para"][0]] = para_A
                ccmParaDict[ccmSimDict["scan_para"][1]] = para_B
                params = [para_A, para_B]

                mesh, domain_all, boundaries = ccm_domain.genMesh(
                    ccm_constants.R_chlor, Nmesh, R_tub, TubuleMarker,
                    PyrMarker, PyrMemMarker, ChlorMemMarker)

                # + volume/surface integration element
                dx = Measure("dx", domain=mesh, subdomain_data=domain_all)
                ds = Measure("ds", domain=mesh, subdomain_data=boundaries)
                dS = Measure("dS", domain=mesh, subdomain_data=boundaries)

                # + cytosolic CO2 and bicarb, change for different cytosolic pH
                C_cyt = ccmParaDict["Ccyt"] / ccm_constants.Conc
                H_cyt = 10 * C_cyt

                k_c = ccmParaDict["kc"] / (ccm_constants.Lnth /
                                           ccm_constants.Time)
                v_tub = ccmParaDict["vTub"] / (ccm_constants.Lnth /
                                               ccm_constants.Time)
                v_chlor = ccmParaDict["vChlor"] / (ccm_constants.Lnth /
                                                   ccm_constants.Time)
                D_c_slow = ccm_main.interpDeff(ccmParaDict["Delta"], k_c,
                                               ccm_constants.D_c)
                D_h_slow = ccm_main.interpDeff(ccmParaDict["Delta"], v_tub,
                                               ccm_constants.D_h_chlor)
                v_cah3 = ccmParaDict["vCAH3"]
                v_lcib = ccmParaDict["vLCIB"]

                # parse geometry input into geo storage
                ccmGeoDict["Starch_Permeability"]["kcPyrMem"], ccmGeoDict[
                    "Starch_Permeability"][
                        "khPyrMem"] = ccm_main.get_permeability(
                            ccmParaDict["PyrMemDiffB"])
                ccmGeoDict["Rubisco_Geometry"]["r_start"] = ccmParaDict[
                    "Rub_s"]
                ccmGeoDict["Rubisco_Geometry"]["r_end"] = ccmParaDict["Rub_e"]
                ccmGeoDict["LCIB_Geometry"]["r_start"] = ccmParaDict["LCIB_s"]
                ccmGeoDict["LCIB_Geometry"]["r_end"] = ccmParaDict["LCIB_e"]
                ccmGeoDict["CAH3_Geometry"]["r_start"] = ccmParaDict["CAH3_s"]
                ccmGeoDict["CAH3_Geometry"]["r_end"] = ccmParaDict["CAH3_e"]
                ccmGeoDict["BST_Geometry"]["r_start"] = ccmParaDict["BST_s"]
                ccmGeoDict["BST_Geometry"]["r_end"] = ccmParaDict["BST_e"]

                # if enzyme domain varies from default, rescale rate
                # Rubisco rate is rescaled in ccm_eqn and ccm_erg
                fv_tub, fs_tub = ccm_eqn.get_geo_factors(
                    ccmParaDict["ThyGeometry"], R_tub)
                if (ccmParaDict["LCIB_s"] != ccm_constants.R_pyr) or (
                        ccmParaDict["LCIB_e"] != ccm_constants.R_chlor):
                    v_lcib = ccmParaDict["vLCIB"] * get_rescale_factor(
                        dx, 1.0 - fv_tub, ccm_constants.R_pyr,
                        ccm_constants.R_chlor, ccmParaDict["LCIB_s"],
                        ccmParaDict["LCIB_e"])
                if (ccmParaDict["BST_s"] != 0.0) or (ccmParaDict["BST_e"] !=
                                                     ccm_constants.R_chlor):
                    v_tub = ccmParaDict["vTub"] * get_rescale_factor(
                        dx, fv_tub * fs_tub, 0.0, ccm_constants.R_chlor,
                        ccmParaDict["BST_s"], ccmParaDict["BST_e"])
                if (ccmParaDict["CAH3_s"] != 0.0) or (ccmParaDict["CAH3_e"] !=
                                                      ccm_constants.R_pyr):
                    v_cah3 = ccmParaDict["vCAH3"] * get_rescale_factor(
                        dx, fv_tub, 0.0, ccm_constants.R_pyr,
                        ccmParaDict["CAH3_s"], ccmParaDict["CAH3_e"])

                # + define functionspace and variable
                P1 = FiniteElement('P', mesh.ufl_cell(), 1)
                element = MixedElement(
                    [P1, P1, P1, P1, P1,
                     P1])  # 1 = tubule CO2, 2 = pyrenoid CO2, 3 = stroma CO2
                # 4 = tubule bicarb, 5 = pyrenoid bicarb, 6 = stroma bicarb
                V = FunctionSpace(mesh, element)
                VDG = FunctionSpace(mesh, "DG", 0)
                u = Function(V)

                # + assemble equations
                L_tot, Jacobian = ccm_eqn.gen_weak_form_eqn(u, VDG, k_c, D_c_slow, D_h_slow,\
                R_tub, C_cyt, H_cyt,\
                ccmParaDict["TubTrans"], ccmParaDict["ChlorTrans"],\
                v_tub, v_chlor, ccmParaDict["gamma"],\
                ccmParaDict["LCIB"], v_lcib, v_cah3, \
                ccmParaDict["StromaDiffB"], ccmGeoDict, ccmParaDict["ThyGeometry"], ccmParaDict["RubGeo"], \
                dx, ds, dS, pyr_marker=PyrMarker, tubule_marker=TubuleMarker, \
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

                # + solve equations
                solver.solve(problem, u.vector())
                (c1p, c2p, c2p_, h1p, h2p, h2p_) = u.split(deepcopy=True)

                # this function can return multiple kinds of output data
                output = ccm_erg.get_carb_flux(u, VDG, k_c, D_c_slow, D_h_slow,\
                            R_tub, C_cyt, H_cyt,\
                            ccmParaDict["TubTrans"], ccmParaDict["ChlorTrans"],\
                            v_tub, v_chlor, ccmParaDict["gamma"],\
                            ccmParaDict["LCIB"], ccmParaDict["vLCIB"],ccmParaDict["vCAH3"],\
                            ccmParaDict["StromaDiffB"], ccmGeoDict, ccmParaDict["ThyGeometry"], ccmParaDict["RubGeo"], \
                            ccmParaDict["output_type"],\
                            dx, ds, dS, pyr_marker = PyrMarker, tubule_marker=TubuleMarker,\
                            pyr_mem_marker=PyrMemMarker, chlor_mem_marker=ChlorMemMarker)

                # + save data
                row = params + output
                for k in np.arange(0, column):
                    data_test[count][k] = row[k]

                count = count + 1
                print("count = %d\n" % count)

        # Write desired output to a file with the variable names in the filename
        ofile = open(txtofilename, 'w')
        for k in np.arange(N_test):
            for kk in np.arange(column):
                ofile.write(str(data_test[k][kk]) + ' ')
                ofile.write('\n')
        ofile.close()

if ccmSimDict["SaveFile"] == True:

    saveTOcsv(ccmSimDict, namefolder + ccmSimDict["SimDictName"])
    saveTOcsv(ccmParaDict, namefolder + ccmSimDict["ParaDictName"])
    saveTOcsv(ccmGeoDict, namefolder + ccmSimDict["GeoDictName"])

    para = []
    if os.path.exists(namefolder + ccmSimDict["oFileName"]):
        fr = open(namefolder + ccmSimDict["oFileName"], 'r')
        for line in fr:
            para_this = float(line)
            para.append(float(para_this))
        fr.close()

    para = para + list(loop_para_val)
    paraset = set(para)
    para = list(paraset)

    parafile = open(namefolder + ccmSimDict["oFileName"], 'w')
    for para_val in para:
        parafile.write(str(para_val) + ' ' + '\n')
    parafile.close()
