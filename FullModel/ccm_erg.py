# ==================================================
# ============= A. import packages =================
from dolfin import *
import importlib.util
import numpy as np
import ufl
import os
import sys
set_log_active(False)

cdir = os.path.dirname(os.path.abspath(__file__))
spec_const = importlib.util.spec_from_file_location("ccm_constants",
                                                    cdir + "/ccm_constants.py")
ccm_constants = importlib.util.module_from_spec(spec_const)
spec_const.loader.exec_module(ccm_constants)

parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["quadrature_degree"] = 4

# ==================================================
# ==================================================


# + compute Rubisco flux
def get_Rubisco_flux(c2, c2_, r, dx, fv_tub, RubGeoDict):

    fv_rub = conditional(
        le(r, RubGeoDict["r_end"]),
        conditional(ge(r, RubGeoDict["r_start"]), 1. - fv_tub, Constant(0.0)),
        Constant(0.0))
    fv_pyr    = conditional(le(r, ccm_constants.R_pyr), 1 - Expression("x[0] > r_tub ? 0.0 : (x[0] > rin ? fv_in*rin*rin/x[0]/x[0] : fv_in)",\
            r_tub=1.0, fv_in=ccm_constants.fv_in, rin=ccm_constants.rin, degree = 1) ,Constant(0.0))

    vol_rub = assemble(r**2 * fv_rub * dx)
    vol_pyr = assemble(r**2 * fv_pyr * dx)

    Vmax_correct = ccm_constants.Vmax * vol_pyr / vol_rub

    f_correction = conditional(
        ge(r, RubGeoDict["r_start"]),
        conditional(le(r, RubGeoDict["r_end"]), Constant(1.0), Constant(0.0)),
        Constant(0.0))

    R_rub = conditional(le(r, ccm_constants.R_pyr), (Vmax_correct * c2) /
                        (ccm_constants.Km_eff + c2) * f_correction,
                        Constant(0.0))
    R_rub_ = conditional(ge(r, ccm_constants.R_pyr), (Vmax_correct * c2_) /
                         (ccm_constants.Km_eff + c2_) * f_correction,
                         Constant(0.0))

    return R_rub, R_rub_, fv_rub


# + compute geometrical factors
def get_geo_factors(thy_geo, R_tub):
    if thy_geo == "plant":
        fv_tub = Expression("x[0] > r_tub ? 0.0 : fv_in",\
            r_tub=R_tub, fv_in=ccm_constants.fv_plant, degree = 1)
        fs_tub = Constant(2. /
                          ccm_constants.a_plant)  # = ds / (4*pi*r^2*dr*fv_tub)
    else:  # use chlamy geometry
        fv_tub = Expression("x[0] > r_tub ? 0.0 : (x[0] > rin ? fv_in*rin*rin/x[0]/x[0] : fv_in)",\
            r_tub=R_tub, fv_in=ccm_constants.fv_in, rin=ccm_constants.rin, degree = 1)
        fs_tub = Constant(2. / ccm_constants.a_tub)
    return fv_tub, fs_tub


# + compute the transportation flux of transporter/channel
def get_tub_mem_flux(h1, h2, h2_, r, type, rate_tub, gamma, bstGeoDict):
    if type == "none":
        return Constant(0.)
    elif type == "channel":
        jbst_null = conditional(
            le(r, ccm_constants.R_pyr),
            rate_tub *
            (h2 * ccm_constants.f_chlor /
             (1 + ccm_constants.f_chlor) - h1 * ccm_constants.f_tub /
             (1 + ccm_constants.f_tub)),
            rate_tub *
            (h2_ * ccm_constants.f_chlor /
             (1 + ccm_constants.f_chlor) - h1 * ccm_constants.f_tub /
             (1 + ccm_constants.f_tub)))
        jbst_new  = conditional(le(r,bstGeoDict["r_end"]),\
                                conditional(ge(r,bstGeoDict["r_start"]),jbst_null, Constant(0.0)),\
                                Constant(0.0))
        return jbst_new
    elif type == "active":
        jbst_null = conditional(
            le(r, ccm_constants.R_pyr),
            rate_tub *
            (h2 * ccm_constants.f_chlor /
             (1 + ccm_constants.f_chlor) - gamma * h1 * ccm_constants.f_tub /
             (1 + ccm_constants.f_tub)),
            rate_tub *
            (h2_ * ccm_constants.f_chlor /
             (1 + ccm_constants.f_chlor) - gamma * h1 * ccm_constants.f_tub /
             (1 + ccm_constants.f_tub)))
        jbst_new  = conditional(le(r,bstGeoDict["r_end"]),\
                                conditional(ge(r,bstGeoDict["r_start"]),jbst_null, Constant(0.0)),\
                                Constant(0.0))
        return jbst_new
    else:
        return Constant(0.)


def get_chlor_mem_flux(h2_, H_cyt, type, rate_chlor, gamma):
    if type == "none":
        return Constant(0.)
    elif type == "channel":
        return rate_chlor * (
            H_cyt - h2_(ccm_constants.R_chlor) * ccm_constants.f_chlor /
            (1 + ccm_constants.f_chlor))
    elif type == "active":
        return rate_chlor * (H_cyt - gamma * h2_(ccm_constants.R_chlor) *
                             ccm_constants.f_chlor /
                             (1 + ccm_constants.f_chlor))
    else:
        return Constant(0.)


def get_lcib_flux(status, r, c2, h2, c2_, h2_, v_lcib, lcibGeoDict):
    kmax_C_lcib = v_lcib / (1 / ccm_constants.Time)
    kmax_H_chlor = kmax_C_lcib * pow(
        10, ccm_constants.pK_eff - ccm_constants.pH_chlor)
    if status == "on":
        jlcib_null = conditional(le(r, ccm_constants.R_pyr),\
                        (kmax_C_lcib * c2  - kmax_H_chlor * h2  * ccm_constants.f_chlor / (1+ccm_constants.f_chlor))  / (1 + c2 /ccm_constants.Km_CO2 + h2  * (ccm_constants.f_chlor / (1+ccm_constants.f_chlor)) / ccm_constants.Km_HCO3),\
                        (kmax_C_lcib * c2_ - kmax_H_chlor * h2_ * ccm_constants.f_chlor / (1+ccm_constants.f_chlor))  / (1 + c2_/ccm_constants.Km_CO2 + h2_ * (ccm_constants.f_chlor / (1+ccm_constants.f_chlor)) / ccm_constants.Km_HCO3))
        jlcib_new  = conditional(le(r,lcibGeoDict["r_end"]),\
                                conditional(ge(r,lcibGeoDict["r_start"]),jlcib_null, Constant(0.0)),\
                                Constant(0.0))
    else:
        jlcib_new = Constant(0.0)

    jsp = conditional(le(r, ccm_constants.R_pyr),\
                  ccm_constants.v_sp * c2  - ccm_constants.v_sp * pow(10, ccm_constants.pK_eff - ccm_constants.pH_chlor) * h2  * ccm_constants.f_chlor / (1+ccm_constants.f_chlor),\
                  ccm_constants.v_sp * c2_ - ccm_constants.v_sp * pow(10, ccm_constants.pK_eff - ccm_constants.pH_chlor) * h2_ * ccm_constants.f_chlor / (1+ccm_constants.f_chlor))
    return (jlcib_new + jsp)


def get_cah3_flux(r, c1, h1, v_cah3, cah3GeoDict):
    r_start = cah3GeoDict["r_start"]
    r_end = cah3GeoDict["r_end"]
    kmax_C = v_cah3 / (1 / ccm_constants.Time)
    kmax_H_tub = kmax_C * pow(10, ccm_constants.pK_eff - ccm_constants.pH_tub)
    jcah3 = conditional(
        le(r, r_end),
        conditional(ge(r, r_start),
                    (kmax_C * c1 - kmax_H_tub * h1 * ccm_constants.f_tub /
                     (1 + ccm_constants.f_tub)) /
                    (1 + c1 / ccm_constants.Km_CO2 + h1 * ccm_constants.f_tub /
                     (1 + ccm_constants.f_tub) / ccm_constants.Km_HCO3),
                    Constant(0.0)), Constant(0.0))
    jsp = (ccm_constants.v_sp * c1 - ccm_constants.v_sp *
           pow(10, ccm_constants.pK_eff - ccm_constants.pH_tub) * h1 *
           ccm_constants.f_tub / (1 + ccm_constants.f_tub))
    return (jcah3 + jsp)




def get_carb_flux(u, VDG, k_c, D_c_slow, D_h_slow, R_tub, C_cyt, H_cyt,\
                  tub_trans_type, chlor_trans_type, v_tub, v_chlor, gamma,\
                  stroma_ca_stat, v_lcib, v_cah3,\
                  stroma_barrier, ccmGeoDict, thy_geo, rub_geo, output_type,\
                  dx, ds, dS, pyr_marker, tubule_marker, pyr_mem_marker, chlor_mem_marker ):

    V = u.function_space()

    # Define trialfunction and testfunction
    du = TrialFunction(V)
    v = TestFunction(V)
    (v0, v1, v1_, v2, v3, v3_) = split(v)

    # Split system functions to access components
    (c1, c2, c2_, h1, h2,
     h2_) = u.split(deepcopy=True)  # c=[CO2], h=[HCO3-], 1=tubule, 2=stroma

    r = Expression("x[0]", degree=1)
    k_c_pyrmem = ccmGeoDict["Starch_Permeability"]["kcPyrMem"] / (
        ccm_constants.Lnth / ccm_constants.Time)
    k_h_pyrmem = ccmGeoDict["Starch_Permeability"]["khPyrMem"] / (
        ccm_constants.Lnth / ccm_constants.Time)

    # + Geometric correction
    fv_tub, fs_tub = get_geo_factors(thy_geo, R_tub)
    fs_stroma = (fv_tub / (1 - fv_tub)) * fs_tub  # for the stroma

    # ===================================================================== #
    # rate of consumption of CO2 by Rubisco
    R_rub, R_rub_, fv_rub = get_Rubisco_flux(c2, c2_, r, dx, fv_tub,
                                             ccmGeoDict["Rubisco_Geometry"])
    # ===================================================================== #

    # ===================================================================== #
    # flux of CO2 from tubules to pyrenoid
    jc_12 = conditional(
        le(r, R_tub),
        conditional(le(r, ccm_constants.R_pyr), k_c * (c1 - c2),
                    k_c * (c1 - c2_)), Constant(0.0))
    # ===================================================================== #

    # ===================================================================== #
    # + flux of bicarb channels/transporters
    j_active_tub = get_tub_mem_flux(h1, h2, h2_, r, tub_trans_type, v_tub,
                                    gamma, ccmGeoDict["BST_Geometry"])
    j_active_chlor = get_chlor_mem_flux(h2_, H_cyt, chlor_trans_type, v_chlor,
                                        gamma)


    jh_passive = conditional(le(r,ccm_constants.R_pyr),\
                            ccm_constants.k_h2co3 * ((h1 / (1 + ccm_constants.f_tub)) - (h2 /(1 + ccm_constants.f_chlor))) + ccm_constants.k_hco3 * ((h1 * ccm_constants.f_tub / (1 + ccm_constants.f_tub)) - (h2 * ccm_constants.f_chlor / (1 + ccm_constants.f_chlor))),\
                            ccm_constants.k_h2co3 * ((h1 / (1 + ccm_constants.f_tub)) - (h2_ /(1 + ccm_constants.f_chlor))) + ccm_constants.k_hco3 * ((h1 * ccm_constants.f_tub / (1 + ccm_constants.f_tub)) - (h2_ * ccm_constants.f_chlor / (1 + ccm_constants.f_chlor))))

    # flux of H2CO3 from tubules to pyrenoid
    jh_12 = conditional(le(r, R_tub), jh_passive - j_active_tub, Constant(0.0))
    # ===================================================================== #

    # ===================================================================== #
    # rate of tubule CO2 conversion catalyzed by CA (full MM kinetics)
    j_ca = get_cah3_flux(r, c1, h1, v_cah3, ccmGeoDict["CAH3_Geometry"])
    # ===================================================================== #

    # ===================================================================== #
    # rate of stromal CO2 conversion catalyzed by CA (assume same kinetics)
    j_lcib = get_lcib_flux(stroma_ca_stat, r, c2, h2, c2_, h2_, v_lcib,
                           ccmGeoDict["LCIB_Geometry"])
    # ===================================================================== #

    # differential diffusion:
    if stroma_barrier == "on" and thy_geo == "chlamy" and rub_geo == "pyrenoid":
        D_c_stroma = Expression("x[0] > R_pyr ? D_c_out : D_c_pyr",
                                R_pyr=ccm_constants.R_pyr,
                                D_c_out=D_c_slow,
                                D_c_pyr=ccm_constants.D_c,
                                degree=0)
        D_h_stroma = Expression("x[0] > R_pyr ? D_h_out : D_h_pyr",
                                R_pyr=ccm_constants.R_pyr,
                                D_h_out=D_h_slow,
                                D_h_pyr=ccm_constants.D_h_chlor,
                                degree=0)
    elif stroma_barrier == "on" and thy_geo == "chlamy" and rub_geo == "diffuse":
        D_c_stroma = D_c_slow
        D_h_stroma = D_h_slow
    else:
        D_c_stroma = ccm_constants.D_c
        D_h_stroma = ccm_constants.D_h_chlor

    # compute CO2 fixation index
    c_pyr_avg = (assemble(c2 * r**2 * fv_rub * dx(pyr_marker)) + assemble(
        c2_ * r**2 * fv_rub *
        (dx(tubule_marker) + dx(0)))) / assemble(r**2 * fv_rub * dx)
    perc_max_flux = (
        assemble(c2 /
                 (c2 + ccm_constants.Km_eff) * r**2 * fv_rub * dx(pyr_marker))
        + assemble(c2_ / (c2_ + ccm_constants.Km_eff) * r**2 * fv_rub *
                   (dx(tubule_marker) + dx(0)))) / assemble(r**2 * fv_rub * dx)

    # compute flux terms
    """
    [Note flux 14 and 15 are non-zero only in the modified geometry]
    
    Flux diagram is given  by
    
    A --------------------- B
    |\                     /|
    | \                   / |
    |  C --------------- D  |
    |  |                 |  |
    |  |                 |  |
    |  |                 |  |
    |  E --------------- F--------- I
    | /                   \ |
    |/                     \|
    G --------------------- H ------ J
    |
    |
    X (CO2 FIXATION)
    
    where the nodes are defined as follows:
    
    A = Inner tubule CO2
    B = Outer tubule CO2
    C = Inner tubule HCO3-
    D = Outer tubule HCO3-
    E = Pyrenoid HCO3-
    F = Stroma   HCO3-
    G = Pyrenoid CO2
    H = Stroma   CO2
    I = Cytosol  HCO3-
    J = Cytosol  CO2
    
    and the edges are defined as follows:
    
    1. A-B = CO2 diffusion in tubule
    2. C-A = Cah3 conversion
    3. A-G = CO2 diffusion to pyrenoid
    4. G-X = CO2 fixation
    5. G-H = CO2 diffusion from pyrenoid to stroma
    6. D-C = HCO3- diffusion in the tubule
    7. E-C = HCO3- diffusion from pyrenoid to tubule
    8. F-E = HCO3- diffusion to pyrenoid from stroma
    9. F-D = HCO3- transport to tubule
    10.B-H = CO2 diffusion to stroma
    11.H-F = LCIB conversion
    12.I-F = HCO3- import
    13.H-J = CO2   leakage
    14.D-B = Cah3 conversion (new)
    15.G-E = LCIB conversion (new)
    
    Thus, flux balance of the 8 nodes inside the chloroplast is given by:
    
    A: -(flux1) - (flux3) + (flux2) = 0
    B: +(flux1) - (flux10)+ [flux14]= 0
    C: -(flux2) + (flux6) + (flux7) = 0
    D: -(flux6) + (flux9) - [flux14]= 0
    E: +(flux8) - (flux7) + [flux15]= 0
    F: +(flux11)+ (flux12)- (flux8) - (flux9) = 0
    G: +(flux3) - (flux4) - (flux5) - [flux15]= 0
    H: +(flux5) + (flux10)- (flux11)- (flux13)= 0
    
    """

    # flux1
    #    dc1dx = project(grad(c1)[0],VDG)
    dc1dx = project(grad(c1)[0], c1.function_space())
    flux1 = -4 * pi * ccm_constants.R_pyr**2 * fv_tub(
        ccm_constants.R_pyr) * ccm_constants.D_c * dc1dx(ccm_constants.R_pyr)

    # * Cah3 conversion flux
    # flux2
    flux2 = -assemble(4 * pi * r**2 * fv_tub * j_ca * dx(pyr_marker))
    # flux14
    flux14 = -assemble(4 * pi * r**2 * fv_tub * j_ca *
                       (dx(tubule_marker) + dx(0)))

    # flux3
    flux3 = assemble(4 * pi * r**2 * fv_tub * fs_tub * jc_12 * dx(pyr_marker))

    # flux4
    flux4 = assemble(R_rub * 4 * pi * r**2 * (1.0 - fv_tub) * dx(pyr_marker))

    # flux5
    flux5 = 4 * pi * ccm_constants.R_pyr**2 * (1 - fv_tub(
        ccm_constants.R_pyr)) * k_c_pyrmem * (c2(ccm_constants.R_pyr) -
                                              c2_(ccm_constants.R_pyr))

    # flux6
    dh1dx = project(grad(h1)[0], c1.function_space())
    flux6 = 4 * pi * ccm_constants.R_pyr**2 * fv_tub(
        ccm_constants.R_pyr) * ccm_constants.D_h_tub * dh1dx(
            ccm_constants.R_pyr)

    # flux7
    flux7 = -assemble(4 * pi * r**2 * fv_tub * fs_tub * jh_12 * dx(pyr_marker))

    # flux 8
    flux8 = 4 * pi * ccm_constants.R_pyr**2 * (1 - fv_tub(
        ccm_constants.R_pyr)) * k_h_pyrmem * (h2_(ccm_constants.R_pyr) -
                                              h2(ccm_constants.R_pyr))

    # flux9
    flux9 = -assemble(4 * pi * r**2 * fv_tub * fs_tub * jh_12 *
                      (dx(tubule_marker) + dx(0)))

    # flux10
    flux10 = assemble(4 * pi * r**2 * fv_tub * fs_tub * jc_12 *
                      (dx(tubule_marker) + dx(0)))

    # * LCIB conversion flux
    # flux11
    flux11 = assemble(4 * pi * r**2 * (1.0 - fv_tub) * j_lcib *
                      (dx(tubule_marker) + dx(0)))
    # flux15
    flux15 = assemble(4 * pi * r**2 * (1.0 - fv_tub) * j_lcib * dx(pyr_marker))

    # flux12
    flux12_a = get_chlor_mem_flux(h2_, H_cyt, chlor_trans_type, v_chlor, gamma)
    flux12_p = (ccm_constants.k_h2co3 * ((H_cyt / ccm_constants.f_cyt) - (h2_(ccm_constants.R_chlor) /(1 + ccm_constants.f_chlor) )) ) \
        + (ccm_constants.k_hco3  * ((H_cyt) - (h2_(ccm_constants.R_chlor) * ccm_constants.f_chlor / (1 + ccm_constants.f_chlor))))

    flux12 = 4 * pi * ccm_constants.R_chlor**2 * (flux12_a + flux12_p)

    # flux13
    flux13 = 4 * pi * ccm_constants.R_chlor**2 * k_c * (
        c2_(ccm_constants.R_chlor) - C_cyt)

    # summary flux
    # (i)   Total bicarb uptake
    J_H_0 = 4 * pi * ccm_constants.R_chlor**2 * (ccm_constants.k_h2co3 * (
        (H_cyt / ccm_constants.f_cyt) - (h2_(ccm_constants.R_chlor) /
                                         (1 + ccm_constants.f_chlor))))
    J_H_diff = 4 * pi * ccm_constants.R_chlor**2 * (ccm_constants.k_hco3 * (
        (H_cyt) - (h2_(ccm_constants.R_chlor) * ccm_constants.f_chlor /
                   (1 + ccm_constants.f_chlor))))
    J_H_tot = flux12
    # (ii)  Total CO2 leakage
    J_C = flux13
    # (iii) Total Rubisco fixation flux
    J_rub = flux4
    # (iv)  Total LCIB flux
    J_LCIB = flux11 + flux15
    # (v)   Total Cah3 flux
    J_CAH3 = flux2 + flux14

    # energy calculation
    Gatp = 51.5 / 2.48
    gamma_chlor = 1
    if chlor_trans_type == "active":
        gamma_chlor = gamma
    Erg_exact = (J_LCIB / J_rub) * np.log(100) / Gatp + (J_C / J_rub) * np.log(
        100 / gamma_chlor) / Gatp + 2.5 + np.log(100 / gamma_chlor) / Gatp + (
            J_H_diff / J_rub) * np.log(gamma_chlor) / Gatp - (
                J_H_0 / J_rub) * np.log(10 / gamma_chlor * pow(
                    10, ccm_constants.pKa1 - ccm_constants.pH_cyt)) / Gatp
    Erg_estimate = (J_LCIB / J_rub) * np.log(100) / Gatp + (
        J_C / J_rub) * np.log(100 / gamma_chlor) / Gatp + 2.5 + np.log(
            100 / gamma_chlor) / Gatp

    # average concentration
    # 1) inner tubule; 2) outer tubule; 3) stroma; 4) pyrenoid
    c_vec_1 = assemble(
        4 * pi * r**2 * fv_tub * c1 * dx(pyr_marker)) / assemble(
            4 * pi * r**2 * fv_tub * dx(pyr_marker))
    c_vec_2 = assemble(4 * pi * r**2 * fv_tub * c1 *
                       (dx(tubule_marker) + dx(0))) / assemble(
                           4 * pi * r**2 * fv_tub *
                           (dx(tubule_marker) + dx(0)))
    c_vec_3 = assemble(4 * pi * r**2 * (1.0 - fv_tub) * c2_ *
                       (dx(tubule_marker) + dx(0))) / assemble(
                           4 * pi * r**2 * (1.0 - fv_tub) *
                           (dx(tubule_marker) + dx(0)))
    c_vec_4 = assemble(4 * pi * r**2 *
                       (1.0 - fv_tub) * c2 * dx(pyr_marker)) / assemble(
                           4 * pi * r**2 * (1.0 - fv_tub) * dx(pyr_marker))

    h_vec_1 = assemble(4 * pi * r**2 * fv_tub * h1 * ccm_constants.f_tub /
                       (1 + ccm_constants.f_tub) * dx(pyr_marker)) / assemble(
                           4 * pi * r**2 * fv_tub * dx(pyr_marker))
    h_vec_2 = assemble(4 * pi * r**2 * fv_tub * h1 * ccm_constants.f_tub /
                       (1 + ccm_constants.f_tub) *
                       (dx(tubule_marker) + dx(0))) / assemble(
                           4 * pi * r**2 * fv_tub *
                           (dx(tubule_marker) + dx(0)))
    h_vec_3 = assemble(
        4 * pi * r**2 * (1.0 - fv_tub) * h2_ * ccm_constants.f_chlor /
        (1 + ccm_constants.f_chlor) *
        (dx(tubule_marker) + dx(0))) / assemble(4 * pi * r**2 *
                                                (1.0 - fv_tub) *
                                                (dx(tubule_marker) + dx(0)))
    h_vec_4 = assemble(
        4 * pi * r**2 * (1.0 - fv_tub) * h2 * ccm_constants.f_chlor /
        (1 + ccm_constants.f_chlor) * dx(pyr_marker)) / assemble(
            4 * pi * r**2 * (1.0 - fv_tub) * dx(pyr_marker))

    if output_type == "screen":

        print("\n\n")
        print("CO2 fixation and its energetic cost")
        print("-----------------------------------")
        print("Pyrenoid CO2 conc (mM): %f" % c_pyr_avg)
        print("Norm. CO2 fixation flux: %f" % perc_max_flux)
        # Chenyi: add net flux calculation?
        #print("Net. CO2 fixation flux: %f" % net_flux)
        print("ATP per CO2 fixed: %.2f" % (Erg_estimate - 2.5))

        print("\n\n")
        print("Flux balance in different compartments")
        print("--------------------------------------")

        print("C1: %f" % (-flux1 - flux3 + flux2))
        print("C2: %f" % (+flux1 - flux10 + flux14))
        print("C3: %f" % (+flux5 + flux10 - flux11 - flux13))
        print("C4: %f" % (+flux3 - flux4 - flux5 - flux15))

        print("H1: %f" % (-flux2 + flux6 + flux7))
        print("H2: %f" % (-flux6 + flux9 - flux14))
        print("H3: %f" % (+flux11 + flux12 - flux8 - flux9))
        print("H4: %f" % (+flux8 - flux7 + flux15))
        print("\n\n")

        return [1.0]

    elif output_type == "flux":
        return [
            flux1, flux2, flux3, flux4, flux5, flux6, flux7, flux8, flux9,
            flux10, flux11, flux12, flux13, flux14, flux15
        ]

    elif output_type == "erg":
        return [
            c_pyr_avg, perc_max_flux, J_H_0, J_H_diff, J_H_tot, J_C, J_rub,
            J_LCIB, J_CAH3, Erg_exact, Erg_estimate
        ]

    elif output_type == "conc":
        return [
            c_vec_1, c_vec_2, c_vec_3, c_vec_4, h_vec_1, h_vec_2, h_vec_3,
            h_vec_4
        ]

    else:
        return [0.0]
