"""
This file contains the key equations used for calculating fluxes, 
generating the weak form of the model equation, and computing 
the geometric factors governing the chloroplast. Functions are

get_Rubisco_flux()
get_geo_factors()
get_tub_mem_flux()
get_chlor_mem_flux()
get_lcib_flux()
get_cah3_flux() and
gen_weak_form_eqn() 

all with appropriate arguments. 
"""
# ==================================================
# ============= A. import packages =================
from dolfin import *
import numpy as np
import importlib.util
import ufl
import os
import sys
set_log_active(False)

cdir          = os.path.dirname(os.path.abspath(__file__))
spec_const    = importlib.util.spec_from_file_location("ccm_constants",\
                                                    cdir + "/ccm_constants.py")
ccm_constants = importlib.util.module_from_spec(spec_const)
spec_const.loader.exec_module(ccm_constants)

# ==================================================
# ==================================================


# compute Rubisco flux
def get_Rubisco_flux(c2, c2_, r, dx, fv_tub, RubGeoDict):
    # c2 : co2 conc. inside  the pyrenoid radius
    # c2_: co2 conc. outside the pyenoid radius
    # r  : radial coordinate
    # fv_tub: tubule volume fraction
    
    # Rubisco volume fraction specified by the input dictionary RubGeoDict
    fv_rub = conditional(
        le(r, RubGeoDict["r_end"]),
        conditional(ge(r, RubGeoDict["r_start"]), 1. - fv_tub, Constant(0.0)),
        Constant(0.0)
        )
    # Rubisco volume fraction of the default geometry
    fv_pyr = conditional(
        le(r, ccm_constants.R_pyr),
        1 - Expression("x[0] > r_tub ? 0.0 : (x[0] > rin ? fv_in*rin*rin/x[0]/x[0] : fv_in)",\
        r_tub=1.0, fv_in=ccm_constants.fv_in, rin=ccm_constants.rin, degree = 1),
        Constant(0.0)
        )

    vol_rub = assemble(r**2 * fv_rub * dx)
    vol_pyr = assemble(r**2 * fv_pyr * dx)

    # correction assuming a fixed total number of Rubisco
    Vmax_correct = ccm_constants.Vmax * vol_pyr / vol_rub
    f_correction = conditional(
        ge(r, RubGeoDict["r_start"]),
        conditional(le(r, RubGeoDict["r_end"]), Constant(1.0), Constant(0.0)),
        Constant(0.0)
        )
    
    # flux of CO2 fixation by Rubisco inside the pyrenoid
    R_rub = conditional(le(r, ccm_constants.R_pyr), (Vmax_correct * c2) /
                        (ccm_constants.Km_eff + c2) * f_correction,
                        Constant(0.0))
    # flux of CO2 fixation by Rubisco outside the pyrenoid
    R_rub_ = conditional(ge(r, ccm_constants.R_pyr), (Vmax_correct * c2_) /
                         (ccm_constants.Km_eff + c2_) * f_correction,
                         Constant(0.0))

    return R_rub, R_rub_




# compute thylakoid geometrical factors, i.e.,
# the volume fraction fv_tub, and
# the surface-to-volume ratio fs_tub
def get_geo_factors(thy_geo, R_tub):
    # (optional) plant geometry assuming a meshwork structure
    # of thylakoids throughout the chloroplast
    if thy_geo == "plant":
        fv_tub = Expression("x[0] > r_tub ? 0.0 : fv_in",\
            r_tub=R_tub, fv_in=ccm_constants.fv_plant, degree = 1)
        fs_tub = Constant(2. /ccm_constants.a_plant)
    else:# use chlamy geometry
        fv_tub = Expression("x[0] > r_tub ? 0.0 : (x[0] > rin ? fv_in*rin*rin/x[0]/x[0] : fv_in)",\
            r_tub=R_tub, fv_in=ccm_constants.fv_in, rin=ccm_constants.rin, degree = 1)
        fs_tub = Constant(2. / ccm_constants.a_tub)
    return fv_tub, fs_tub



# compute the bst mediated hco3- flux across the thylakoid membranes
def get_tub_mem_flux(h1, h2, h2_, r, type, rate_tub, gamma, bstGeoDict):
    # h1 : concentration of hco3- + h2co3 in the thylakoid lumen
    # h2 : concentration of hco3- + h2co3 in the pyrenoid matrix
    # h2_: concentration of hco3- + h2co3 in the chloroplast stroma
    # rate_tub: velocity of the bst channel/transporter
    # gamma:    irreversibility of bst transporters
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



# compute LCIA mediated hco3- flux across the chloroplast envelope
def get_chlor_mem_flux(h2_, H_cyt, type, rate_chlor, gamma):
    # H_cyt   : cytosolic hco3- concentration
    # rate_tub: velocity of the LCIA channel/transporter
    # gamma   : irreversibility of LCIA transporters
    if type == "none":
        return Constant(0.)
    elif type == "channel":
        return rate_chlor * (H_cyt - h2_ * ccm_constants.f_chlor /
                             (1 + ccm_constants.f_chlor))
    elif type == "active":
        return rate_chlor * (H_cyt - gamma * h2_ * ccm_constants.f_chlor /
                             (1 + ccm_constants.f_chlor))
    else:
        return Constant(0.)



# compute the co2 to hco3- conversion flux mediated by LCIB
def get_lcib_flux(status, r, c2, h2, c2_, h2_, v_lcib, lcibGeoDict):
    # first-order rate constant of co2 -> hco3- conversion
    kmax_C_lcib = v_lcib / (1 / ccm_constants.Time)
    # first-order rate constant of hco3- -> co2 conversion
    kmax_H_chlor = kmax_C_lcib * pow(10, ccm_constants.pK_eff - ccm_constants.pH_chlor)
    if status == "on":
        jlcib_null = conditional(le(r, ccm_constants.R_pyr),\
                        (kmax_C_lcib * c2  - kmax_H_chlor * h2  * ccm_constants.f_chlor / (1+ccm_constants.f_chlor))  / (1 + c2 /ccm_constants.Km_CO2 + h2  * (ccm_constants.f_chlor / (1+ccm_constants.f_chlor)) / ccm_constants.Km_HCO3),\
                        (kmax_C_lcib * c2_ - kmax_H_chlor * h2_ * ccm_constants.f_chlor / (1+ccm_constants.f_chlor))  / (1 + c2_/ccm_constants.Km_CO2 + h2_ * (ccm_constants.f_chlor / (1+ccm_constants.f_chlor)) / ccm_constants.Km_HCO3))
        jlcib_new  = conditional(le(r,lcibGeoDict["r_end"]),\
                                conditional(ge(r,lcibGeoDict["r_start"]), jlcib_null, Constant(0.0)),\
                                Constant(0.0))
    else:
        jlcib_new = Constant(0.0)
    jsp = conditional(le(r, ccm_constants.R_pyr),\
                  ccm_constants.v_sp * c2  - ccm_constants.v_sp * pow(10, ccm_constants.pK_eff - ccm_constants.pH_chlor) * h2  * ccm_constants.f_chlor / (1+ccm_constants.f_chlor),\
                  ccm_constants.v_sp * c2_ - ccm_constants.v_sp * pow(10, ccm_constants.pK_eff - ccm_constants.pH_chlor) * h2_ * ccm_constants.f_chlor / (1+ccm_constants.f_chlor))
    return (jlcib_new + jsp)


# compute the co2 to hco3- conversion flux mediated by CAH3
def get_cah3_flux(r, c1, h1, v_cah3, cah3GeoDict):
    r_start = cah3GeoDict["r_start"]
    r_end = cah3GeoDict["r_end"]
    # first-order rate constant of co2 -> hco3- conversion
    kmax_C = v_cah3 / (1 / ccm_constants.Time)
    # first-order rate constant of hco3- -> co2 conversion
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



# compute the weak form of the PDEs describing the reaction-diffusion-transport processes
# involved in the CCM
def gen_weak_form_eqn(u, VDG, k_c, D_c_slow, D_h_slow, R_tub, C_cyt, H_cyt,\
                      tub_trans_type, chlor_trans_type, v_tub, v_chlor, gamma,\
                      stroma_ca_stat, v_lcib, v_cah3,\
                      stroma_barrier, ccmGeoDict, thy_geo, rub_geo,\
                      dx, ds, dS, pyr_marker, tubule_marker, pyr_mem_marker, chlor_mem_marker):

    V = u.function_space()

    # Define trialfunction and testfunction
    du = TrialFunction(V)
    v = TestFunction(V)
    (v0, v1, v1_, v2, v3, v3_) = split(v)

    # Split system functions to access components
    # c=[CO2], h=[HCO3-], 1=tubule, 2=matrix, 2_ = stroma
    (c1, c2, c2_, h1, h2, h2_) = split(u)

    r = Expression("x[0]", degree=1)
    k_c_pyrmem = ccmGeoDict["Starch_Permeability"]["kcPyrMem"] / (
        ccm_constants.Lnth / ccm_constants.Time)
    k_h_pyrmem = ccmGeoDict["Starch_Permeability"]["khPyrMem"] / (
        ccm_constants.Lnth / ccm_constants.Time)

    # ===================================================================== #
    # + Geometric correction
    fv_tub, fs_tub = get_geo_factors(thy_geo, R_tub)
    
    fs_stroma = (fv_tub / (1 - fv_tub)) * fs_tub  # for the stroma
    # ===================================================================== #

    # ===================================================================== #
    # rate of consumption of CO2 by Rubisco
    R_rub, R_rub_ = get_Rubisco_flux(c2, c2_, r, dx, fv_tub,
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
    # bst mediated hco3- flux
    j_active_tub = get_tub_mem_flux(h1, h2, h2_, r, tub_trans_type, v_tub,
                                    gamma, ccmGeoDict["BST_Geometry"])
    # lcia mediated hco3- flux
    j_active_chlor = get_chlor_mem_flux(h2_, H_cyt, chlor_trans_type, v_chlor,
                                        gamma)

    # passive diffusion of hco3- and h2co3 across the thylakoid membranes
    jh_passive = conditional(le(r,ccm_constants.R_pyr),\
                             ccm_constants.k_h2co3 * ((h1 / (1 + ccm_constants.f_tub)) - (h2 /(1 + ccm_constants.f_chlor))) + ccm_constants.k_hco3 * ((h1 * ccm_constants.f_tub / (1 + ccm_constants.f_tub)) - (h2 * ccm_constants.f_chlor / (1 + ccm_constants.f_chlor))),\
                             ccm_constants.k_h2co3 * ((h1 / (1 + ccm_constants.f_tub)) - (h2_ /(1 + ccm_constants.f_chlor))) + ccm_constants.k_hco3 * ((h1 * ccm_constants.f_tub / (1 + ccm_constants.f_tub)) - (h2_ * ccm_constants.f_chlor / (1 + ccm_constants.f_chlor))))

    # total flux of hco3- and h2co3 from tubules to pyrenoid
    jh_12 = conditional(le(r, R_tub), jh_passive - j_active_tub, Constant(0.0))
    # ===================================================================== #

    # ===================================================================== #
    # rate of co2 -> hco3- conversion catalyzed by CAH3
    j_ca = get_cah3_flux(r, c1, h1, v_cah3, ccmGeoDict["CAH3_Geometry"])
    # ===================================================================== #

    # ===================================================================== #
    # rate of sco2 -> hco3- conversion catalyzed by LCIB
    j_lcib = get_lcib_flux(stroma_ca_stat, r, c2, h2, c2_, h2_, v_lcib,
                           ccmGeoDict["LCIB_Geometry"])
    # ===================================================================== #

    # differential diffusion (potentially slowed by thylakoid stacks in the stroma)
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

    # ===================================================================== #
    # reaction-diffusion eqns
    # CO2 in the thylakoid lumen
    f_c1 = (ccm_constants.D_c * r * r * fv_tub * grad(c1)[0] * grad(v0)[0]) \
         + (jc_12 * r * r * fs_tub * fv_tub * v0 ) \
         + (r * r * j_ca * fv_tub * v0)                 # CA equilibration term
    L_c1 = f_c1 * (dx(tubule_marker) + dx(pyr_marker))

    # CO2 in the pyernoid matrix
    f_c2 = (D_c_stroma * r * r * (1.0-fv_tub) * grad(c2)[0] * grad(v1)[0] ) \
         - (jc_12 * r * r * fs_stroma * (1.0-fv_tub) * v1 ) \
         + (R_rub * r * r * (1.0-fv_tub) * v1)\
         + (r * r * j_lcib * (1.0-fv_tub) * v1)
    L_c2 = (f_c2 * dx(pyr_marker)
            ) + k_c_pyrmem * ccm_constants.R_pyr * ccm_constants.R_pyr * (
                1.0 - fv_tub) * (avg(c2) - avg(c2_)) * (
                    avg(v1) - avg(v1_)) * dS(pyr_mem_marker)

    # CO2 in the stroma
    f_c2_ = (D_c_stroma * r * r * (1.0-fv_tub) * grad(c2_)[0] * grad(v1_)[0] ) \
          - (jc_12 * r * r * fs_stroma * (1.0-fv_tub) * v1_ ) \
          + (r * r * j_lcib * (1.0-fv_tub) * v1_)\
          + (R_rub_ * r * r * (1.0-fv_tub) * v1_)
    L_c2_ = (f_c2_ * (dx(tubule_marker) + dx(0))) - (
        v1_ * ccm_constants.R_chlor * ccm_constants.R_chlor * k_c *
        (C_cyt - c2_) * ds(chlor_mem_marker))

    # HCO3- + H2CO3 in the thylakoid lumen
    f_h1 = (ccm_constants.D_h_tub * r * r * fv_tub * grad(h1)[0] * grad(v2)[0]) \
         + (jh_12   * r * r * fv_tub * fs_tub * v2) \
         - (r * r * j_ca * fv_tub * v2)                 # CA equilibration term
    L_h1 = f_h1 * (dx(tubule_marker) + dx(pyr_marker))

    # HCO3- + H2CO3 in the pyrenoid matrix
    f_h2 = (D_h_stroma * r * r * (1.0 - fv_tub) * grad(h2)[0] * grad(v3)[0])\
            - (jh_12 * r * r * fs_stroma * (1.0 - fv_tub) * v3)\
            - (r * r * j_lcib * (1.0-fv_tub) * v3)
    L_h2 = (f_h2 * dx(pyr_marker)
            ) + k_h_pyrmem * ccm_constants.R_pyr * ccm_constants.R_pyr * (
                1.0 - fv_tub) * (avg(h2) - avg(h2_)) * (
                    avg(v3) - avg(v3_)) * dS(pyr_mem_marker)

    # HCO3- + H2CO3 in the stroma
    f_h2_ = (D_h_stroma * r * r * (1.0 - fv_tub) * grad(h2_)[0] * grad(v3_)[0]) \
          - (jh_12 * r * r * fs_stroma * (1.0 - fv_tub) * v3_) \
          - (r * r * j_lcib * (1.0-fv_tub) * v3_)
    L_h2_ = (f_h2_ * (dx(tubule_marker) + dx(0)) ) \
        - (v3_ * ccm_constants.R_chlor * ccm_constants.R_chlor * ccm_constants.k_h2co3 * ((H_cyt / ccm_constants.f_cyt) - (h2_ /(1 + ccm_constants.f_chlor) )) * ds(chlor_mem_marker)) \
        - (v3_ * ccm_constants.R_chlor * ccm_constants.R_chlor * ccm_constants.k_hco3  * ((H_cyt) - (h2_ * ccm_constants.f_chlor / (1 + ccm_constants.f_chlor))) * ds(chlor_mem_marker)) \
        - (v3_ * ccm_constants.R_chlor * ccm_constants.R_chlor * j_active_chlor
           * ds(chlor_mem_marker))
    # ===================================================================== #

    L_tot = L_c1 + L_c2 + L_c2_ + L_h1 + L_h2 + L_h2_
    Jacobian = derivative(L_tot, u, du)

    return L_tot, Jacobian
