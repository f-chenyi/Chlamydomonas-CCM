"""
This file contains all the biological constants of Chlamy CCM.

--------------------------------------------------------------
"""

# ============================================================= #
# + Define relevant length scale, time scale, and conc. scale

Lnth = 3.14 * pow(10, -6)                                           # Typical length scale; same as the chloroplast radius
Time = Lnth * Lnth / pow(10, -9)                                    # Characteristic time for small molecules to diffuse over the
                                                                    # chloroplast assuming diffusion constant ~ 1E-9 m^2/s
Conc = 0.001                                                        # Typical concentration 1 mM
# ============================================================= #

# ============================================================= #
# + Passive transport: permeability of different species

k_h2co3 = 3 * pow(10, -5) / (Lnth / Time)                           # permeability of membrane to H2CO3, m/s
k_hco3  = 5 * pow(10, -8) / (Lnth / Time)                           # permeability of membrane to HCO3-, m/s
# ============================================================= #

# ============================================================= #
# + pH and pKa [o]
pKa1     = 3.4                                                      # 1st pKa of H2CO3
pH_tub   = 6.0                                                      # tubule pH
pH_cyt   = 7.1                                                      # cytosolic pH
pH_chlor = 8.0                                                      # stromal pH
pK_eff   = 6.1                                                      # effective pKa of CO2 <--> HCO3- conversion

# + pH dependent partition factors
# The fraction of HCO3- is f / (1+f), and
# the fraction of H2CO3 is 1 / (1+f) in a given compartment
f_cyt   = pow(10, -pKa1 + pH_cyt)
f_chlor = pow(10, -pKa1 + pH_chlor)
f_tub   = pow(10, -pKa1 + pH_tub)
# ============================================================= #

# ============================================================= #
# + Diffusion constant
D_c     = 1.88 * pow(10, -9) / (Lnth * Lnth / Time)                 # diffusion constant of CO2 in water, m^2/s
D_h     = 1.15 * pow(10, -9) / (Lnth * Lnth / Time)                 # diffusion constant of HCO3- in water, m^2/s
D_h2co3 = 1.15 * pow(10, -9) / (Lnth * Lnth / Time)                 # diffusion constant of H2CO3 in water, m^2/s (UPDATE)

# for the sum of HCO3- + H2CO3 in each compartment
D_h_chlor = f_chlor / (1 + f_chlor) * D_h + 1 / (1 + f_chlor) * D_h2co3
D_h_tub   = f_tub   / (1 + f_tub)   * D_h + 1 / (1 + f_tub)   * D_h2co3

# ============================================================= #

# ============================================================= #
# + Geometry factors
R_chlor = 3.14 * pow(10, -6) / Lnth                                 # radius of Chlamy chloroplast
R_pyr   = 0.3 * R_chlor                                             # radius of Chlamy pyrenoid

N_tub = 40                                                          # number of thylkaoid tubules
rin   = 0.4 * pow(10, -6)  / Lnth                                   # radius of the tubule meshwork
a_tub = 0.05 * pow(10, -6) / Lnth                                   # cylindrical radius of the thylakoid tubules
fv_in = N_tub / 4 * (a_tub / rin)**2                                # volume fraction

# for plant thylakoid geometry (!)
# (optional, assuming a meshwork structure of thylakoids
#  throughout the chloroplast)
fv_plant = 0.35
a_plant = 0.25 * pow(10, -6) / Lnth
# ============================================================= #

# ============================================================= #
# + Rubsico reaction kinetics
C_rub  = 0.005                                                      # conc. of Rubisco active sites, M
kcat   = 3.0                                                        # kcat(turnover number) of Rubisco, s^-1
Km_c   = 3.0 * pow(10, -5) / Conc                                   # concentration of CO2 that half-maximizes rate of Rubisco
Km_o   = 15.0 * pow(10, -5) / Conc                                  # concentration of O2 that half-maximizes rate of Rubisco
O      = 0.00023 / Conc                                             # concentration of O2 in the pyrenoid

Vmax   = (kcat * C_rub) / (Conc / Time)                             # maximum velocity of Rubisco carboxylation
Km_eff = Km_c * (1 + O / Km_o)                                      # effective Km of Rubisco for CO2
# ============================================================= #

# ============================================================= #
# + Carbonic anhydrase kinetics
v_sp       = 0.036 / (1 / Time)                                     # velocity of spontaneous interconversion
#kmax_C     = 1E4 / (1 / Time)
#kmax_H_tub = kmax_C * pow(10, pK_eff - pH_tub)
#kmax_H_chlor  = kmax_C_lcib * pow(10, pK_eff - pH_chlor)

Km_CO2     = 0.005 / Conc                                           # estimates of Km for carbonic anhydrases
Km_HCO3    = 0.005 / Conc
# ============================================================= #

#Km_active_chlor = 0.005 / Conc
#Km_active_tub   = 0.005 / Conc
