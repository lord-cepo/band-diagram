import material

"""
Description
-----------
Datasheet of some of the most used semiconductor and metals. A list of all
semiconductors and all metals is defined.

Semiconductors
--------------
Semiconductors are called by their chemical formula, as a function of doping 
type and doping density. Ternary compounds (e.g. Al(x) Ga(1-x) As) are functions
with the additional parameter x.
Example: Silicon = Si('n', 1e15)

Metals
------
Metals are directly initializated, without any user-defined parameter.
Example: Gold = Au
"""

# https://www.ioffe.ru/SVA/NSM/Semicond/

def Si(doping_type=None, doping=0):
    SI_EM_LONG = 0.98
    SI_EM_TRAN = 0.19
    SI_HM_HEAV = 0.49
    SI_HM_LIGH = 0.16
    SI_N_V = 6
    return material.semiconductor(
        Eg=1.12, 
        electron_affinity=4.05, 
        epsilon=11.7, 
        doping_type=doping_type,
        doping=doping,
        effective_m_e=(SI_N_V**2*SI_EM_LONG*SI_EM_TRAN**2)**(1/3),
        effective_m_h=(SI_HM_HEAV**1.5+SI_HM_LIGH**1.5)**(2/3)
    )

def GaAs(doping_type=None, doping=0):
    GAAS_EM_LONG = 0.063 
    GAAS_EM_TRAN = 0.063 # seems isotropic
    GAAS_HM_HEAV = 0.51
    GAAS_HM_LIGH = 0.082
    GAAS_N_V = 1 # direct gap, only one valley at Gamma
    return material.semiconductor(
        Eg=1.424, 
        electron_affinity=4.07, 
        epsilon=12.9, 
        doping_type=doping_type,
        doping=doping,
        effective_m_e=(GAAS_N_V**2*GAAS_EM_LONG*GAAS_EM_TRAN**2)**(1/3),
        effective_m_h=(GAAS_HM_HEAV**1.5+GAAS_HM_LIGH**1.5)**(2/3)
    )   

def Ge(doping_type=None, doping=0):
    GE_EM_LONG = 1.6
    GE_EM_TRAN = 0.08
    GE_HM_HEAV = 0.33
    GE_HM_LIGH = 0.043
    GE_N_V = 4 
    return material.semiconductor(
        Eg=0.661, 
        electron_affinity=4.0, 
        epsilon=16.2, 
        doping_type=doping_type,
        doping=doping,
        effective_m_e=(GE_N_V**2*GE_EM_LONG*GE_EM_TRAN**2)**(1/3),
        effective_m_h=(GE_HM_HEAV**1.5+GE_HM_LIGH**1.5)**(2/3)
    )   

def GaP(doping_type=None, doping=0):
    GAP_EM_LONG = 1.12
    GAP_EM_TRAN = 0.22
    GAP_HM_HEAV = 0.79
    GAP_HM_LIGH = 0.14
    GAP_N_V = 3
    return material.semiconductor(
        Eg=2.26, 
        electron_affinity=3.8, 
        epsilon=11.1, 
        doping_type=doping_type,
        doping=doping,
        effective_m_e=(GAP_N_V**2*GAP_EM_LONG*GAP_EM_TRAN**2)**(1/3),
        effective_m_h=(GAP_HM_HEAV**1.5+GAP_HM_LIGH**1.5)**(2/3)
    ) 

def InAs(doping_type=None, doping=0):
    INAS_EM_LONG = 0.023
    INAS_EM_TRAN = 0.023
    INAS_HM_HEAV = 0.41
    INAS_HM_LIGH = 0.026
    INAS_N_V = 1
    return material.semiconductor(
        Eg=0.354, 
        electron_affinity=4.9, 
        epsilon=12.3, 
        doping_type=doping_type,
        doping=doping,
        effective_m_e=(INAS_N_V**2*INAS_EM_LONG*INAS_EM_TRAN**2)**(1/3),
        effective_m_h=(INAS_HM_HEAV**1.5+INAS_HM_LIGH**1.5)**(2/3)
    ) 

def GaSb(doping_type=None, doping=0):
    GASB_EM_LONG = 0.041
    GASB_EM_TRAN = 0.041
    GASB_HM_HEAV = 0.4
    GASB_HM_LIGH = 0.05
    GASB_N_V = 1
    return material.semiconductor(
        Eg=0.726, 
        electron_affinity=4.06, 
        epsilon=15.7, 
        doping_type=doping_type,
        doping=doping,
        effective_m_e=(GASB_N_V**2*GASB_EM_LONG*GASB_EM_TRAN**2)**(1/3),
        effective_m_h=(GASB_HM_HEAV**1.5+GASB_HM_LIGH**1.5)**(2/3)
    ) 

def InSb(doping_type=None, doping=0):
    INSB_EM_LONG = 0.014
    INSB_EM_TRAN = 0.014
    INSB_HM_HEAV = 0.43
    INSB_HM_LIGH = 0.015
    INSB_N_V = 1
    return material.semiconductor(
        Eg=0.17, 
        electron_affinity=4.59, 
        epsilon=16.7, 
        doping_type=doping_type,
        doping=doping,
        effective_m_e=(INSB_N_V**2*INSB_EM_LONG*INSB_EM_TRAN**2)**(1/3),
        effective_m_h=(INSB_HM_HEAV**1.5+INSB_HM_LIGH**1.5)**(2/3)
    ) 

def InP(doping_type=None, doping=0):
    INP_EM_LONG = 0.08
    INP_EM_TRAN = 0.08
    INP_HM_HEAV = 0.6
    INP_HM_LIGH = 0.089
    INP_N_V = 1
    return material.semiconductor(
        Eg=1.344, 
        electron_affinity=4.38, 
        epsilon=12.5, 
        doping_type=doping_type,
        doping=doping,
        effective_m_e=(INP_N_V**2*INP_EM_LONG*INP_EM_TRAN**2)**(1/3),
        effective_m_h=(INP_HM_HEAV**1.5+INP_HM_LIGH**1.5)**(2/3)
    ) 

def AlN(doping_type=None, doping=0):
    return material.semiconductor(
        Eg=6.1, 
        electron_affinity=0.6, 
        epsilon=9.14, 
        doping_type=doping_type,
        doping=doping,
        effective_m_e=0.4,
        effective_m_h=7.26
    ) 

def InN(doping_type=None, doping=0):
    return material.semiconductor(
        Eg=1.97, # It has a narrower fundamental BG !!!!
        electron_affinity=4.7, 
        epsilon=8.4, 
        doping_type=doping_type,
        doping=doping,
        effective_m_e=0.12,
        effective_m_h=1.65
    ) 

def GaN(doping_type=None, doping=0):
    return material.semiconductor(
        Eg=3.2, 
        electron_affinity=4.1, 
        epsilon=8.9, 
        doping_type=doping_type,
        doping=doping,
        effective_m_e=0.20,
        effective_m_h=1.5
    ) 

# Intended as Al(x) Ga(1-x) As
def AlGaAs(x, doping_type=None, doping=0):
    return material.semiconductor(
        Eg=1.424+1.247*x if x<0.45 else 1.9+0.125*x+0.143*x**2,
        electron_affinity=4.07-1.1*x if x<0.45 else 3.64-0.14*x,
        epsilon=12.9-2.84*x,
        doping_type=doping_type,
        doping=doping,
        effective_m_e=(0.063+0.083*x)**1.5 if x<0.41 else (0.85-0.14*x)**1.5,
        effective_m_h=(0.51+0.25*x)**1.5 if x<0.41 else (0.85-0.14*x)**1.5
    )

"""
Some thought on effective masses. We want density of states masses, given (for
electrons) as effective_m_e = (n_v^2 mx my mz)^(1/3). Silicon has (as an 
example) two transversal masses and one longitudinal mass and six equivalent
direcions n_v (six cigars). Regarding hole masses, I have not found anything on
the web, but since DOS scales with m^3/2 and since most of the semiconductors
have p level splitting (and therefore light and heavy holes) I assumed that the
formula (for isotropic valence band (???)) is effective_m_h^3/2 = m_heavy^3/2 +
m_light^3/2. Maybe split-off band does something to DOS at different energy gap,
but it's not that important (we are tolerating worse approximations)
"""

list_of_semiconductors = [
    Si, GaAs, Ge, GaP, InAs, GaSb, InSb, InP, AlN, InN, GaN 
]

# https://en.wikipedia.org/wiki/Work_function#Work_functions_of_elements

# TODO: find electron density at Fermi energy
Ag = material.metal(work_function=4.26) # up to 4.74
Al = material.metal(work_function=4.06) # up tp 4.26
As = material.metal(work_function=3.75)
Au = material.metal(work_function=5.10) # up to 5.47
B  = material.metal(work_function=4.45)
C  = material.metal(work_function=5)
Ca = material.metal(work_function=2.87)
Cd = material.metal(work_function=4.08)
Cu = material.metal(work_function=4.53) # up to 5.10
Fe = material.metal(work_function=4.67) # up to 4.81
Ga = material.metal(work_function=4.32)
In = material.metal(work_function=4.09)
K  = material.metal(work_function=2.29)
Li = material.metal(work_function=2.9)
Mg = material.metal(work_function=3.66)
Mn = material.metal(work_function=4.1)
Na = material.metal(work_function=2.36)
Pb = material.metal(work_function=4.25)
Pd = material.metal(work_function=5.22) # up to 5.60
Pt = material.metal(work_function=5.12) # up to 5.93
Sb = material.metal(work_function=4.55) # up to 4.70
Se = material.metal(work_function=5.9)
Sn = material.metal(work_function=4.42)
Te = material.metal(work_function=4.95)
Ti = material.metal(work_function=4.33)
Zn = material.metal(work_function=3.63) # up to 4.9

list_of_metals = [
    Ag, Al, As, Au, B , C , Ca, Cd, Cu, Fe, Ga, In, K , 
    Li, Mg, Mn, Na, Pb, Pd, Pt, Sb, Se, Sn, Te, Ti, Zn
]