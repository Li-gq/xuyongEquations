import math

# Values retrieved from physics.nist.gov/constants
BOLTZMANN_SI = 1.38064852e-23
FINE_STRUCT = 7.2973525664e-3
AVOGADRO = 6.022140857e23  # [1/mol]
Clight_SI = 299792458.  # [m/s]
a0_SI = 5.2917721067e-11  # [m]
me_SI = 9.10938356e-31  # [kg]
qe_SI = 1.6021766208e-19  # Electron Charge [C]

Clight = 1. / FINE_STRUCT
EPSILON0 = 1. / (4. * math.pi)
MU0 = 1. / (EPSILON0 * Clight**2)
MUB = 1. / 2.

Meter = 1. / a0_SI
Kilogram = 1. / me_SI
Coulomb = 1. / qe_SI
Second = Clight_SI / Clight * Meter

Newton = Kilogram * Meter / Second**2
Joule = Newton * Meter
Volt = Joule / Coulomb
Ampere = Coulomb / Second
Ohm = Volt / Ampere
Siemens = Ampere / Volt

hbar_SI = 1 / (Joule * Second)
Angstrom = Meter * 1e-10
eV = Volt
BOLTZMANN = BOLTZMANN_SI * Joule
Ha = 1.0

# 实验参数
params = {
    'probe':{
    'a': 6,
    'EM': 1.0,
    'fermi': 0.37*eV,
    'workfunction': 4.5*eV
},
'superconductor':{
    'a': 7.3*Angstrom,
    'gap': 1.8e-3*eV,
    'tc': 14.5,
    'EM': 12,  # 表面电子有效质量
    'Dirac_fermi': 0.2e-3*eV,
    'Dirac_workfunction': 4.0*eV
}}
# print(params['superconductor']['a']*a0_SI)
# print(Joule)