"""

Constants for SimQDStates

by Christian Heyn

"""
# %%
import math
import sys

from dataclasses import dataclass
from scipy.constants import (
    pi,
    h,
    hbar,
    electron_mass,
    epsilon_0,
    elementary_charge,
)


@dataclass
class Material:
    relative_electron_mass: float
    relative_hole_mass: float
    Kane_energy: float
    refractive_index: float

    band_gap_T_0: float
    Varshni_a: float
    Varshni_b: float

    temperature: float = 4

    def __post_init__(self):
        # Note mass here means effective mass
        self.electron_mass = self.relative_electron_mass * electron_mass
        self.hole_mass = self.relative_hole_mass * electron_mass
        self.permittivity = self.refractive_index**2 * epsilon_0

        self.exciton_mass = self._exciton_mass(self.electron_mass, self.hole_mass)
        self.bohr_radius = self._bohr_radius(self.permittivity, self.exciton_mass)
        self.rydberg_constant = self._rydberg_constant(
            self.permittivity, self.exciton_mass
        )

        self.band_gap = self._band_gap(
            self.band_gap_T_0, self.Varshni_a, self.Varshni_b, self.temperature
        )

    def _exciton_mass(self, electron_mass, hole_mass):
        exciton_mass = 1 / ((1 / electron_mass) + (1 / hole_mass))

        return exciton_mass

    def _bohr_radius(self, permittivity, exciton_mass):
        bohr_radius = (4 * pi * permittivity * hbar**2) / (
            exciton_mass * elementary_charge**2
        )

        return bohr_radius

    def _rydberg_constant(self, permittivity, exciton_mass):
        rydberg_constant = (exciton_mass * elementary_charge**4) / (
            2 * (2 * permittivity * h) ** 2
        )

        return rydberg_constant

    def _band_gap(self, band_gap_T_0, a, b, T):
        band_gap = band_gap_T_0 - ((a * T**2) / (T + b))

        return band_gap


@dataclass
class Alloy:
    host: Material
    guest: Material
    x: float  # composition
    C: float  # bowing parameter
    name: str = "InGaAs"

    def _bandgap(self):
        alloy_bandgap = (
            (1 - self.x) * self.host.band_gap
            + self.x * self.guest.band_gap
            - self.x * (1 - self.x) * self.C
        )

        return alloy_bandgap

    def _effective_mass(self):
        if self.name == "AlGaAs":
            alloy_electron_mass = 0.067 + 0.057 * self.x
            alloy_hole_mass = 0.51 + 0.25 * self.x
        if self.name == "InGaAs":
            alloy_electron_mass = (
                0.024 + 0.035 * (1 - self.x) + 0.008 * (1 - self.x) ** 2
            )
            alloy_hole_mass = 0.53  # for 30% Indium
        #else:
        #    raise ValueError("No such alloy name is recognized")

        return alloy_electron_mass, alloy_hole_mass

    def _bandgap_discontuity(self):
        if self.name == "AlGaAs":
            dE = self._bandgap() - self.host.band_gap
            f = 63 / 37
            dEc = dE * f / (1 + f)
            dEv = dE - dEc
            return dEc, dEv
        if self.name == "InGaAs":
            dE = self._bandgap() - self.host.band_gap
            f = 60 / 40  # GUESS from Yu et al
            dEc = dE * f / (1 + f)
            dEv = dE - dEc
            return dEc, dEv
        #else:
        #    raise ValueError("No such alloy name is recognized")


# #####

gallium_arsenide = Material(
    relative_electron_mass=0.067,
    relative_hole_mass=0.51,
    Kane_energy=28.8,
    refractive_index=3.4,
    band_gap_T_0=1.519,
    Varshni_a=0.5405e-3,
    Varshni_b=204,
    temperature=300,
)

aluminium_arsenide = Material(
    relative_electron_mass=0.15,
    relative_hole_mass=0.76,
    Kane_energy=21.1,
    refractive_index=2.9,
    band_gap_T_0=3.099,
    Varshni_a=0.885e-3,
    Varshni_b=503,
    temperature=300,
)

indium_arsenide = Material(
    relative_electron_mass=0.026,
    relative_hole_mass=0.41,
    Kane_energy=21.5,
    refractive_index=3.45,
    band_gap_T_0=0.417,
    Varshni_a=0.276e-3,
    Varshni_b=93,
    temperature=300,
)

#indium_arsenide.calculate_derived_properties()
#gallium_arsenide.calculate_derived_properties()

indium_gallium_arsenide = Alloy(
    host=indium_arsenide, guest=gallium_arsenide, x=0.3, C=0.477, name="InGaAs"
)

aluminium_gallium_arsenide = Alloy(
    host=gallium_arsenide,
    guest=indium_arsenide,
    x=0.34,
    C=-0.127 + 1.310 * 0.34,
    name="InGaAs",
)

####

h = 6.62606876e-34;             # Planck's constant [J*s]
hbar = h/(2*math.pi);           # Planck's constant reduced [J*s]
epsilon0 = 8.854187817e-12;     # Vacuum permittivity [F/m] = [A*s/(V*m)]
me0 = 9.10938188e-31;           # Electron mass [kg]
c = 299792458;                  # speed of light in m/s
e = 1.60217662e-19;             # elementary charge in coulombs
JeV = 1/1.60217646e-19;         # Joule in eV

# GaAs QDs related constants
me = 0.067
me_GaAs = me*me0                                        # GaAs effective electron mass [kg]
mh = 0.51
mhh_GaAs = mh*me0                                       # GaAs effective heavy hole mass [kg]
epsilon_GaAs = 13.1*epsilon0;                           # GaAs permittivity
m_Ex = 1/(1/me_GaAs + 1/mhh_GaAs);                      # Exciton effective mass
r_Ex_GaAs = 4*math.pi*epsilon_GaAs*hbar**2/(m_Ex*e**2); # GaAs exciton Bohr radius in [m]
E_Ry_GaAs = m_Ex*e**4/(2*(2*epsilon_GaAs*h)**2);        # Rydberg constant in GaAs bulk in [J]
E_p = 28.8;                                             # Kane energy in eV
n = 3.6;                                                # refractive index

def GaAs_bandgap(T): # Varshni
    # Vurgaftman, Thurmond
    EGaAs0 = 1.519; a = 5.405e-4; b = 204;
    EGaAs = EGaAs0 - a*T*T/(T+b);
    #gallium_arsenide.calculate_derived_properties()
    #EGaAs = indium_arsenide.band_gap
    return EGaAs

def AlAs_bandgap(T): # Varshni
    # Vurgaftman
    EGaAs0 = 3.099; a = 8.85e-4; b = 530;
    EGaAs = EGaAs0 - a*T*T/(T+b);
    #aluminium_arsenide.calculate_derived_properties()
    #EGaAs = indium_arsenide.band_gap
    return EGaAs

def AlGaAs_bandgap(T, x): # Varshni
    # Vurgaftman
    C = -0.127 + 1.310*x
    EGaAs = GaAs_bandgap(T)
    EAlAs = AlAs_bandgap(T)
    EAlGaAs = EGaAs + (EAlAs-EGaAs-C)*x+C*x*x

    #EAlGaAs = indium_gallium_arsenide._bandgap()
    return EAlGaAs

def AlGaAs_effMass(x):
    me = 0.067+0.057*x
    mhh = 0.51+0.25*x

    #me, mhh = indium_gallium_arsenide._effective_mass()
    return round(me,4), round(mhh,4)

# GaAs - AlGaAs valence band edge discontinuity
# from Adachi, Properties of semiconductor alloys, 2009: dEc/dEv = 63/37 = f
# dEc+dEv = dE = EAlGaAs-EGaAs
# dEc = dEv*f = (dE-dEc)*f = dE*f-dEc*f
# dEc+dE*f = dEc*(1+f) = dE*f
# dEc = dE*f/(1+f)
def bandgapDiscontinuity(T, x):
    EGaAs = GaAs_bandgap(T)
    EAlGaAs = AlGaAs_bandgap(T, x)
    dE = EAlGaAs-EGaAs
    f = 63/37
    dEc = dE*f/(1+f)
    dEv = dE-dEc

    #dEc, dEv = indium_gallium_arsenide._bandgap_discontuity()
    return dEc, dEv

def QDlifetime(overlapp, E_PL):
    tau = (1/2) * 3 * h**2 * c**3 * epsilon0 * me0 / (math.pi * n * e**2 * E_p/JeV * E_PL/JeV * overlapp)
    return tau
