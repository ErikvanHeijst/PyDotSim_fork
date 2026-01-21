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

    def calculate_derived_properties(self):
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


if __name__ == "__main__":
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

    indium_arsenide.calculate_derived_properties()
    gallium_arsenide.calculate_derived_properties()
