import numpy as np
from pymagdisc.tracer.si_constants import e, eV, c, mp, me, amu


class Particle:
    """Particle used in Boris particle tracer."""

    def __init__(self, partype, Ep):
        if partype == "p":
            self.parname = "proton"
            self.mp = mp  # proton mass kg
            self.qp = e
        elif partype == "e":
            self.parname = "electron"
            self.mp = me  # electron mass kg
            self.qp = -e
        elif partype == "O+":
            self.parname = "O+"
            self.mp = 15.999 * amu
            self.qp = e
        elif partype == "O++":
            self.parname = "O++"
            self.mp = 15.999 * amu
            self.qp = 2 * e
        elif partype == "S+":
            self.parname = "S+"
            self.mp = 32.07 * amu
            self.qp = e
        elif partype == "S++":
            self.parname = "S++"
            self.mp = 32.07 * amu
            self.qp = 2 * e
        elif partype == "S+++":
            self.parname = "S+++"
            self.mp = 32.07 * amu
            self.qp = 3 * e

        # Basic velocity and energy calculations
        self.EMeV = Ep * 1e6 * eV  # particle energy in MeV
        self.vpc = np.sqrt(2 / self.mp * self.EMeV)  # classic velocity
        print(
            f"classic v = {self.vpc:.4g} m/s; E={0.5 * self.mp * self.vpc**2 / (1e6 * eV):.4g} MeV"
        )

        self.vpr = (
            c
            * np.sqrt(self.EMeV * (self.EMeV + 2 * self.mp * c**2))
            / (self.EMeV + self.mp * c**2)
        )  # relativistic velocity

        print(
            f"relativ v = {self.vpr:.4g} m/s; E={(1 / np.sqrt(1 - self.vpr**2 / c**2) - 1) * self.mp * c**2 / (1e6 * eV):.4g} MeV"
        )
        self.vp = self.vpr
        # velocity squared
        self.V2 = self.vp**2

        # relativistic factor
        self.gamma = 1 / np.sqrt(1 - self.V2 / c**2)

        # particle charge-to-mass ratio
        self.qOverM = self.qp / self.mp
