import numpy as np
from pymagdisc.tracer.si_constants import e, eV, c, mp, me


class Particle:
    """Particle used in Boris particle tracer."""

    def __init__(self, partype, Ep):
        if partype == "p":
            self.mp = mp  # proton mass kg
            self.qp = e
        elif partype == "e":
            self.mp = me  # electron mass kg
            self.qp = -e
        elif partype == "O+":
            self.mp = 16 * mp
            self.qp = e
        elif partype == "O++":
            self.mp = 16 * mp
            self.qp = 2 * e
        elif partype == "S+":
            self.mp = 32 * mp
            self.qp = e
        elif partype == "S++":
            self.mp = 32 * mp
            self.qp = 2 * e
        elif partype == "S+++":
            self.mp = 32 * mp
            self.qp = 3 * e

        # Basic velocity and energy calculations
        self.EMeV = Ep * 1e6 * eV  # particle energy in MeV
        self.vp = np.sqrt(2 / self.mp * self.EMeV)  # classic velocity
        print(
            f"classic v = {self.vp:.4e} m/s; E={0.5*self.mp*self.vp**2/(1e6*eV):.4e} MeV \n"
        )

        self.vp = (
            c
            * np.sqrt(self.EMeV * (self.EMeV + 2 * self.mp * c ** 2))
            / (self.EMeV + self.mp * c ** 2)
        )  # relativistic velocity

        print(
            f"relativ v = {self.vp:.4e} m/s; E={(1/np.sqrt(1-self.vp**2/c**2)-1)*self.mp*c**2/(1e6*eV):.4e} MeV \n"
        )
        # velocity squared
        self.V2 = self.vp ** 2

        # relativistic factor
        self.gamma = 1 / np.sqrt(1 - self.V2 / c ** 2)

        # particle charge-to-mass ratio
        self.qOverM = self.qp / self.mp

        self.rho = 0
