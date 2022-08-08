# Description: Plot contours of magnetic potential for mdisc and dipole model.

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import argparse
from pymagdisc.data.load_data import load_model
from pymagdisc.analysis.physics import calc_Pcold, calc_beta, calc_Jphi
from pymagdisc.bfield.magneticField3D import (
    mdiscMagneticField3D,
    dipoleMagneticField3D,
    MDiscField,
    _calc_posn,
)
from pymagdisc import config
import warnings
from scipy.interpolate import RectBivariateSpline as interp2

warnings.filterwarnings("ignore")


def plot_compare_contours_alpha(
    radDist, colatDeg, MDFile="jup_mdisc_kh3e7_rmp90fix.mat"
):
    """Plot the contours for magnetic potential given colatitude (Deggrees)
    for mdisc and dipole models.

    Args:
        radDist (int): Cylindrical radial distance [RJ].
        colatDeg (array): Colatitude [Degrees].
        MDFile (str, optional): Magnetodisc model output file. Defaults to 'jup_mdisc_kh3e7_rmp90fix.mat'.

    """
    MD = load_model(f"{config.PATH_TO_DATA}{MDFile}")

    # Ionosphere footpoint
    mu = np.cos(np.deg2rad(colatDeg))

    # Calculate disc and dipole values:
    mdiscVal = interp2(
        MD["c2d"]["r"][0, :], MD["c2d"]["mu"][:, 0], MD["v2d"]["alpha"].T
    ).ev(radDist, mu)
    dipVal = interp2(
        MD["c2d"]["r"][0, :], MD["c2d"]["mu"][:, 0], MD["v2d"]["alphadip"].T
    ).ev(radDist, mu)

    # Plot the contours:
    c1 = plt.contour(
        MD["c2d"]["rho"],
        MD["c2d"]["z"],
        MD["v2d"]["alpha"],
        levels=np.unique(sorted(mdiscVal)),
        linestyles="solid",
        colors="k",
    )
    c2 = plt.contour(
        MD["c2d"]["rho"],
        MD["c2d"]["z"],
        MD["v2d"]["alphadip"],
        levels=np.unique(sorted(dipVal)),
        linestyles="dashed",
        colors="k",
    )
    h1, _ = c1.legend_elements()
    h2, _ = c2.legend_elements()
    plt.legend([h1[0], h2[0]], ["MDisc", "Dipole"])

    pc = plt.pcolormesh(
        MD["c2d"]["rho"],
        MD["c2d"]["z"],
        np.log10(MD["v2d"]["alpha"]),
        shading="gouraud",
        vmin=-4,
        vmax=0,
        cmap="jet",
    )
    plt.colorbar(label="Magnetic potential")
    plt.xlabel(r"$\rho$ [$R_J$]")
    plt.ylabel(r"z [$R_J$]")
    plt.xlim([-2, 32])
    plt.ylim([-17, 17])
    plt.show()


if __name__ == "__main__":
    plot_compare_contours_alpha(
        radDist=1, colatDeg=[6, 12, 17, 20, 25], MDFile="jup_mdisc_kh3e7_rmp90fix.mat"
    )
