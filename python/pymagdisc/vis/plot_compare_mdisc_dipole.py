# Description: Plot magnetodisc vs dipole model.

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import argparse
from pymagdisc.data.load_data import load_model
from pymagdisc.analysis.physics import calc_Pcold, calc_beta, calc_Jphi
from pymagdisc.bfield.magneticField3D import (
    mdiscMagneticField3D,
    dipoleMagneticField3D,
    magneticFieldInterpolator)
from pymagdisc import config
import warnings

warnings.filterwarnings("ignore")


def plot_compare_mdisc_dipole(MDFile="jup_mdisc_kh3e7_rmp90fix.mat"):
    """Example to compare Jupiter's magnetodisc to a dipole.

    Args:
        MDFile (str, optional): Magnetodisc model output file. Defaults to 'jup_mdisc_kh3e7_rmp90fix.mat'.
    """

    MD = load_model(f"{config.PATH_TO_DATA}{MDFile}")

    # equatorial radius in m
    Re = MD["planet"]["r0"]

    # magnetic equator field strength in T
    B0 = MD["planet"]["B0"]

    # corrected to match dipole value in the magnetodisc file
    B0 = MD["planet"]["B0"] * MD["v2d"]["Bthdip"][MD["dims"]["imu0"] - 1, 0]
    Md = [0, 0, B0 * Re ** 3]  # Magnetic moment (z-aligned)
    Rd = [0, 0, 0]  # Centered moment (no offset)

    p = 0
    t = np.linspace(-np.pi / 2, np.pi / 2, 100)

    # unit half circle
    x = np.cos(t) * np.cos(p)
    y = np.cos(t) * np.sin(p)
    z = np.sin(t)
    interpolator = magneticFieldInterpolator(MD)
    B, gradB = mdiscMagneticField3D([x, y, z], MD, interpolator=interpolator)
    Bd, gradBd = dipoleMagneticField3D(Md, Rd, [Re * x, Re * y, Re * z])

    xc = x
    yc = y
    zc = z

    # plot B vs theta
    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    ax.plot(t, B[0], "r-", label="$B_x$")
    ax.plot(t, Bd[0], "r--", label="$B_{xd}$")
    ax.plot(t, B[2], "k-", label="$B_z$")
    ax.plot(t, Bd[2], "k--", label="$B_{zd}$")
    ax.set_xlabel(r"$\theta$")
    ax.set_ylabel("B [T]")
    ax.legend()

    fig, ax = plt.subplots(2, 1, figsize=(10, 10), sharex=False, sharey=False)

    for L in [5, 10, 20]:
        r = L * np.cos(t) ** 2
        ir = np.argwhere(r >= 1).reshape(-1)
        x = r[ir] * np.cos(t[ir]) * np.cos(p)
        y = r[ir] * np.cos(t[ir]) * np.sin(p)
        z = r[ir] * np.sin(t[ir])

        # plot mdisc field
        B, gradB = mdiscMagneticField3D([x, y, z], MD=MD, interpolator=interpolator, compute_gradB=True)
        b = np.sqrt(B[0] ** 2 + B[1] ** 2 + B[2] ** 2)
        db = np.sqrt(gradB[0] ** 2 + gradB[1] ** 2 + gradB[2] ** 2)
        ax[0].quiver(
            x,
            z,
            B[0] / b,
            B[2] / b,
            scale=40,
            width=0.001,
            headwidth=10,
            color="r",
            label=r"$B$",
        )  # N.B. a smaller scale parameter makes the arrow longer
        ax[1].quiver(
            x,
            z,
            gradB[0] / db,
            gradB[2] / db,
            scale=40,
            width=0.001,
            headwidth=10,
            color="r",
            label=r"$\nabla B$",
        )

        # plot dipole field
        Bd, gradBd = dipoleMagneticField3D(Md, Rd, [Re * x, Re * y, Re * z])
        bd = np.sqrt(Bd[0] ** 2 + Bd[1] ** 2 + Bd[2] ** 2)
        dbd = np.sqrt(gradBd[0] ** 2 + gradBd[1] ** 2 + gradBd[2] ** 2)
        ax[0].quiver(
            x,
            z,
            Bd[0] / bd,
            Bd[2] / bd,
            scale=40,
            width=0.001,
            headwidth=10,
            color="k",
            label=r"$B_d$",
        )
        ax[1].quiver(
            x,
            z,
            gradBd[0] / dbd,
            gradBd[2] / dbd,
            scale=40,
            width=0.001,
            headwidth=10,
            color="k",
            label=r"$\nabla B_d$",
        )

    ax[1].set_xlabel("x [$R_J$]")
    [axi.set_ylabel("z [$R_J$]") for axi in ax]
    [axi.axis("equal") for axi in ax]
    ax[0].legend([r"$B$", r"$B_d$"])
    ax[1].legend([r"$\nabla B$", r"$\nabla B_d$"])
    [axi.plot(xc, zc, "k-") for axi in ax]

    plt.show()


if __name__ == "__main__":
    plot_compare_mdisc_dipole(MDFile="jup_mdisc_kh3e7_rmp90fix.mat")
