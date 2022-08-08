# Description: Plot magnetodisc model.

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import argparse
from pymagdisc.data.load_data import load_model
from pymagdisc.analysis.physics import calc_Pcold, calc_beta, calc_Jphi
from pymagdisc import config
import warnings

warnings.filterwarnings("ignore")


def plot_2d_cmp_exp_Pcold(MD_CMP, MD_EXP):
    """Plot cold plasma pressure color map overlayed with B-field configuration for compressed and expanded MD.

    Args:
        MD_CMP (dict): Compressed magnetodisc model.
        MD_EXP (dict): Expanded magnetodisc model.

    """

    # Pressure calc
    pcold2dcmp = calc_Pcold(MD_CMP)
    pcold2dexp = calc_Pcold(MD_EXP)

    betacmp = calc_beta(MD_CMP)
    betaexp = calc_beta(MD_EXP)

    ## Plot
    fig, ax = plt.subplots(1, 2, sharey=True, figsize=(8, 5), constrained_layout=True)

    # Compressed MD
    beta_level = 0.2
    pc = ax[0].pcolormesh(
        MD_CMP["c2d"]["rho"],
        MD_CMP["c2d"]["z"],
        np.log10(pcold2dcmp),
        vmin=-11,
        vmax=-7,
        cmap="jet",
    )
    ax[0].contour(
        MD_CMP["c2d"]["rho"],
        MD_CMP["c2d"]["z"],
        np.log10(MD_CMP["v2d"]["alpha"]),
        levels=[-3, -2] + list(np.arange(-1.5, 0, 0.1)),
        linestyles="solid",
        colors="m",
    )  # specify levels to plot contour lines at specific z function values
    ax[0].contour(
        MD_CMP["c2d"]["rho"],
        MD_CMP["c2d"]["z"],
        betacmp,
        levels=[beta_level],
        linestyles="solid",
        colors="w",
    )
    ax[0].set_title(rf"COMPRESSED: $\beta$={beta_level}")
    ax[0].set_ylabel(r"$Z/R_J$")
    ax[0].set_xlabel(r"$\rho_{CYL}/R_J$")

    # Expanded MD
    pc = ax[1].pcolormesh(
        MD_EXP["c2d"]["rho"],
        MD_EXP["c2d"]["z"],
        np.log10(pcold2dexp),
        vmin=-11,
        vmax=-7,
        cmap="jet",
    )
    ax[1].contour(
        MD_EXP["c2d"]["rho"],
        MD_EXP["c2d"]["z"],
        np.log10(MD_EXP["v2d"]["alpha"]),
        levels=[-3, -2] + list(np.arange(-1.5, 0, 0.1)),
        linestyles="solid",
        colors="m",
    )  # specify levels to plot contour lines at specific z function values
    ax[1].contour(
        MD_EXP["c2d"]["rho"],
        MD_EXP["c2d"]["z"],
        betaexp,
        levels=[beta_level],
        linestyles="solid",
        colors="w",
    )
    ax[1].set_title(rf"EXPANDED: $\beta$={beta_level}")
    ax[1].set_xlabel(r"$\rho_{CYL}/R_J$")

    # Format
    [axi.set_aspect("equal", "box") for axi in ax]
    [axi.set_xlim([0.01, max(MD_EXP["v1d"]["rEq"])]) for axi in ax]
    [axi.set_ylim([-14.99, 14.99]) for axi in ax]
    plt.colorbar(pc, ax=ax[1], label="$log_{10} P_{COLD}$", shrink=0.23)


def plot_2d_cmp_exp_jphi(MD_CMP, MD_EXP):
    """Plot azimuthal current density color map overlayed with B-field configuration for compressed and expanded MD.

    Args:
        MD_CMP (dict): Compressed magnetodisc model.
        MD_EXP (dict): Expanded magnetodisc model.
    """
    # Pressure calc
    pcold2dcmp = calc_Pcold(MD_CMP)
    pcold2dexp = calc_Pcold(MD_EXP)

    betacmp = calc_beta(MD_CMP)
    betaexp = calc_beta(MD_EXP)

    # current calc
    jphi2dcmp = calc_Jphi(MD_CMP)

    jphi2dexp = calc_Jphi(MD_EXP)

    ma_per_rjsqr_si = 1e6 / (MD_EXP["scales"]["length"] ** 2)

    # Plot
    fig, ax = plt.subplots(1, 2, sharey=True, figsize=(8, 5), constrained_layout=True)

    # Compressed MD
    beta_level = 0.2
    pc = ax[0].pcolormesh(
        MD_CMP["c2d"]["rho"],
        MD_CMP["c2d"]["z"],
        jphi2dcmp / ma_per_rjsqr_si,
        vmin=-0.43,
        vmax=5,
        cmap="jet",
    )
    ax[0].contour(
        MD_CMP["c2d"]["rho"],
        MD_CMP["c2d"]["z"],
        np.log10(MD_CMP["v2d"]["alpha"]),
        levels=[-3, -2] + list(np.arange(-1.5, 0, 0.1)),
        linestyles="solid",
        colors="m",
    )  # specify levels to plot contour lines at specific z function values
    ax[0].contour(
        MD_CMP["c2d"]["rho"],
        MD_CMP["c2d"]["z"],
        betacmp,
        levels=[0.2],
        linestyles="solid",
        colors="w",
    )
    ax[0].set_title(rf"COMPRESSED: $\beta$={beta_level}")
    ax[0].set_ylabel(r"$Z/R_J$")
    ax[0].set_xlabel(r"$\rho_{CYL}/R_J$")

    # Expanded MD
    pc = ax[1].pcolormesh(
        MD_EXP["c2d"]["rho"],
        MD_EXP["c2d"]["z"],
        jphi2dexp / ma_per_rjsqr_si,
        vmin=-0.43,
        vmax=5,
        cmap="jet",
    )
    ax[1].contour(
        MD_EXP["c2d"]["rho"],
        MD_EXP["c2d"]["z"],
        np.log10(MD_EXP["v2d"]["alpha"]),
        levels=[-3, -2] + list(np.arange(-1.5, 0, 0.1)),
        linestyles="solid",
        colors="m",
    )  # specify levels to plot contour lines at specific z function values
    ax[1].contour(
        MD_EXP["c2d"]["rho"],
        MD_EXP["c2d"]["z"],
        betaexp,
        levels=[0.2],
        linestyles="solid",
        colors="w",
    )
    ax[1].set_title(rf"EXPANDED: $\beta$={beta_level}")
    ax[1].set_xlabel(r"$\rho_{CYL}/R_J$")

    # Format
    [axi.set_aspect("equal", "box") for axi in ax]
    [axi.set_xlim([0.01, max(MD_EXP["v1d"]["rEq"])]) for axi in ax]
    [axi.set_ylim([-14.99, 14.99]) for axi in ax]
    plt.colorbar(pc, ax=ax[1], label=r"$J_{\Phi} (MA R_J^{  -2})$", shrink=0.23)


def plot_model_cut(rho1, rho2, z1, z2, MD_INP):
    """'LINEAR CUT' PROFILES

    Args:
        rho1 ([type]): [description]
        rho2 ([type]): [description]
        z1 ([type]): [description]
        z2 ([type]): [description]
        MD_INP ([type]): [description]
    """
    return


def plot_command():
    pass
