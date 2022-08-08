import numpy as np


def calc_Pcold(MD_INP):
    """Calculate cold plasma pressure.

    Args:
        MD_INP (dict): Magnetodisc model input.
    Returns:
        pcold2d (array): 2d Cold plasma pressure.
    """
    pcold2d = MD_INP["v2d"]["PcEq"] * np.exp(
        -1.0
        * (MD_INP["v2d"]["rho0"] ** 2 - MD_INP["c2d"]["rho"] ** 2)
        / (2.0 * MD_INP["v2d"]["confLen"] ** 2)
    )

    return pcold2d


def calc_Jphi(MD_INP):
    """Calculate azimuthal current density.

    Args:
        MD_INP (dict): Magnetodisc model input.

    Returns:
        jphi2d (array): 2d azimuthal current density.
    """
    jphi2d = (
        MD_INP["v2d"]["jphi"]
        * MD_INP["scales"]["pressure"]
        / (MD_INP["scales"]["length"] * MD_INP["scales"]["bfield"])
    )
    return jphi2d


def calc_beta(MD_INP):
    """Calculate cold plasma beta.

    Args:
        MD_INP (dict): Magnetodisc model input.

    Returns:
        beta (array): 2d cold plasma beta.
    """
    pcold2d = calc_Pcold(MD_INP)
    beta = pcold2d / (0.5 * MD_INP["v2d"]["B"] ** 2)
    return beta
