from scipy.interpolate import RectBivariateSpline as interp2


def magneticFieldInterpolator(MD):
    """Create interpolate object for evaluating magnetodisc magnetic field at given positions (r, mu).

    Args:
        MD (dict): Magnetodisc model data.

    Returns:
        interpolator (dict): Interpolators for Br, Bth, and B.
    """
    MD_Br = MD["v2d"]["Br"] * MD["scales"]["bfield"]
    MD_Bth = MD["v2d"]["Bth"] * MD["scales"]["bfield"]
    MD_B = MD["v2d"]["B"] * MD["scales"]["bfield"]

    interpolator = {key: None for key in ["Br", "Bth", "B"]}
    interpolator["Br"] = interp2(MD["c2d"]["r"][0, :], MD["c2d"]["mu"][:, 0], MD_Br.T)
    interpolator["Bth"] = interp2(MD["c2d"]["r"][0, :], MD["c2d"]["mu"][:, 0], MD_Bth.T)
    interpolator["B"] = interp2(MD["c2d"]["r"][0, :], MD["c2d"]["mu"][:, 0], MD_B.T)
    return interpolator
