from scipy.interpolate import RectBivariateSpline as interp2


def magneticFieldInterpolator(MD):
    """Create interpolate object for evaluating magnetodisc magnetic field at given positions (r, mu).

    Args:
        MD (dict): Magnetodisc model data.

    Returns:
        interpolator (dict): Interpolators for Br, Bth, and B.
    """
    opts = ("cubic", "cubic")

    print(f"Generating field functions ({opts[0]},{opts[1]}) from ", end="")
    print("magnetodisk field ... ", end="")

    R = MD["c2d"]["r"][0, :]
    Mu = MD["c2d"]["mu"][:, 0]

    interpolator = {key: None for key in ["Br", "Bth", "B"]}

    Br = MD["v2d"]["Br"] * MD["scales"]["bfield"]
    interpolator["Br"] = interp2(R, Mu, Br.T, kx=3, ky=3)

    Bth = MD["v2d"]["Bth"] * MD["scales"]["bfield"]
    interpolator["Bth"] = interp2(R, Mu, Bth.T)

    B = MD["v2d"]["B"] * MD["scales"]["bfield"]
    interpolator["B"] = interp2(R, Mu, B.T)

    print("done.")

    return interpolator
