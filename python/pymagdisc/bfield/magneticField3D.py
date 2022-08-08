import numpy as np
from scipy.interpolate import RectBivariateSpline as interp2
from pymagdisc.data.load_data import load_model
from pymagdisc import config


def _calc_posn(m, rm, r):
    """Calculate position related variables.

    Args:
        m (list, array): Magnetic moment vector. Shape = (3,)
        rm (list, array): Position of the magnetic moment. shape = (3,)
        r (list): Cartesion coordinates where to calculate magnetic field. Shape = (3, -1)

    Returns:
        (tuple): Position related variables used by dipoleMagneticField3D().
    """

    # Convert to array
    m = np.array(m).reshape(3, 1)
    rm = np.array(rm).reshape(3, 1)
    r = np.array(r).reshape(3, -1)

    # Subtract moment position
    X, Y, Z = r - rm
    r0 = np.array([X, Y, Z])

    # modulus square
    R2 = X ** 2 + Y ** 2 + Z ** 2

    # modulus
    R = R2 ** 0.5

    # modulus cube
    R3 = R2 ** 1.5

    # normalised direciton
    rnorm = np.array([X / R, Y / R, Z / R])

    return m, rm, r, X, Y, Z, r0, rnorm, R, R2, R3


def dipoleMagneticField3D(m, rm, r) -> tuple:
    """Calculate the 3D dipole magnetic field at position r.

    Args:
        m (list, array): Magnetic moment vector. Shape = (3,)
        rm (list, array): Position of the magnetic moment. shape = (3,)
        r (list): Cartesion coordinates where to calculate magnetic field. Shape = (3, -1)

    Returns:
        B (array): Magnetic field at position r.
        gradB (array): Gradient of dipole magnetic field at position r. Shape (len(m), len(r)).

    """

    # Get position variables
    m, rm, r, X, Y, Z, r0, rnorm, R, R2, R3 = _calc_posn(m, rm, r)

    # scalar product M.r; shape=(-1,)
    mDotR = (m.T @ rnorm).sum(axis=0)

    # Calculate Bx, By and Bz
    B = np.zeros(shape=(4, len(mDotR)))
    B[0:3, :] = (3 * mDotR * rnorm - m) / R3

    # |B|^2
    B[3, :] = (3 * mDotR ** 2 + np.sum(m ** 2)) / R3 ** 2

    # d|B|/dx, d|B|/dy, d|B|/dz
    fR6 = 3 / R3 ** 2
    u = np.sqrt(3 * mDotR ** 2 + np.sum(m ** 2))
    fac1 = R2 * mDotR / u
    fac2 = R * u

    gradB = np.zeros(shape=(3, len(mDotR)))
    gradB[0:3, :] = fR6 * (fac1 * m * (1 - rnorm ** 2) - r0 * fac2)

    return B, gradB


def mdiscMagneticField3D(r, MD, interpolator, compute_gradB=False):
    """Calculate the 3D magnetodisc magnetic field at position r.

    \vec{B}(\vec{m},\vec{r} ) = \frac{\mu_0}{4\pi r^3}
    \left(3(\vec{m}\cdot\vec{\hat{r}}) \vec{\hat{r}} -\vec{m} \right) +
    \frac{2\mu_0}{3} \vec{m} \delta^3(\vec{r})

    Args:
        r (list): Cartesion coordinates where to calculate magnetic field. Shape = (3, -1)
        MD (dict): Magnetodisc model data.
        interpolator (dict): Interpolate object for evaluating B-field at given points (r, mu).
        compute_gradB (bool, optional): Default is False.

    Returns:
        B (array): Magnetic field at position r.
        gradB (array): Gradient of magnetodisc magnetic field at position r. Shape (3, len(r)).
    """
    # Separate coordinates R=(X,Y,Z)
    X, Y, Z = r

    # R modulus
    R = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)

    # Projection on (x,y) plane
    Rcyl = np.sqrt(X ** 2 + Y ** 2)

    # colatitude in radians
    theta = np.pi / 2 - np.arctan2(Z, Rcyl)

    # cos and sin of latitude
    coslat = Rcyl / R
    sinlat = Z / R
    assert (
        coslat ** 2 + sinlat ** 2 - 1 < 1e-6
    ).all(), "Incorrect cos and sin of latitude. Sum of squares not equal to 1."

    # cos and sin of longitude
    coslon = X / Rcyl
    sinlon = Y / Rcyl
    assert (
        coslon ** 2 + sinlon ** 2 - 1 < 1e-6
    ).all(), "Incorrect cos and sin of longitude. Sum of squares not equal to 1."

    # Calculate Br, Btheta, gradB based on the MD file
    Br, Bt, gradB = MDiscField(
        R, theta, MD, compute_gradB=compute_gradB, interpolator=interpolator
    )

    # |B|^2
    B2 = Br ** 2 + Bt ** 2

    # Calculate corresponding Bx, By and Bz
    Bx = (Br * coslat + Bt * sinlat) * coslon
    By = (Br * coslat + Bt * sinlat) * sinlon
    Bz = Br * sinlat - Bt * coslat
    B = np.array([Bx, By, Bz, B2])

    if compute_gradB:
        # d|B|/dr, d|B|/dtheta
        dBr = gradB[0]
        dBt = gradB[1]

        # d|B|/dx, d|B|/dy, d|B|/dz
        gradBx = (dBr * coslat + dBt * sinlat) * coslon
        gradBy = (dBr * coslat + dBt * sinlat) * sinlon
        gradBz = dBr * sinlat - dBt * coslat
        gradB = np.array([gradBx, gradBy, gradBz])
    else:
        gradB = None

    return B, gradB


def MDiscField(r, theta, MD, interpolator, compute_gradB=False):
    """Calculate Br and Btheta at position r and theta based on the MD file.

    Args:
        r (array): radial distance in [RJ]
        theta (array): colatitude in [radians] (0 degrees from z-axis)
        MD (dict): Magnetodisc model data.
        interpolator (dict): Interpolate object for evaluating B-field at given points (r, mu).
        compute_gradB (bool, optional): Default is False.

    Returns:
        Br (array): radial field in [T]
        Bt (array): meridional field in [T]
        gradB (array): gradient(r, theta) of magnetic field magnitude in [T/m]

    """

    # Define mu = cos(theta)
    mu = np.cos(theta)

    # Field components in model (tilted, rotating) frame
    # fast gridded data interpolant

    Br = interpolator["Br"].ev(r, mu)
    Bth = interpolator["Bth"].ev(r, mu)

    EPS = 1e-6

    # radial component
    if compute_gradB:
        dr = 2 * r * EPS * MD["scales"]["length"]
        rp = r * (1 + EPS)
        rm = r * (1 - EPS)
        dBdr = (interpolator["B"].ev(rp, mu) - interpolator["B"].ev(rm, mu)) / dr

        # tangential component
        dmudt = -np.sqrt(1 - mu ** 2)

        # add EPS to deal with mu=0
        mup = mu * (1 + EPS) + EPS
        mum = mu * (1 - EPS) - EPS
        rdmu = r * (mup - mum) * MD["scales"]["length"]
        dBdt = (
            (interpolator["B"].ev(r, mup) - interpolator["B"].ev(r, mum)) / rdmu * dmudt
        )
        gradB = np.array([dBdr, dBdt])
    else:
        gradB = None

    return Br, Bth, gradB


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
