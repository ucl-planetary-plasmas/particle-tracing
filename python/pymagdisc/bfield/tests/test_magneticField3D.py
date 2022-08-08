from ..magneticField3D import dipoleMagneticField3D, _calc_posn
import numpy as np
import pytest

params = {
    "z_aligned_dipole_no_offset": (
        [0, 0, 1],
        [0, 0, 0],
        [
            [
                4.00087565e-02,
                5.13388429e-02,
                6.44476969e-02,
                7.93932920e-02,
                9.62145308e-02,
                1.14930717e-01,
                1.35541209e-01,
                1.58025253e-01,
                1.82342006e-01,
                2.08430742e-01,
            ],
            [
                -9.79931913e-18,
                -1.25743899e-17,
                -1.57851331e-17,
                -1.94457482e-17,
                -2.35657634e-17,
                -2.81499070e-17,
                -3.31980215e-17,
                -3.87050240e-17,
                -4.46609108e-17,
                -5.10508082e-17,
            ],
            [
                -1.09923155e-01,
                -1.28238230e-01,
                -1.47260124e-01,
                -1.66791901e-01,
                -1.86630153e-01,
                -2.06566865e-01,
                -2.26391332e-01,
                -2.45892093e-01,
                -2.64858879e-01,
                -2.83084555e-01,
            ],
        ],
    ),
}


@pytest.mark.parametrize("m, rm, r", list(params.values()), ids=list(params.keys()))
def test_dipoleMagneticField3D_B_calc(m, rm, r):
    """Test calculation of B in Cartesian from polar coordinates

    Args:
        m (list, array): Magnetic moment vector. Shape = (3,)
        rm (list, array): Position of the magnetic moment. shape = (3,)
        r (list): Cartesion coordinates where to calculate magnetic field. Shape = (3, -1)
    """
    # Get position variables
    m, rm, r, X, Y, Z, r0, rnorm, R, R2, R3 = _calc_posn(m, rm, r)

    # B from dipoleMagneticField3D calculation
    B, _ = dipoleMagneticField3D(m, rm, r)

    # radial distance from z-axis in cylindrical coord
    Rcyl = np.sqrt(X ** 2 + Y ** 2)

    # colatitude
    theta = np.pi / 2 - np.arctan2(Z, Rcyl)

    # cos and sin of colatitude
    cost = Z / R
    sint = np.sqrt(1 - cost ** 2)

    # B radial and tangential
    Br = m[2] / R3 * 2 * cost
    Bt = m[2] / R3 * sint

    # cos and sin of latitude
    coslat = Rcyl / R
    sinlat = Z / R

    # cos and sin of longitude
    coslon = X / Rcyl
    sinlon = Y / Rcyl

    # B in polar to cartesian
    Bx = (Br * coslat + Bt * sinlat) * coslon
    By = (Br * coslat + Bt * sinlat) * sinlon
    Bz = Br * sinlat - Bt * coslat

    # Check the dipoleMagneticField3D is the same as B in Carteisan from polar coordinates.
    np.testing.assert_allclose(
        B[0],
        Bx,
        atol=1e-3,
        err_msg="Error in calculation of B in Cartesian from polar coordinates.",
    )
    np.testing.assert_allclose(
        B[2],
        Bz,
        atol=1e-3,
        err_msg="Error in calculation of B in Cartesian from polar coordinates.",
    )


params = {
    "z_aligned_dipole_no_offset": (
        [0, 0, 1],
        [0, 0, 0],
        [
            [
                4.00087565e-02,
                5.13388429e-02,
                6.44476969e-02,
                7.93932920e-02,
                9.62145308e-02,
                1.14930717e-01,
                1.35541209e-01,
                1.58025253e-01,
                1.82342006e-01,
                2.08430742e-01,
            ],
            [
                -9.79931913e-18,
                -1.25743899e-17,
                -1.57851331e-17,
                -1.94457482e-17,
                -2.35657634e-17,
                -2.81499070e-17,
                -3.31980215e-17,
                -3.87050240e-17,
                -4.46609108e-17,
                -5.10508082e-17,
            ],
            [
                -1.09923155e-01,
                -1.28238230e-01,
                -1.47260124e-01,
                -1.66791901e-01,
                -1.86630153e-01,
                -2.06566865e-01,
                -2.26391332e-01,
                -2.45892093e-01,
                -2.64858879e-01,
                -2.83084555e-01,
            ],
        ],
    ),
}


@pytest.mark.parametrize("m, rm, r", list(params.values()), ids=list(params.keys()))
def test_dipoleMagneticField3D_gradB_calc(m, rm, r):
    """Test calculation of gradB in Cartesian from polar coordinates

    Args:
        m (list, array): Magnetic moment vector. Shape = (3,)
        rm (list, array): Position of the magnetic moment. shape = (3,)
        r (list): Cartesion coordinates where to calculate magnetic field. Shape = (3, -1)
    """
    # Get position variables
    m, rm, r, X, Y, Z, r0, rnorm, R, R2, R3 = _calc_posn(m, rm, r)

    # B from dipoleMagneticField3D calculation
    B, gradB = dipoleMagneticField3D(m, rm, r)

    # radial distance from z-axis in cylindrical coord
    Rcyl = np.sqrt(X ** 2 + Y ** 2)

    # colatitude
    theta = np.pi / 2 - np.arctan2(Z, Rcyl)

    # cos and sin of colatitude
    cost = Z / R

    # cos and sin of latitude
    coslat = Rcyl / R
    sinlat = Z / R

    # cos and sin of longitude
    coslon = X / Rcyl
    sinlon = Y / Rcyl

    EPS = 1e-6
    rp = R * (1 + EPS)
    rm = R * (1 - EPS)
    dr = 2 * R * EPS

    # gradient with respect to r
    dBdr = m[2] * (1 / rp ** 3 - 1 / rm ** 3) / dr * np.sqrt(1 + 3 * cost ** 2)

    mup = cost * (1 + EPS) + EPS
    mum = cost * (1 - EPS) - EPS
    rdmu = R * (mup - mum)
    dmudt = -np.sqrt(1 - cost ** 2)

    # gradient with respect to theta
    dBdt = (
        m[2]
        / R3
        * (np.sqrt(1 + 3 * mup ** 2) - np.sqrt(1 + 3 * mum ** 2))
        / rdmu
        * dmudt
    )
    gradBx = (dBdr * coslat + dBdt * sinlat) * coslon
    gradBy = (dBdr * coslat + dBdt * sinlat) * sinlon
    gradBz = dBdr * sinlat - dBdt * coslat

    # Check the dipoleMagneticField3D is the same as B in Carteisan from polar coordinates.
    np.testing.assert_allclose(
        gradB[0],
        gradBx,
        atol=3e3,
        err_msg="Error in calculation of B in Cartesian from polar coordinates.",
    )
    np.testing.assert_allclose(
        gradB[2],
        gradBz,
        atol=3e3,
        err_msg="Error in calculation of B in Cartesian from polar coordinates.",
    )
