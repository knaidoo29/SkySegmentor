import numpy as np
import magpie


def rotate_usphere(phi, the, angles):
    """Rotates spherical coordinates by Euler angles performed along the z-axis,
    then y-axis and then z-axis.

    Parameters
    ----------
    phi, the : float or array
        Spherical angular coordinates.
    angles : list
        Euler angles defining rotations about the z-axis, y-axis then z-axis.
    """
    if np.isscalar(phi):
        r = 1.
    else:
        r = np.ones(len(phi))
    x, y, z = magpie.coords.sphere2cart(r, phi, the)
    x, y, z = magpie.coords.rotate3d_Euler(x, y, z, angles, axes='zyz', center=[0., 0., 0.])
    _, phi, the = magpie.coords.cart2sphere(x, y, z)
    return phi, the


def midpoint_usphere(phi1, phi2, the1, the2):
    """Finds the spherical angular coordinates of the midpoint between two points
    on a unit sphere.

    Parameters
    ----------
    phi1, phi2 : float
        Longitude coordinates for both points.
    the1, the2 : float
        Latitude coordinates for both points.

    Returns
    -------
    midphi, midthe : float
        Midpoint along the longitude phi and latitude theta.
    """
    x1, y1, z1 = magpie.coords.sphere2cart(1., phi1, the1)
    x2, y2, z2 = magpie.coords.sphere2cart(1., phi2, the2)
    xm = 0.5*(x1+x2)
    ym = 0.5*(y1+y2)
    zm = 0.5*(z1+z2)
    _, midphi, midthe = magpie.coords.cart2sphere(xm, ym, zm)
    return midphi, midthe


def rotate2plane(c1, c2):
    """Finds the rotation angles to place the two coordinates c1 and c2 along a
    latitude = pi/2 (i.e. equator of sphere) and with a midpoint of longitude = pi.

    Parameters
    ----------
    c1, c2 : float
        Coordinates of two points where c1 = [phi1, theta1] and c2 = [phi2, theta2].

    Returns
    -------
    a1, a2, a3 : lists
        Euler angles of rotation.
    """
    _c1 = np.copy(c1)
    _c2 = np.copy(c2)
    cen1_phi, cen1_the = _c1[0], _c1[1]
    cen2_phi, cen2_the = _c2[0], _c2[1]
    cenm_phi, cenm_the = midpoint_usphere(cen1_phi, cen2_phi, cen1_the, cen2_the)
    a1 = np.copy([-cenm_phi, -cenm_the, 0.])
    cen1_phi, cen1_the = rotate_usphere(cen1_phi, cen1_the, a1)
    cen2_phi, cen2_the = rotate_usphere(cen2_phi, cen2_the, a1)
    a2 = np.copy([np.pi-cen1_phi, 0., 0.])
    cen1_phi, cen1_the = rotate_usphere(cen1_phi, cen1_the, a2)
    cen2_phi, cen2_the = rotate_usphere(cen2_phi, cen2_the, a2)
    a3 = np.copy([np.pi/2., np.pi/2., np.pi])
    cen1_phi, cen1_the = rotate_usphere(cen1_phi, cen1_the, a3)
    cen2_phi, cen2_the = rotate_usphere(cen2_phi, cen2_the, a3)
    c1 = [cen1_phi, cen1_the]
    c2 = [cen2_phi, cen2_the]
    return a1, a2, a3


def forward_rotate(phi, the, a1, a2, a3):
    """Applies a forward rotation of spherical angular coordinates phi and theta
    using the forward Euler angles of rotation a1, a2 and a3.

    Parameters
    ----------
    phi, the : float or array
        Spherical angular coordinates.
    a1, a2, a3 : lists
        Euler angles of rotation.
    """
    _phi, _the = np.copy(phi), np.copy(the)
    _phi, _the = rotate_usphere(_phi, _the, a1)
    _phi, _the = rotate_usphere(_phi, _the, a2)
    _phi, _the = rotate_usphere(_phi, _the, a3)
    return _phi, _the


def backward_rotate(phi, the, a1, a2, a3):
    """Applies a backward rotation of spherical angular coordinates phi and theta
    using the forward Euler angles of rotation a1, a2 and a3.

    Parameters
    ----------
    phi, the : float or array
        Spherical angular coordinates.
    a1, a2, a3 : lists
        Euler angles of rotation.
    """
    _phi, _the = np.copy(phi), np.copy(the)
    ra3 = -np.copy(a3[::-1])
    ra2 = -np.copy(a2[::-1])
    ra1 = -np.copy(a1[::-1])
    _phi, _the = rotate_usphere(_phi, _the, ra3)
    _phi, _the = rotate_usphere(_phi, _the, ra2)
    _phi, _the = rotate_usphere(_phi, _the, ra1)
    return _phi, _the
