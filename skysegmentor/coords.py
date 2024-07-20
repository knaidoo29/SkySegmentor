import numpy as np

from . import maths, utils


def cart2sphere(x, y, z, center=[0., 0., 0.]):
    """Return polar coordinates for a given set of cartesian coordinates.

    Parameters
    ----------
    x : array
        x coordinate.
    y : array
        y coordinate.
    center : list
        Center point of polar coordinate grid.

    Returns
    -------
    r : array
        Radial coordinates.
    phi : array
        Phi coordinates [0, 2pi].
    theta : array
        Theta coordinates [0, pi].
    """
    r = np.sqrt((x-center[0])**2. + (y-center[1])**2. + (z-center[2])**2.)
    phi = np.arctan2(y-center[1], x-center[0])
    if utils.isscalar(phi) is True:
        if phi < 0.:
            phi += 2.*np.pi
        if r != 0.:
            theta = np.arccos((z-center[2])/r)
        else:
            theta = 0.
    else:
        condition = np.where(phi < 0.)
        phi[condition] += 2.*np.pi
        theta = np.zeros(len(phi))
        condition = np.where(r != 0.)[0]
        theta[condition] = np.arccos((z[condition]-center[2])/r[condition])
    return r, phi, theta


def sphere2cart(r, phi, theta, center=[0., 0., 0.]):
    """Converts spherical polar coordinates into cartesian coordinates.

    Parameters
    ----------
    r : array
        Radial distance.
    phi : array
        Longitudinal coordinates (radians = [0, 2pi]).
    theta : array
        Latitude coordinates (radians = [0, pi]).
    center : list
        Center point of spherical polar coordinate grid.

    Returns
    -------
    x, y, z : array
        Euclidean coordinates.
    """
    phi = np.copy(phi)
    theta = np.copy(theta)
    x = r * np.cos(phi) * np.sin(theta)
    y = r * np.sin(phi) * np.sin(theta)
    z = r * np.cos(theta)
    x += center[0]
    y += center[1]
    z += center[2]
    return x, y, z


def distusphere(phi1, theta1, phi2, theta2):
    """Compute angular (great-arc) distance between two points on a unit
    sphere.

    Parameters
    ----------
    phi1, theta1 : float or array
        Location of first points on the unit sphere.
    phi2, theta2 : float or array
        Location of second points on the unit sphere.

    Returns
    -------
    dist : float or array
        Angular great-arc distance.
    """
    if utils.isscalar(phi1) is True:
        r = 1.
    else:
        r = np.ones(len(phi1))
    x1, y1, z1 = sphere2cart(r, phi1, theta1)
    x2, y2, z2 = sphere2cart(r, phi2, theta2)
    a = np.array([x1, y1, z1])
    b = np.array([x2, y2, z2])
    cross = maths.vector_cross(a, b)
    normcross = maths.vector_norm(cross)
    dot = maths.vector_dot(a, b)
    dist = np.arctan2(normcross, dot)
    return dist
