import numpy as np

from . import coords, maths


def _rotmat_x(angle):
    """Rotation matrix for rotation around the x-axis.

    Parameters
    ----------
    angle : float
        Angle of rotation.

    Returns
    -------
    rotx : array
        Rotation matrix.
    """
    rotx = np.array([1., 0., 0., 0., np.cos(angle), -np.sin(angle), 0., np.sin(angle), np.cos(angle)])
    return rotx


def _rotmat_y(angle):
    """Rotation matrix for rotation around the y-axis.

    Parameters
    ----------
    angle : float
        Angle of rotation.

    Returns
    -------
    roty : array
        Rotation matrix.
    """
    roty = np.array([np.cos(angle), 0., np.sin(angle), 0., 1., 0., -np.sin(angle), 0., np.cos(angle)])
    return roty


def _rotmat_z(angle):
    """Rotation matrix for rotation around the z-axis.

    Parameters
    ----------
    angle : float
        Angle of rotation.

    Returns
    -------
    rotz : array
        Rotation matrix.
    """
    rotz = np.array([np.cos(angle), -np.sin(angle), 0., np.sin(angle), np.cos(angle), 0., 0., 0., 1.])
    return rotz


def _rotmat_axis(angle, axis):
    """Rotation matrix for any given axis.

    Parameters
    ----------
    angle : float
        Angle of rotation.
    axis : int
        Axis of rotation, 0 = x-axis, 1 = y-axis and 2 = z-axis.
    """
    if axis == 0:
        return _rotmat_x(angle)
    elif axis == 1:
        return _rotmat_y(angle)
    elif axis == 2:
        return _rotmat_z(angle)


def _rotmat_euler(angles, axes):
    """Euler rotational matrix.

    Parameters
    ----------
    angles : float
        Angles of rotation.
    axes : int
        Axes of rotation.

    Returns
    -------
    _rot : array
        Euler rotational matrix.
    """
    _rot2 = _rotmat_axis(angles[2], axes[2])
    _rot1 = _rotmat_axis(angles[1], axes[1])
    _rot0 = _rotmat_axis(angles[0], axes[0])
    _rot21 = maths.matrix_dot_3by3(_rot2, _rot1)
    _rot = maths.matrix_dot_3by3(_rot21, _rot0)
    return _rot


def _rotate_3d(x, y, z, rot):
    """Rotates cartesian coordinates given an input rotational matrix.

    Parameters
    ----------
    x, y, z : array
        3D cartesian coordinates to rotate.
    rot : array
        Rotational matrix.

    Returns
    -------
    xrot, yrot, zrot : array
        Rotated 3D cartesian coordinates.
    """
    xrot = rot[0]*x + rot[1]*y + rot[2]*z
    yrot = rot[3]*x + rot[4]*y + rot[5]*z
    zrot = rot[6]*x + rot[7]*y + rot[8]*z
    return xrot, yrot, zrot


def rotate3d_Euler(x, y, z, angles, axes='zyz', center=[0., 0., 0.]):
    """Rotates points in 3D cartesian coordinates by Euler angle around
    specified axes.

    Parameters
    ----------
    x, y, z : float or array
        Cartesian coordinates.
    angles : array
        Euler angles.
    axes : str, optional
        Euler angle axes, default z-y-z rotations.
    center : list[float], optional
        Center of rotation, default=[0., 0., 0.].
    k : array, optional
        If k is specified, points are rotated around a unit vector k.

    Returns
    -------
    xrot, yrot, zrot : float or array
        Rotated x, y and z coordinates.
    """
    assert len(angles) == len(axes), "Length of angles and axes must be consistent."
    assert len(angles) == 3, "Length of angles and axes must be == 3."
    axes_int = []
    for i in range(0, len(axes)):
        if axes[i] == 'x':
            axes_int.append(0)
        elif axes[i] == 'y':
            axes_int.append(1)
        elif axes[i] == 'z':
            axes_int.append(2)
    _x, _y, _z = x-center[0], y-center[1], z-center[2]
    rot = _rotmat_euler(angles, axes_int)
    _x, _y, _z = _rotate_3d(_x, _y, _z, rot)
    xrot = _x + center[0]
    yrot = _y + center[1]
    zrot = _z + center[2]
    return xrot, yrot, zrot


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
    x, y, z = coords.sphere2cart(r, phi, the)
    x, y, z = rotate3d_Euler(x, y, z, angles, axes='zyz', center=[0., 0., 0.])
    _, phi, the = coords.cart2sphere(x, y, z)
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
    x1, y1, z1 = coords.sphere2cart(1., phi1, the1)
    x2, y2, z2 = coords.sphere2cart(1., phi2, the2)
    xm = 0.5*(x1+x2)
    ym = 0.5*(y1+y2)
    zm = 0.5*(z1+z2)
    _, midphi, midthe = coords.cart2sphere(xm, ym, zm)
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
