import numpy as np


def vector_norm(a):
    """Returns the magnitude a vector.

    Parameters
    ----------
    a : array
        Vector a.
    """
    a1, a2, a3 = a[0], a[1], a[2]
    return np.sqrt(a1**2. + a2**2. + a3**2.)


def vector_dot(a, b):
    """Returns the vector dot product.

    Parameters
    ----------
    a : array
        Vector a.
    b : array
        Vector b.
    """
    a1, a2, a3 = a[0], a[1], a[2]
    b1, b2, b3 = b[0], b[1], b[2]
    return a1*b1 + a2*b2 + a3*b3


def vector_cross(a, b):
    """Returns the vector cross product.

    Parameters
    ----------
    a : array
        Vector a.
    b : array
        Vector b.

    Returns
    -------
    s : array
        Cross product vector.
    """
    a1, a2, a3 = a[0], a[1], a[2]
    b1, b2, b3 = b[0], b[1], b[2]
    s1 = a2*b3 - a3*b2
    s2 = a3*b1 - a1*b3
    s3 = a1*b2 - a2*b1
    s = np.array([s1, s2, s3])
    return s


def matrix_dot_3by3(mat1, mat2):
    """Dot product of 2 3x3 matrices.

    Parameters
    ----------
    mat1, mat2 : array
        Flattened 3by3 matrices.

    Returns
    -------
    mat3 : array
        Matrix dot product output.
    """
    mat3 = np.array([
        mat1[0]*mat2[0] + mat1[1]*mat2[3] + mat1[2]*mat2[6],
        mat1[0]*mat2[1] + mat1[1]*mat2[4] + mat1[2]*mat2[7],
        mat1[0]*mat2[2] + mat1[1]*mat2[5] + mat1[2]*mat2[8],
        mat1[3]*mat2[0] + mat1[4]*mat2[3] + mat1[5]*mat2[6],
        mat1[3]*mat2[1] + mat1[4]*mat2[4] + mat1[5]*mat2[7],
        mat1[3]*mat2[2] + mat1[4]*mat2[5] + mat1[5]*mat2[8],
        mat1[6]*mat2[0] + mat1[7]*mat2[3] + mat1[8]*mat2[6],
        mat1[6]*mat2[1] + mat1[7]*mat2[4] + mat1[8]*mat2[7],
        mat1[6]*mat2[2] + mat1[7]*mat2[5] + mat1[8]*mat2[8]
    ])
    return mat3
