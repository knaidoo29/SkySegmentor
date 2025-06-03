import numpy as np
import skysegmentor


def test_zero_vector():
    v = np.array([0.0, 0.0, 0.0])
    result = skysegmentor.vector_norm(v)
    assert np.isclose(result, 0.0)


def test_unit_vectors():
    v = np.array([1.0, 0.0, 0.0])
    assert np.isclose(skysegmentor.vector_norm(v), 1.0)
    v = np.array([0.0, 1.0, 0.0])
    assert np.isclose(skysegmentor.vector_norm(v), 1.0)
    v = np.array([0.0, 0.0, 1.0])
    assert np.isclose(skysegmentor.vector_norm(v), 1.0)


def test_positive_components():
    v = np.array([3.0, 4.0, 0.0])
    result = skysegmentor.vector_norm(v)
    expected = 5.0  # 3-4-5 triangle
    assert np.isclose(result, expected)


def test_negative_components():
    v = np.array([-3.0, -4.0, 0.0])
    result = skysegmentor.vector_norm(v)
    expected = 5.0
    assert np.isclose(result, expected)


def test_mixed_components():
    v = np.array([1.0, -2.0, 2.0])
    result = skysegmentor.vector_norm(v)
    expected = np.sqrt(1.0**2 + (-2.0)**2 + 2.0**2)
    assert np.isclose(result, expected)


def test_non_integer_components():
    v = np.array([1.5, 2.5, -3.5])
    result = skysegmentor.vector_norm(v)
    expected = np.sqrt(1.5**2 + 2.5**2 + (-3.5)**2)
    assert np.isclose(result, expected)


def test_large_values():
    v = np.array([1e10, 1e10, 1e10])
    result = skysegmentor.vector_norm(v)
    expected = np.sqrt(3 * (1e10)**2)
    assert np.isclose(result, expected)


def test_small_values():
    v = np.array([1e-10, -1e-10, 1e-10])
    result = skysegmentor.vector_norm(v)
    expected = np.sqrt(3 * (1e-10)**2)
    assert np.isclose(result, expected)


def test_vector_dot_basic():
    a = np.array([1, 0, 0])
    b = np.array([0, 1, 0])
    assert np.isclose(skysegmentor.vector_dot(a, b), 0.0)

    a = np.array([1, 2, 3])
    b = np.array([4, -5, 6])
    expected = 1*4 + 2*(-5) + 3*6
    assert np.isclose(skysegmentor.vector_dot(a, b), expected)


def test_vector_dot_commutative():
    a = np.array([1, 2, 3])
    b = np.array([4, 5, 6])
    assert np.isclose(skysegmentor.vector_dot(a, b), skysegmentor.vector_dot(b, a))


def test_vector_dot_zero_vector():
    a = np.array([0, 0, 0])
    b = np.array([1, 2, 3])
    assert np.isclose(skysegmentor.vector_dot(a, b), 0.0)


def test_vector_cross_orthogonality():
    a = np.array([1, 0, 0])
    b = np.array([0, 1, 0])
    c = skysegmentor.vector_cross(a, b)
    expected = np.array([0, 0, 1])
    assert np.allclose(c, expected)


def test_vector_cross_antiparallel():
    a = np.array([1, 0, 0])
    b = np.array([-1, 0, 0])
    c = skysegmentor.vector_cross(a, b)
    expected = np.array([0, 0, 0])
    assert np.allclose(c, expected)


def test_vector_cross_parallel():
    a = np.array([2, 2, 2])
    b = a * 3
    c = skysegmentor.vector_cross(a, b)
    expected = np.array([0, 0, 0])
    assert np.allclose(c, expected)


def test_vector_cross_general():
    a = np.array([1, 2, 3])
    b = np.array([4, 5, 6])
    c = skysegmentor.vector_cross(a, b)
    # Manual cross product calculation
    expected = np.array([
        2*6 - 3*5,
        3*4 - 1*6,
        1*5 - 2*4
    ])
    assert np.allclose(c, expected)


def test_matrix_dot_3by3_identity():
    identity = np.array([
        1, 0, 0,
        0, 1, 0,
        0, 0, 1,
    ])
    mat = np.array([
        2, 3, 4,
        5, 6, 7,
        8, 9, 10,
    ])
    result = skysegmentor.matrix_dot_3by3(identity, mat)
    # Should return mat unchanged
    assert np.allclose(result, mat)

    result2 = skysegmentor.matrix_dot_3by3(mat, identity)
    assert np.allclose(result2, mat)


def test_matrix_dot_3by3_known_product():
    mat1 = np.array([
        1, 2, 3,
        0, 1, 4,
        5, 6, 0,
    ])
    mat2 = np.array([
        -1, 0, 0,
        0, 1, 0,
        0, 0, 1,
    ])
    expected = np.dot(mat1.reshape(3, 3), mat2.reshape(3, 3)).flatten()
    result = skysegmentor.matrix_dot_3by3(mat1, mat2)
    assert np.allclose(result, expected)


def test_matrix_dot_3by3_zero_matrix():
    zero = np.zeros(9)
    mat = np.array([
        1, 2, 3,
        4, 5, 6,
        7, 8, 9,
    ])
    result = skysegmentor.matrix_dot_3by3(mat, zero)
    assert np.allclose(result, zero)
    result = skysegmentor.matrix_dot_3by3(zero, mat)
    assert np.allclose(result, zero)