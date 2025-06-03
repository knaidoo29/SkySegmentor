import numpy as np
import skysegmentor

def test_cart2sphere_scalar_origin():
    r, phi, theta = skysegmentor.cart2sphere(0.0, 0.0, 0.0)
    assert r == 0.0
    assert phi == 0.0
    assert theta == 0.0

def test_cart2sphere_scalar_x_axis():
    r, phi, theta = skysegmentor.cart2sphere(1.0, 0.0, 0.0)
    assert np.isclose(r, 1.0)
    assert np.isclose(phi, 0.0)
    assert np.isclose(theta, np.pi / 2)

def test_cart2sphere_scalar_y_axis():
    r, phi, theta = skysegmentor.cart2sphere(0.0, 1.0, 0.0)
    assert np.isclose(phi, np.pi / 2)

def test_cart2sphere_scalar_z_axis():
    r, phi, theta = skysegmentor.cart2sphere(0.0, 0.0, 1.0)
    assert np.isclose(theta, 0.0)

def test_cart2sphere_scalar_negative_phi():
    r, phi, theta = skysegmentor.cart2sphere(-1.0, -1.0, 0.0)
    assert 0.0 <= phi <= 2 * np.pi

def test_cart2sphere_vectorized_input():
    x = np.array([1.0, 0.0])
    y = np.array([0.0, 1.0])
    z = np.array([0.0, 0.0])
    r, phi, theta = skysegmentor.cart2sphere(x, y, z)

    assert np.allclose(r, [1.0, 1.0])
    assert np.allclose(theta, [np.pi / 2, np.pi / 2])
    assert np.all(phi >= 0.0) and np.all(phi <= 2 * np.pi)

def test_cart2sphere_custom_center():
    x = np.array([2.0])
    y = np.array([2.0])
    z = np.array([2.0])
    center = [1.0, 1.0, 1.0]
    r, phi, theta = skysegmentor.cart2sphere(x, y, z, center=center)
    assert np.isclose(r, np.sqrt(3))
    assert np.isclose(phi, np.pi / 4)
    assert np.isclose(theta, np.arccos(1/np.sqrt(3)))

def test_cart2sphere_zero_radius_vector():
    x = np.array([1.0])
    y = np.array([1.0])
    z = np.array([1.0])
    center = [1.0, 1.0, 1.0]
    r, phi, theta = skysegmentor.cart2sphere(x, y, z, center=center)
    assert r == 0.0
    assert phi == 0.0
    assert theta == 0.0

def test_single_point_origin_center():
    # radius 1, phi=0, theta=pi/2 should be (1,0,0)
    r = np.array([1.0])
    phi = np.array([0.0])
    theta = np.array([np.pi / 2])
    x, y, z = skysegmentor.sphere2cart(r, phi, theta)
    assert np.isclose(x, 1.0)
    assert np.isclose(y, 0.0)
    assert np.isclose(z, 0.0)


def test_single_point_different_phi_theta():
    # radius 1, phi=pi/2, theta=pi/2 should be (0,1,0)
    r = np.array([1.0])
    phi = np.array([np.pi / 2])
    theta = np.array([np.pi / 2])
    x, y, z = skysegmentor.sphere2cart(r, phi, theta)
    assert np.isclose(x, 0.0)
    assert np.isclose(y, 1.0)
    assert np.isclose(z, 0.0)


def test_single_point_pole():
    # radius 1, theta=0 (north pole), any phi should be (0,0,1)
    r = np.array([1.0])
    phi = np.array([123.45])  # arbitrary phi
    theta = np.array([0.0])
    x, y, z = skysegmentor.sphere2cart(r, phi, theta)
    assert np.isclose(x, 0.0)
    assert np.isclose(y, 0.0)
    assert np.isclose(z, 1.0)
    

def test_single_point_south_pole():
    # radius 1, theta=pi (south pole), any phi should be (0,0,-1)
    r = np.array([1.0])
    phi = np.array([0.0])
    theta = np.array([np.pi])
    x, y, z = skysegmentor.sphere2cart(r, phi, theta)
    assert np.isclose(x, 0.0)
    assert np.isclose(y, 0.0)
    assert np.isclose(z, -1.0)


def test_multiple_points_array_input():
    r = np.array([1.0, 2.0])
    phi = np.array([0.0, np.pi / 2])
    theta = np.array([np.pi / 2, np.pi / 2])
    x, y, z = skysegmentor.sphere2cart(r, phi, theta)
    expected_x = np.array([1.0, 0.0])
    expected_y = np.array([0.0, 2.0])
    expected_z = np.array([0.0, 0.0])
    assert np.allclose(x, expected_x)
    assert np.allclose(y, expected_y)
    assert np.allclose(z, expected_z)


def test_center_offset():
    r = np.array([1.0])
    phi = np.array([0.0])
    theta = np.array([np.pi / 2])
    center = [1.0, 2.0, 3.0]
    x, y, z = skysegmentor.sphere2cart(r, phi, theta, center=center)
    assert np.isclose(x, 1.0 + 1.0)  # 1 + center_x
    assert np.isclose(y, 0.0 + 2.0)  # 0 + center_y
    assert np.isclose(z, 0.0 + 3.0)  # 0 + center_z


def test_zero_radius():
    r = np.array([0.0])
    phi = np.array([0.0])
    theta = np.array([0.0])
    x, y, z = skysegmentor.sphere2cart(r, phi, theta)
    assert np.isclose(x, 0.0)
    assert np.isclose(y, 0.0)
    assert np.isclose(z, 0.0)


def test_phi_range():
    # phi = 2*pi should wrap around to phi = 0
    r = np.array([1.0])
    phi_0 = np.array([0.0])
    phi_2pi = np.array([2 * np.pi])
    theta = np.array([np.pi / 2])
    x0, y0, z0 = skysegmentor.sphere2cart(r, phi_0, theta)
    x2pi, y2pi, z2pi = skysegmentor.sphere2cart(r, phi_2pi, theta)
    assert np.isclose(x0, x2pi)
    assert np.isclose(y0, y2pi)
    assert np.isclose(z0, z2pi)


def test_theta_range():
    # theta at boundaries 0 and pi
    r = np.array([1.0])
    phi = np.array([np.pi / 4])  # arbitrary
    theta_0 = np.array([0.0])
    theta_pi = np.array([np.pi])
    x0, y0, z0 = skysegmentor.sphere2cart(r, phi, theta_0)
    xpi, ypi, zpi = skysegmentor.sphere2cart(r, phi, theta_pi)
    # at poles x,y = 0
    assert np.isclose(x0, 0.0)
    assert np.isclose(y0, 0.0)
    assert np.isclose(z0, 1.0)
    assert np.isclose(xpi, 0.0)
    assert np.isclose(ypi, 0.0)
    assert np.isclose(zpi, -1.0)


def test_zero_distance_scalar():
    # Distance between same point should be 0
    phi = 0.5
    theta = 1.0
    dist = skysegmentor.distusphere(phi, theta, phi, theta)
    assert np.isclose(dist, 0.0)


def test_opposite_points_scalar():
    # Distance between opposite points should be pi
    phi1 = 0.0
    theta1 = np.pi / 2
    phi2 = np.pi
    theta2 = np.pi / 2
    dist = skysegmentor.distusphere(phi1, theta1, phi2, theta2)
    assert np.isclose(dist, np.pi, atol=1e-14)


def test_orthogonal_points_scalar():
    # Points 90 degrees apart on equator
    phi1 = 0.0
    theta1 = np.pi / 2
    phi2 = np.pi / 2
    theta2 = np.pi / 2
    dist = skysegmentor.distusphere(phi1, theta1, phi2, theta2)
    assert np.isclose(dist, np.pi / 2, atol=1e-14)


def test_zero_distance_array():
    # Array input with identical points returns zeros array
    phi = np.array([0.0, 0.1, 1.0])
    theta = np.array([np.pi / 2, np.pi / 4, 0.5])
    dist = skysegmentor.distusphere(phi, theta, phi, theta)
    assert np.allclose(dist, 0.0)


def test_mixed_array_scalar():
    # One point scalar, one point array - expect vector of distances
    phi1 = 0.0
    theta1 = np.pi / 2
    phi2 = np.array([0.0, np.pi / 2, np.pi])
    theta2 = np.array([np.pi / 2, np.pi / 2, np.pi / 2])
    dist = skysegmentor.distusphere(phi1, theta1, phi2, theta2)
    expected = np.array([0.0, np.pi / 2, np.pi])
    assert np.allclose(dist, expected, atol=1e-14)


def test_distance_symmetry():
    # distusphere should be symmetric: dist(a,b) == dist(b,a)
    phi1 = np.array([0.1, 1.5])
    theta1 = np.array([np.pi / 4, np.pi / 3])
    phi2 = np.array([1.0, 2.0])
    theta2 = np.array([np.pi / 3, np.pi / 6])
    dist_ab = skysegmentor.distusphere(phi1, theta1, phi2, theta2)
    dist_ba = skysegmentor.distusphere(phi2, theta2, phi1, theta1)
    assert np.allclose(dist_ab, dist_ba)


def test_distance_with_phi_wraparound():
    # Distance between points where phi differs by 2*pi should be zero
    phi1 = 0.0
    theta1 = np.pi / 2
    phi2 = 2 * np.pi
    theta2 = np.pi / 2
    dist = skysegmentor.distusphere(phi1, theta1, phi2, theta2)
    assert np.isclose(dist, 0.0, atol=1e-14)


def test_distance_poles():
    # Distance from north pole to south pole
    phi1 = 0.0
    theta1 = 0.0  # north pole
    phi2 = 0.0
    theta2 = np.pi  # south pole
    dist = skysegmentor.distusphere(phi1, theta1, phi2, theta2)
    assert np.isclose(dist, np.pi, atol=1e-14)