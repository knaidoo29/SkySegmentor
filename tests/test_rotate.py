import numpy as np
import pytest
import skysegmentor


def test_rotmat_x_shape_and_type():
    angle = 0.5
    R = skysegmentor._rotmat_x(angle)
    assert isinstance(R, np.ndarray)
    assert R.shape == (9,)

def test_rotmat_y_shape_and_type():
    angle = 0.5
    R = skysegmentor._rotmat_y(angle)
    assert isinstance(R, np.ndarray)
    assert R.shape == (9,)

def test_rotmat_z_shape_and_type():
    angle = 0.5
    R = skysegmentor._rotmat_z(angle)
    assert isinstance(R, np.ndarray)
    assert R.shape == (9,)

def test_rotmat_x_zero_angle():
    R = skysegmentor._rotmat_x(0.0)
    expected = np.array([
        1., 0., 0.,
        0., 1., 0.,
        0., 0., 1.
    ])
    assert np.allclose(R, expected)

def test_rotmat_y_zero_angle():
    R = skysegmentor._rotmat_y(0.0)
    expected = np.array([
        1., 0., 0.,
        0., 1., 0.,
        0., 0., 1.
    ])
    assert np.allclose(R, expected)

def test_rotmat_z_zero_angle():
    R = skysegmentor._rotmat_z(0.0)
    expected = np.array([
        1., 0., 0.,
        0., 1., 0.,
        0., 0., 1.
    ])
    assert np.allclose(R, expected)

def test_rotmat_x_90_degrees():
    angle = np.pi / 2
    R = skysegmentor._rotmat_x(angle)
    expected = np.array([
        1., 0., 0.,
        0., 0., -1.,
        0., 1., 0.
    ])
    assert np.allclose(R, expected, atol=1e-8)

def test_rotmat_y_90_degrees():
    angle = np.pi / 2
    R = skysegmentor._rotmat_y(angle)
    expected = np.array([
        0., 0., 1.,
        0., 1., 0.,
        -1., 0., 0.
    ])
    assert np.allclose(R, expected, atol=1e-8)

def test_rotmat_z_90_degrees():
    angle = np.pi / 2
    R = skysegmentor._rotmat_z(angle)
    expected = np.array([
        0., -1., 0.,
        1., 0., 0.,
        0., 0., 1.
    ])
    assert np.allclose(R, expected, atol=1e-8)

@pytest.mark.parametrize("axis, expected_func", [
    (0, skysegmentor._rotmat_x),
    (1, skysegmentor._rotmat_y),
    (2, skysegmentor._rotmat_z),
])
def test_rotmat_axis_dispatch(axis, expected_func):
    angle = 0.5
    R_axis = skysegmentor._rotmat_axis(angle, axis)
    R_expected = expected_func(angle)
    assert np.allclose(R_axis, R_expected)


def test_rotmat_euler_identity():
    # zero angles -> identity matrix
    angles = [0.0, 0.0, 0.0]
    axes = [0, 1, 2]  # arbitrary axes order
    R = skysegmentor._rotmat_euler(angles, axes)
    I = np.array([1, 0, 0, 0, 1, 0, 0, 0, 1])
    assert np.allclose(R, I), "Euler rotation with zero angles should be identity matrix"

def test_rotmat_euler_rotation_order():
    # Test if multiplication order matters
    angles = [np.pi/2, 0.0, 0.0]
    axes = [0, 1, 2]
    R = skysegmentor._rotmat_euler(angles, axes)
    R_manual = skysegmentor.matrix_dot_3by3(
        skysegmentor.matrix_dot_3by3(
            skysegmentor._rotmat_axis(angles[2], axes[2]), skysegmentor._rotmat_axis(angles[1], axes[1])
        ), skysegmentor._rotmat_axis(angles[0], axes[0])
    )
    assert np.allclose(R, R_manual), "Euler rotation matrix multiplication mismatch"

def test_rotate_3d_identity():
    # Rotation matrix is identity; output should equal input
    x = np.array([1.0, 0.0, 0.0])
    y = np.array([0.0, 1.0, 0.0])
    z = np.array([0.0, 0.0, 1.0])
    I = np.array([1, 0, 0, 0, 1, 0, 0, 0, 1])
    xrot, yrot, zrot = skysegmentor._rotate_3d(x, y, z, I)
    assert np.allclose(xrot, x)
    assert np.allclose(yrot, y)
    assert np.allclose(zrot, z)

def test_rotate_3d_90deg_z_axis():
    # Rotate point (1,0,0) by 90 degrees around z-axis -> should become (0,1,0)
    angle = np.pi/2
    rotz = skysegmentor._rotmat_z(angle)
    x = np.array([1.0])
    y = np.array([0.0])
    z = np.array([0.0])
    xrot, yrot, zrot = skysegmentor._rotate_3d(x, y, z, rotz)
    assert np.allclose(xrot, 0.0, atol=1e-8)
    assert np.allclose(yrot, 1.0, atol=1e-8)
    assert np.allclose(zrot, 0.0, atol=1e-8)

def test_rotate_3d_shape_preservation():
    # Confirm input and output shapes match
    x = np.random.rand(5)
    y = np.random.rand(5)
    z = np.random.rand(5)
    angle = np.pi / 3
    rot = skysegmentor._rotmat_x(angle)
    xrot, yrot, zrot = skysegmentor._rotate_3d(x, y, z, rot)
    assert xrot.shape == x.shape
    assert yrot.shape == y.shape
    assert zrot.shape == z.shape

def test_rotate3d_Euler_identity():
    # Zero angles, rotation around default z-y-z axes, should return original coords
    x, y, z = 1.0, 2.0, 3.0
    angles = np.array([0.0, 0.0, 0.0])
    xrot, yrot, zrot = skysegmentor.rotate3d_Euler(x, y, z, angles)
    assert np.isclose(xrot, x)
    assert np.isclose(yrot, y)
    assert np.isclose(zrot, z)

def test_rotate3d_Euler_known_rotation_z_axis():
    # Rotate point (1,0,0) by 90 degrees around z-axis (single axis)
    x, y, z = np.array([1.0]), np.array([0.0]), np.array([0.0])
    angles = np.array([np.pi/2, 0.0, 0.0])
    axes = "zzz"  # rotate around z three times; only first angle matters here
    xrot, yrot, zrot = skysegmentor.rotate3d_Euler(x, y, z, angles, axes=axes)
    assert np.allclose(xrot, 0.0, atol=1e-8)
    assert np.allclose(yrot, 1.0, atol=1e-8)
    assert np.allclose(zrot, 0.0, atol=1e-8)

def test_rotate3d_Euler_center_shift():
    # Rotation around point other than origin: rotate (1,0,0) by 90 degrees around z, center at (1,0,0)
    x, y, z = 1.0, 0.0, 0.0
    angles = np.array([np.pi/2, 0.0, 0.0])
    center = [1.0, 0.0, 0.0]
    xrot, yrot, zrot = skysegmentor.rotate3d_Euler(x, y, z, angles, axes="zzz", center=center)
    # Since point coincides with center, rotation should leave it unchanged
    assert np.isclose(xrot, 1.0)
    assert np.isclose(yrot, 0.0)
    assert np.isclose(zrot, 0.0)

def test_rotate3d_Euler_array_input():
    # Multiple points rotated by zero angles returns original arrays
    x = np.array([1.0, 0.0, 0.0])
    y = np.array([0.0, 1.0, 0.0])
    z = np.array([0.0, 0.0, 1.0])
    angles = np.array([0.0, 0.0, 0.0])
    xrot, yrot, zrot = skysegmentor.rotate3d_Euler(x, y, z, angles)
    assert np.allclose(xrot, x)
    assert np.allclose(yrot, y)
    assert np.allclose(zrot, z)

def test_rotate3d_Euler_invalid_input_length():
    # Angles and axes length mismatch raises assertion
    x, y, z = 1.0, 1.0, 1.0
    angles = np.array([0.0, 0.0])  # length 2 instead of 3
    axes = "zyz"
    with pytest.raises(AssertionError):
        skysegmentor.rotate3d_Euler(x, y, z, angles, axes)

def test_rotate3d_Euler_axes_conversion():
    # Confirm axes string 'xyz' converts to [0,1,2]
    x, y, z = 1.0, 0.0, 0.0
    angles = np.array([np.pi/2, 0.0, 0.0])
    axes = "xyz"
    xrot, yrot, zrot = skysegmentor.rotate3d_Euler(x, y, z, angles, axes=axes)
    # Rotating 90 deg around x should leave x unchanged, rotate y,z accordingly
    # Here y and z are zero so output yrot, zrot should remain 0
    assert np.isclose(xrot, x)
    assert np.isclose(yrot, 0.0)
    assert np.isclose(zrot, 0.0)

def test_rotate_usphere_identity():
    # Rotating by zero angles should return original angles (within numerical tolerance)
    phi = 0.5
    the = 1.0
    angles = [0.0, 0.0, 0.0]
    phi_rot, the_rot = skysegmentor.rotate_usphere(phi, the, angles)
    assert np.isclose(phi_rot, phi)
    assert np.isclose(the_rot, the)

def test_rotate_usphere_array_identity():
    # Arrays of angles rotated by zero produce same outputs
    phi = np.array([0.1, 0.2, 0.3])
    the = np.array([0.4, 0.5, 0.6])
    angles = [0.0, 0.0, 0.0]
    phi_rot, the_rot = skysegmentor.rotate_usphere(phi, the, angles)
    assert np.allclose(phi_rot, phi)
    assert np.allclose(the_rot, the)

def test_rotate_usphere_known_rotation():
    # Rotate point at phi=0, the=pi/2 by 90 deg about z (should shift phi)
    phi = 0.0
    the = np.pi / 2
    angles = [np.pi/2, 0.0, 0.0]  # rotate 90 deg about first z-axis
    phi_rot, the_rot = skysegmentor.rotate_usphere(phi, the, angles)
    # After rotation, phi should be pi/2 (approx) and the remains pi/2
    assert np.isclose(the_rot, the, atol=1e-8)
    assert np.isclose(phi_rot, np.pi/2, atol=1e-8)

def test_rotate_usphere_output_types_and_shapes():
    # Scalar input returns scalars
    phi, the = 0.1, 0.2
    angles = [0.1, 0.2, 0.3]
    phi_rot, the_rot = skysegmentor.rotate_usphere(phi, the, angles)
    assert isinstance(phi_rot, float)
    assert isinstance(the_rot, float)
    
    # Array input returns arrays of same shape
    phi_arr = np.array([0.1, 0.2])
    the_arr = np.array([0.3, 0.4])
    phi_rot_arr, the_rot_arr = skysegmentor.rotate_usphere(phi_arr, the_arr, angles)
    assert isinstance(phi_rot_arr, np.ndarray)
    assert isinstance(the_rot_arr, np.ndarray)
    assert phi_rot_arr.shape == phi_arr.shape
    assert the_rot_arr.shape == the_arr.shape

def test_midpoint_usphere_basic():
    # Midpoint of identical points is the point itself
    phi, the = 0.5, 1.0
    midphi, midthe = skysegmentor.midpoint_usphere(phi, phi, the, the)
    assert np.isclose(midphi, phi)
    assert np.isclose(midthe, the)

def test_midpoint_usphere_antipodal():
    # Midpoint between opposite points on the unit sphere
    phi1, the1 = 0.0, 0.0
    phi2, the2 = np.pi, 0.0
    midphi, midthe = skysegmentor.midpoint_usphere(phi1, phi2, the1, the2)
    # Midpoint of antipodal points on sphere should be ambiguous but test it returns valid floats
    assert isinstance(midphi, float)
    assert isinstance(midthe, float)
    # The midpoint should lie roughly on equator (the=0), but longitude can be anything
    assert np.isclose(midthe, 0.0, atol=1e-8)

def test_rotate2plane_angles_and_properties():
    # Provide two points and verify returned rotations are lists of floats with length 3
    c1 = [0.0, 0.0]
    c2 = [np.pi/2, np.pi/4]
    a1, a2, a3 = skysegmentor.rotate2plane(c1, c2)
    for angles in [a1, a2, a3]:
        assert isinstance(angles, np.ndarray)
        assert len(angles) == 3
        assert all(isinstance(x, (float, np.floating)) for x in angles)
    
    # After all rotations, points c1 and c2 should lie approximately on equator and midpoint longitude pi
    # We can test this by applying rotations manually and checking
    phi1, the1 = c1
    phi2, the2 = c2
    
    # Apply first rotation
    phi1r, the1r = skysegmentor.rotate_usphere(phi1, the1, a1)
    phi2r, the2r = skysegmentor.rotate_usphere(phi2, the2, a1)
    
    # Apply second rotation
    phi1r, the1r = skysegmentor.rotate_usphere(phi1r, the1r, a2)
    phi2r, the2r = skysegmentor.rotate_usphere(phi2r, the2r, a2)
    
    # Apply third rotation
    phi1r, the1r = skysegmentor.rotate_usphere(phi1r, the1r, a3)
    phi2r, the2r = skysegmentor.rotate_usphere(phi2r, the2r, a3)
    
    # Check latitudes are approximately pi/2 (equator)
    assert np.isclose(the1r, np.pi / 2, atol=1e-6)
    assert np.isclose(the2r, np.pi / 2, atol=1e-6)
    
    # Check midpoint longitude approximately pi
    midphi, midthe = skysegmentor.midpoint_usphere(phi1r, phi2r, the1r, the2r)
    assert np.isclose(midphi, np.pi, atol=1e-6)


def test_forward_backward_rotate_scalar():
    phi = 0.5
    the = 1.0
    a1 = [0.1, 0.2, 0.3]
    a2 = [0.4, 0.5, 0.6]
    a3 = [0.7, 0.8, 0.9]
    
    phi_rot, the_rot = skysegmentor.forward_rotate(phi, the, a1, a2, a3)
    phi_back, the_back = skysegmentor.backward_rotate(phi_rot, the_rot, a1, a2, a3)
    
    assert np.isclose(phi, phi_back, atol=1e-8)
    assert np.isclose(the, the_back, atol=1e-8)

def test_forward_backward_rotate_array():
    phi = np.array([0.1, 0.5, 1.0])
    the = np.array([0.2, 0.6, 1.1])
    a1 = [0.1, 0.2, 0.3]
    a2 = [0.4, 0.5, 0.6]
    a3 = [0.7, 0.8, 0.9]
    
    phi_rot, the_rot = skysegmentor.forward_rotate(phi, the, a1, a2, a3)
    phi_back, the_back = skysegmentor.backward_rotate(phi_rot, the_rot, a1, a2, a3)
    
    assert phi_back.shape == phi.shape
    assert the_back.shape == the.shape
    assert np.allclose(phi, phi_back, atol=1e-8)
    assert np.allclose(the, the_back, atol=1e-8)

def test_output_types_and_shapes():
    phi = np.array([0.3, 0.7])
    the = np.array([0.4, 0.8])
    a1 = [0.0, 0.0, 0.0]
    a2 = [0.0, 0.0, 0.0]
    a3 = [0.0, 0.0, 0.0]
    
    phi_rot, the_rot = skysegmentor.forward_rotate(phi, the, a1, a2, a3)
    assert isinstance(phi_rot, np.ndarray)
    assert isinstance(the_rot, np.ndarray)
    assert phi_rot.shape == phi.shape
    assert the_rot.shape == the.shape
    
    phi_rot, the_rot = skysegmentor.forward_rotate(0.5, 0.5, a1, a2, a3)
    assert isinstance(phi_rot, float) or isinstance(phi_rot, np.floating)
    assert isinstance(the_rot, float) or isinstance(the_rot, np.floating)