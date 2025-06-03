import numpy as np
import healpy as hp
import skysegmentor
import pytest

def test_get_partition_IDs_basic():
    partition = np.array([0, 1, 1, 2, 2, 0])
    expected = np.array([0, 1, 2])
    result = skysegmentor.get_partition_IDs(partition)
    assert np.array_equal(result, expected)

def test_get_partition_IDs_no_zero():
    partition = np.array([1, 3, 3, 2, 2, 1])
    expected = np.array([1, 2, 3])
    result = skysegmentor.get_partition_IDs(partition)
    assert np.array_equal(result, expected)

def test_get_partition_IDs_single_value():
    partition = np.array([0, 0, 0])
    expected = np.array([0])
    result = skysegmentor.get_partition_IDs(partition)
    assert np.array_equal(result, expected)

def test_get_partition_IDs_empty():
    partition = np.array([])
    expected = np.array([])
    result = skysegmentor.get_partition_IDs(partition)
    assert np.array_equal(result, expected)

def test_get_partition_IDs_large_array():
    partition = np.array([5, 2, 3, 5, 0, 2, 2, 3])
    expected = np.array([0, 2, 3, 5])
    result = skysegmentor.get_partition_IDs(partition)
    assert np.array_equal(result, expected)

def test_total_partition_weights_basic():
    partition = np.array([0, 1, 1, 2, 2, 0])
    weights = np.array([1, 2, 3, 4, 5, 6])
    expected_ids = np.array([0, 1, 2])
    expected_weights = np.array([1 + 6, 2 + 3, 4 + 5])
    result_ids, result_weights = skysegmentor.total_partition_weights(partition, weights)
    assert np.array_equal(result_ids, expected_ids)
    assert np.allclose(result_weights, expected_weights)

def test_total_partition_weights_single_partition():
    partition = np.array([3, 3, 3])
    weights = np.array([1, 2, 3])
    expected_ids = np.array([3])
    expected_weights = np.zeros(4)  # max ID + 1 = 4
    expected_weights[3] = 6
    result_ids, result_weights = skysegmentor.total_partition_weights(partition, weights)
    assert np.array_equal(result_ids, expected_ids)
    assert np.allclose(result_weights, expected_weights)

def test_total_partition_weights_with_zeros():
    partition = np.array([0, 0, 1])
    weights = np.array([1.0, 2.0, 3.0])
    expected_ids = np.array([0, 1])
    expected_weights = np.zeros(2)  # max ID + 1 = 2
    expected_weights[0] = 3.0
    expected_weights[1] = 3.0
    result_ids, result_weights = skysegmentor.total_partition_weights(partition, weights)
    assert np.array_equal(result_ids, expected_ids)
    assert np.allclose(result_weights, expected_weights)

def test_total_partition_weights_noncontiguous_ids():
    partition = np.array([0, 2, 5, 2])
    weights = np.array([1.0, 2.0, 3.0, 4.0])
    expected_ids = np.array([0, 2, 5])
    expected_weights = np.zeros(6)  # max ID + 1 = 6
    expected_weights[0] = 1.0
    expected_weights[2] = 6.0
    expected_weights[5] = 3.0
    result_ids, result_weights = skysegmentor.total_partition_weights(partition, weights)
    assert np.array_equal(result_ids, expected_ids)
    assert np.allclose(result_weights, expected_weights)

def test_remove_val4array_basic():
    array = np.array([1, 2, 3, 4, 2])
    val = 2
    result = skysegmentor.remove_val4array(array, val)
    expected = np.array([1, 3, 4])
    assert np.array_equal(result, expected)

def test_remove_val4array_no_match():
    array = np.array([1, 3, 5])
    val = 2
    result = skysegmentor.remove_val4array(array, val)
    expected = array
    assert np.array_equal(result, expected)

def test_remove_val4array_all_match():
    array = np.array([2, 2, 2])
    val = 2
    result = skysegmentor.remove_val4array(array, val)
    expected = np.array([])
    assert np.array_equal(result, expected)

def test_remove_val4array_empty_input():
    array = np.array([])
    val = 1
    result = skysegmentor.remove_val4array(array, val)
    expected = np.array([])
    assert np.array_equal(result, expected)

def test_remove_val4array_float_precision():
    array = np.array([1.0, 2.000000001, 3.0])
    val = 2.0
    result = skysegmentor.remove_val4array(array, val)
    expected = array  # 2.000000001 != 2.0 due to precision
    assert np.array_equal(result, expected)

def test_remove_val4array_negative_values():
    array = np.array([-1, -2, 0, 1])
    val = -2
    result = skysegmentor.remove_val4array(array, val)
    expected = np.array([-1, 0, 1])
    assert np.array_equal(result, expected)

def test_fill_map_basic():
    nside = 4
    npix = hp.nside2npix(nside)
    pixID = np.array([0, 1, 2])
    val = 2.5
    result = skysegmentor.fill_map(pixID, nside, val)
    expected = np.zeros(npix)
    expected[pixID] = val
    assert np.allclose(result, expected)

def test_fill_map_default_value():
    nside = 4
    pixID = np.array([5, 6])
    result = skysegmentor.fill_map(pixID, nside)  # default val=1.0
    expected = np.zeros(hp.nside2npix(nside))
    expected[pixID] = 1.0
    assert np.allclose(result, expected)

def test_fill_map_empty_pixID():
    nside = 4
    pixID = np.array([], dtype=int)
    result = skysegmentor.fill_map(pixID, nside, val=3.0)
    expected = np.zeros(hp.nside2npix(nside))
    assert np.allclose(result, expected)

def test_fill_map_multiple_overwrites():
    nside = 4
    pixID = np.array([1, 1, 1])  # duplicate entries
    result = skysegmentor.fill_map(pixID, nside, val=9.0)
    expected = np.zeros(hp.nside2npix(nside))
    expected[1] = 9.0  # value should still be set
    assert np.allclose(result, expected)

def test_fill_map_invalid_pixID_raises():
    nside = 4
    npix = hp.nside2npix(nside)
    pixID = np.array([npix + 1])  # invalid pixel index
    with pytest.raises(IndexError):
        skysegmentor.fill_map(pixID, nside, val=1.0)

def test_barycenter_single_pixel():
    nside = 4
    npix = hp.nside2npix(nside)
    bnmap = np.zeros(npix)
    bnmap[10] = 1.0
    phic, thec, themax = skysegmentor.find_map_barycenter(bnmap)
    assert isinstance(phic, float)
    assert isinstance(thec, float)
    assert isinstance(themax, float)

def test_barycenter_uniform_map():
    nside = 4
    npix = hp.nside2npix(nside)
    bnmap = np.ones(npix)
    phic, thec, themax = skysegmentor.find_map_barycenter(bnmap)
    assert 0 <= phic < 2 * np.pi
    assert 0 <= thec <= np.pi
    assert themax > 0

def test_barycenter_with_weights():
    nside = 4
    npix = hp.nside2npix(nside)
    bnmap = np.zeros(npix)
    bnmap[[5, 10, 15]] = 1.0
    weights = np.zeros(npix)
    weights[5] = 0.1
    weights[10] = 1.0
    weights[15] = 0.5
    phic, thec, themax = skysegmentor.find_map_barycenter(bnmap, weights)
    assert 0 <= phic < 2 * np.pi
    assert 0 <= thec <= np.pi
    assert themax > 0

def test_barycenter_empty_map():
    nside = 4
    npix = hp.nside2npix(nside)
    bnmap = np.zeros(npix)
    with pytest.raises(ValueError):
        # This would cause division by zero or empty max
        skysegmentor.find_map_barycenter(bnmap)

def test_barycenter_partial_filled():
    nside = 4
    npix = hp.nside2npix(nside)
    bnmap = np.zeros(npix)
    bnmap[5] = 1.0
    phic, thec, themax = skysegmentor.find_map_barycenter(bnmap)
    assert isinstance(phic, float)
    assert isinstance(thec, float)
    assert isinstance(themax, float)

def test_find_points_barycenter_basic():
    phi = np.array([0.0, np.pi])
    the = np.array([np.pi / 2, np.pi / 2])
    phic, thec, themax = skysegmentor.find_points_barycenter(phi, the)
    
    assert np.isclose(phic, np.pi / 2)
    assert np.isclose(thec, np.pi / 4)
    assert np.isclose(themax, np.pi / 2)


def test_find_points_barycenter_with_weights():
    phi = np.array([0.0, np.pi])
    the = np.array([np.pi / 2, np.pi / 2])
    weights = np.array([1.0, 3.0])
    phic, thec, themax = skysegmentor.find_points_barycenter(phi, the, weights)
    
    assert 0 <= phic <= 2 * np.pi
    assert 0 <= thec <= np.pi
    assert np.isclose(themax, np.pi)


def test_find_points_barycenter_empty_input():
    with pytest.raises(ValueError, match="Input arrays phi and the must not be empty."):
        skysegmentor.find_points_barycenter(np.array([]), np.array([]))


def test_find_points_barycenter_mismatched_input_lengths():
    phi = np.array([0.0, np.pi])
    the = np.array([np.pi / 2])
    with pytest.raises(ValueError, match="Input arrays phi and the must have the same length."):
        skysegmentor.find_points_barycenter(phi, the)


def test_find_points_barycenter_mismatched_weights_length():
    phi = np.array([0.0, np.pi])
    the = np.array([np.pi / 2, np.pi / 2])
    weights = np.array([1.0])
    with pytest.raises(ValueError, match="Weights array must be the same length as phi and the."):
       skysegmentor. find_points_barycenter(phi, the, weights)


def test_get_map_border_basic():
    nside = 8
    npix = hp.nside2npix(nside)

    # Create a circular region around theta=π/2, phi=π
    bnmap = np.zeros(npix)
    the_center = np.pi / 2
    phi_center = np.pi
    vec = hp.ang2vec(the_center, phi_center)
    disk_pix = hp.query_disc(nside, vec, radius=np.radians(10))
    bnmap[disk_pix] = 1.0

    phi_border, the_border = skysegmentor.get_map_border(bnmap)

    assert isinstance(phi_border, np.ndarray)
    assert isinstance(the_border, np.ndarray)
    assert phi_border.shape == the_border.shape
    assert phi_border.ndim == 1
    assert np.all((phi_border >= 0) & (phi_border <= 2 * np.pi))
    assert np.all((the_border >= 0) & (the_border <= np.pi))


def test_get_map_border_with_weights():
    nside = 8
    npix = hp.nside2npix(nside)
    bnmap = np.zeros(npix)
    weights = np.ones(npix)

    disk_pix = hp.query_disc(nside, hp.ang2vec(np.pi / 2, 0), np.radians(8))
    bnmap[disk_pix] = 1.0
    weights[disk_pix] = np.linspace(1, 10, len(disk_pix))

    phi_border, the_border = skysegmentor.get_map_border(bnmap, wmap=weights)

    assert isinstance(phi_border, np.ndarray)
    assert isinstance(the_border, np.ndarray)
    assert phi_border.shape == the_border.shape
    assert len(phi_border) > 0


def test_get_map_border_custom_resolution():
    nside = 8
    npix = hp.nside2npix(nside)
    bnmap = np.zeros(npix)

    disk_pix = hp.query_disc(nside, hp.ang2vec(np.pi / 2, np.pi / 2), np.radians(5))
    bnmap[disk_pix] = 1.0

    phi_border_low, the_border_low = skysegmentor.get_map_border(bnmap, res=[50, 25])
    phi_border_high, the_border_high = skysegmentor.get_map_border(bnmap, res=[200, 100])

    # Higher resolution should return more detailed border (not necessarily more points)
    assert len(phi_border_low) > 0
    assert len(phi_border_high) > 0


def test_get_map_border_empty_map():
    nside = 8
    npix = hp.nside2npix(nside)
    bnmap = np.zeros(npix)

    with pytest.raises(ValueError, match="Binary map must contain at least one non-zero pixel."):
        skysegmentor.get_map_border(bnmap)


def generate_circle_patch(center_phi, center_the, radius, num_points=200):
    """Generate points in a circular patch on the sphere."""
    phis = np.linspace(0, 2 * np.pi, num_points)
    thetas = np.full_like(phis, radius)
    # Convert to Cartesian
    x, y, z = np.sin(thetas) * np.cos(phis), np.sin(thetas) * np.sin(phis), np.cos(thetas)
    # Rotate to desired center
    from scipy.spatial.transform import Rotation as R

    rot = R.from_euler("yz", [-center_the, -center_phi])
    coords = rot.apply(np.vstack([x, y, z]).T)

    x_new, y_new, z_new = coords.T
    r = np.sqrt(x_new**2 + y_new**2 + z_new**2)
    theta = np.arccos(z_new / r)
    phi = np.arctan2(y_new, x_new) % (2 * np.pi)
    return phi, theta


def test_get_points_border_basic():
    phi, the = generate_circle_patch(np.pi, np.pi / 2, np.radians(10), num_points=300)
    phi_border, the_border = skysegmentor.get_points_border(phi, the)

    assert isinstance(phi_border, np.ndarray)
    assert isinstance(the_border, np.ndarray)
    assert phi_border.shape == the_border.shape
    assert phi_border.ndim == 1
    assert len(phi_border) > 0
    assert np.all((phi_border >= 0) & (phi_border <= 2 * np.pi))
    assert np.all((the_border >= 0) & (the_border <= np.pi))


def test_get_points_border_with_weights():
    phi, the = generate_circle_patch(0.0, np.pi / 3, np.radians(8), num_points=150)
    weights = np.linspace(1, 2, len(phi))

    phi_border, the_border = skysegmentor.get_points_border(phi, the, weights=weights)

    assert isinstance(phi_border, np.ndarray)
    assert isinstance(the_border, np.ndarray)
    assert len(phi_border) > 0


def test_get_points_border_custom_resolution():
    phi, the = generate_circle_patch(np.pi / 2, np.pi / 2, np.radians(15), num_points=500)

    border_low_res = skysegmentor.get_points_border(phi, the, res=20)
    border_high_res = skysegmentor.get_points_border(phi, the, res=200)

    assert len(border_low_res[0]) < len(border_high_res[0])


def test_get_points_border_empty():
    phi = np.array([])
    the = np.array([])

    with pytest.raises(ValueError, match="Input point set is empty"):
        skysegmentor.get_points_border(phi, the)


def generate_disk_mask(nside, center, radius_deg):
    """Create a binary map with a circular disk mask."""
    npix = hp.nside2npix(nside)
    vec = hp.ang2vec(center[1], center[0])  # healpy uses theta, phi
    ipix = hp.query_disc(nside, vec, np.radians(radius_deg))
    mask = np.zeros(npix)
    mask[ipix] = 1.0
    return mask


def test_get_map_most_dist_points_basic():
    nside = 32
    center = (np.pi, np.pi / 2)  # (phi, theta)
    radius_deg = 10
    bnmap = generate_disk_mask(nside, center, radius_deg)

    p1, t1, p2, t2 = skysegmentor.get_map_most_dist_points(bnmap)

    assert isinstance(p1, float)
    assert isinstance(t1, float)
    assert isinstance(p2, float)
    assert isinstance(t2, float)

    # Check that outputs are within bounds
    assert 0 <= p1 <= 2 * np.pi
    assert 0 <= p2 <= 2 * np.pi
    assert 0 <= t1 <= np.pi
    assert 0 <= t2 <= np.pi

    # Should be close to twice the radius (approx. max separation on the border)
    max_sep = skysegmentor.distusphere(np.array([p1]), np.array([t1]), np.array([p2]), np.array([t2]))[0]
    assert np.isclose(max_sep, 2 * np.radians(radius_deg), atol=0.1)


def test_get_map_most_dist_points_with_weights():
    nside = 16
    center = (np.pi / 2, np.pi / 3)
    radius_deg = 5
    bnmap = generate_disk_mask(nside, center, radius_deg)
    wmap = np.random.rand(len(bnmap))

    p1, t1, p2, t2 = skysegmentor.get_map_most_dist_points(bnmap, wmap=wmap)
    assert isinstance(p1, float) and isinstance(p2, float)


def test_get_map_most_dist_points_custom_res():
    nside = 16
    bnmap = generate_disk_mask(nside, (np.pi, np.pi / 2), 8)
    p1_low, t1_low, p2_low, t2_low = skysegmentor.get_map_most_dist_points(bnmap, res=[10, 5])
    p1_high, t1_high, p2_high, t2_high = skysegmentor.get_map_most_dist_points(bnmap, res=[200, 100])

    # Not checking values, just ensuring execution and that higher res works
    assert all(isinstance(x, float) for x in [p1_high, t1_high, p2_high, t2_high])


def test_get_map_most_dist_points_empty_map():
    nside = 16
    bnmap = np.zeros(hp.nside2npix(nside))

    with pytest.raises(ValueError, match="Binary map must contain at least one non-zero pixel."):
        skysegmentor.get_map_most_dist_points(bnmap)


def generate_circle_points(center_phi, center_the, radius_deg, n_points=500):
    """Generate points roughly on a spherical circle."""
    radius_rad = np.radians(radius_deg)
    phis = np.linspace(0, 2 * np.pi, n_points)
    thes = np.full(n_points, center_the - radius_rad)
    return phis, thes


def test_get_points_most_dist_points_basic_circle():
    phi, the = generate_circle_points(np.pi, np.pi / 2, 10)
    p1, t1, p2, t2 = skysegmentor.get_points_most_dist_points(phi, the)

    assert all(isinstance(x, float) for x in [p1, t1, p2, t2])
    assert 0 <= p1 <= 2 * np.pi
    assert 0 <= p2 <= 2 * np.pi
    assert 0 <= t1 <= np.pi
    assert 0 <= t2 <= np.pi

    max_dist = skysegmentor.distusphere(np.array([p1]), np.array([t1]), np.array([p2]), np.array([t2]))[0]
    assert np.isclose(max_dist, 2.7920213882461278, atol=0.1)


def test_get_points_most_dist_points_with_weights():
    phi, the = generate_circle_points(np.pi / 2, np.pi / 3, 15)
    weights = np.random.rand(len(phi))
    p1, t1, p2, t2 = skysegmentor.get_points_most_dist_points(phi, the, weights=weights)
    assert all(isinstance(x, float) for x in [p1, t1, p2, t2])


def test_get_points_most_dist_points_high_resolution():
    phi, the = generate_circle_points(np.pi, np.pi / 4, 8, n_points=1000)
    p1_low, t1_low, p2_low, t2_low = skysegmentor.get_points_most_dist_points(phi, the, res=10)
    p1_high, t1_high, p2_high, t2_high = skysegmentor.get_points_most_dist_points(phi, the, res=200)
    # Ensure execution and output types
    assert isinstance(p1_high, float) and isinstance(p2_high, float)


def test_get_points_most_dist_points_empty_input():
    phi = np.array([])
    the = np.array([])
    with pytest.raises(ValueError, match="Input coordinate arrays are empty."):
        skysegmentor.get_points_most_dist_points(phi, the)

def test_weight_dif_balanced():
    phi = np.array([0.1, 1.0, 2.0, 3.0])
    weights = np.array([1.0, 1.0, 1.0, 1.0])
    phi_split = 1.5  # Splits into [0.1, 1.0] and [2.0, 3.0]
    assert skysegmentor.weight_dif(phi_split, phi, weights) == 0.0

def test_weight_dif_unbalanced():
    phi = np.array([0.1, 1.0, 2.0, 3.0])
    weights = np.array([1.0, 2.0, 3.0, 4.0])
    phi_split = 2.0  # Splits into [0.1, 1.0, 2.0] and [3.0]
    expected = abs((1 + 2 + 3) - 4)
    assert skysegmentor.weight_dif(phi_split, phi, weights) == expected

def test_weight_dif_with_balance():
    phi = np.array([0.1, 1.0, 2.0, 3.0])
    weights = np.array([1.0, 2.0, 3.0, 4.0])
    phi_split = 2.0
    balance = 0.5
    expected = abs((1 + 2 + 3) * balance - 4)
    assert np.isclose(skysegmentor.weight_dif(phi_split, phi, weights, balance), expected)

def test_weight_dif_all_left():
    phi = np.array([0.1, 0.5, 1.0])
    weights = np.array([1.0, 2.0, 3.0])
    phi_split = 2.0
    expected = abs((1 + 2 + 3) - 0)
    assert skysegmentor.weight_dif(phi_split, phi, weights) == expected

def test_weight_dif_all_right():
    phi = np.array([2.1, 2.5, 3.0])
    weights = np.array([1.0, 2.0, 3.0])
    phi_split = 2.0
    expected = abs(0 - (1 + 2 + 3))
    assert skysegmentor.weight_dif(phi_split, phi, weights) == expected

def test_find_dphi_basic():
    phi = np.linspace(0., 4., 100)
    weights = np.ones(100)
    dphi = skysegmentor.find_dphi(phi, weights)
    # Since weights are equal, split roughly in the middle of the range
    assert isinstance(dphi, float)
    assert 1.9 < dphi < 2.1

def test_find_dphi_unbalanced():
    phi = np.linspace(0., 4., 100)
    weights = np.linspace(0., 4., 100)
    balance = 2
    dphi = skysegmentor.find_dphi(phi, weights, balance=balance)
    assert isinstance(dphi, float)
    # The splitting point should shift to balance the weighted sums
    # So it should be greater than the midpoint
    assert dphi > 2

def test_find_dphi_single_point():
    phi = np.array([5])
    weights = np.array([10])
    dphi = skysegmentor.find_dphi(phi, weights)
    # If only one point, dphi should be exactly that point
    assert dphi == 5

def test_find_dphi_zero_weights():
    phi = np.array([0, 1, 2, 3, 4])
    weights = np.zeros_like(phi)
    with pytest.raises(ValueError, match="Weights must contain at least one non-zero value."):
        skysegmentor.find_dphi(phi, weights)

def test_find_dphi_non_uniform_phi():
    phi = np.array([0, 10, 20, 30, 40])
    weights = np.array([1, 2, 3, 4, 5])
    dphi = skysegmentor.find_dphi(phi, weights)
    assert isinstance(dphi, float)
    assert 19 <= dphi <= 21

def test_segmentmap2_basic_partition():
    nside = 8  # low res for test, npix=768
    npix = hp.nside2npix(nside)
    weightmap = np.ones(npix)
    part_map = skysegmentor.segmentmap2(weightmap)
    # partitionmap should only contain 1 or 2 (two partitions)
    unique_parts = np.unique(part_map)
    assert set(unique_parts).issubset({0, 1, 2})
    assert 1 in unique_parts  # original partition present
    assert 2 in unique_parts  # new partition created

def test_segmentmap2_with_partitionmap():
    nside = 8
    npix = hp.nside2npix(nside)
    weightmap = np.ones(npix)
    partitionmap = np.zeros(npix)
    partitionmap[:npix//2] = 1
    partitionmap[npix//2:] = 2
    # Partition only the partition=1
    new_part_map = skysegmentor.segmentmap2(weightmap, partitionmap=partitionmap, partition=1)
    # The partition 1 should be split further (new partition number = max + 1 = 3)
    assert np.any(new_part_map == 3)
    # Other partitions remain unchanged
    assert np.all(new_part_map[partitionmap == 2] == 2)

def test_segmentmap2_balance_effect():
    nside = 8
    npix = hp.nside2npix(nside)
    weightmap = np.linspace(1, 10, npix)
    part_map = skysegmentor.segmentmap2(weightmap, balance=2)
    assert np.any(part_map > 1)  # ensure partitioning happened

def test_segmentmap2_no_weights():
    nside = 8
    npix = hp.nside2npix(nside)
    weightmap = np.zeros(npix)
    with pytest.raises(ValueError, match="Binary map must contain at least one non-zero pixel."):
        skysegmentor.segmentmap2(weightmap)

def test_segmentmap2_single_pixel():
    nside = 1
    npix = hp.nside2npix(nside)  # 12 pixels for nside=1
    weightmap = np.zeros(npix)
    weightmap[5] = 1
    part_map = skysegmentor.segmentmap2(weightmap)
    # partitionmap should be 1 where weight present
    assert part_map[5] == 1
    assert np.sum(part_map) == 1  # only one pixel assigned

def test_segmentpoints2_basic():
    n = 10
    phi = np.linspace(0, np.pi, n)
    the = np.linspace(0, np.pi/2, n)
    weights = np.ones(n)
    part_id = skysegmentor.segmentpoints2(phi, the, weights=weights)
    assert isinstance(part_id, np.ndarray)
    assert set(np.unique(part_id)).issubset({1, 2})
    # Check partitionID length unchanged
    assert len(part_id) == n

def test_segmentpoints2_with_partitionID():
    n = 20
    phi = np.linspace(0, 2*np.pi, n)
    the = np.linspace(0, np.pi, n)
    weights = np.linspace(1, 10, n)
    partitionID = np.ones(n)
    # split partition 1 into 1 and 2
    out_partition = skysegmentor.segmentpoints2(phi, the, weights, partitionID=partitionID, partition=1)
    assert 2 in out_partition
    # other partitions unchanged if any
    assert np.all(out_partition[partitionID != 1] == partitionID[partitionID != 1])

def test_segmentpoints2_balance_effect():
    n = 15
    phi = np.linspace(0, np.pi, n)
    the = np.linspace(0, np.pi/3, n)
    weights = np.linspace(1, 5, n)
    part1 = skysegmentor.segmentpoints2(phi, the, weights, balance=1)
    part2 = skysegmentor.segmentpoints2(phi, the, weights, balance=3)
    # Partition results differ due to balance
    assert not np.array_equal(part1, part2)

def test_segmentpoints2_weights_none():
    n = 5
    phi = np.array([0, 1, 2, 3, 4])
    the = np.array([0, 0.5, 1, 1.5, 2])
    # weights None should default to ones
    part_id = skysegmentor.segmentpoints2(phi, the, weights=None)
    assert set(np.unique(part_id)).issubset({1, 2})

def test_segmentpoints2_single_point():
    phi = np.array([1.0])
    the = np.array([0.5])
    weights = np.array([2.0])
    part_id = skysegmentor.segmentpoints2(phi, the, weights)
    # Single point, partitionID should be [1]
    assert np.array_equal(part_id, np.array([1]))

def test_segmentmapN_basic_partitioning():
    nside = 8
    npix = 12 * nside**2
    weightmap = np.ones(npix)
    partitionmap = skysegmentor.segmentmapN(weightmap, Npartitions=4)
    unique_parts = np.unique(partitionmap)
    assert len(unique_parts) == 4
    assert set(unique_parts).issubset({0,1,2,3,4})

def test_segmentmapN_single_partition():
    nside = 4
    npix = 12 * nside**2
    weightmap = np.ones(npix)
    with pytest.raises(ValueError, match="Npartitions must be > 1."):
        skysegmentor.segmentmapN(weightmap, Npartitions=1)

def test_segmentmapN_no_weights():
    nside = 4
    npix = 12 * nside**2
    weightmap = np.zeros(npix)
    with pytest.raises(ValueError, match="Binary map must contain at least one non-zero pixel."):
        skysegmentor.segmentmapN(weightmap, Npartitions=3)

def test_segmentmapN_partition_numbers_increasing():
    nside = 8
    npix = 12 * nside**2
    weightmap = np.ones(npix)
    partitionmap = skysegmentor.segmentmapN(weightmap, Npartitions=8)
    unique_parts = np.unique(partitionmap)
    # Partition IDs should be between 1 and 8 (and possibly 0)
    assert all((unique_parts >= 0) & (unique_parts <= 8))

def test_segmentmapN_large_Npartitions():
    nside = 4
    npix = 12 * nside**2
    weightmap = np.ones(npix)
    partitionmap = skysegmentor.segmentmapN(weightmap, Npartitions=16)
    unique_parts = np.unique(partitionmap)
    assert len(unique_parts) == 16

def test_segmentpointsN_basic_partitioning():
    n = 20
    phi = np.linspace(0, np.pi, n)
    the = np.linspace(0, np.pi/2, n)
    weights = np.ones(n)
    partitionID = skysegmentor.segmentpointsN(phi, the, Npartitions=4, weights=weights)
    unique_parts = np.unique(partitionID)
    assert len(unique_parts) == 4
    assert set(unique_parts).issubset({1, 2, 3, 4})

def test_segmentpointsN_weights_none():
    n = 10
    phi = np.linspace(0, 2*np.pi, n)
    the = np.linspace(0, np.pi, n)
    partitionID = skysegmentor.segmentpointsN(phi, the, Npartitions=3, weights=None)
    unique_parts = np.unique(partitionID)
    assert len(unique_parts) == 3

def test_segmentpointsN_invalid_Npartitions():
    n = 5
    phi = np.zeros(n)
    the = np.zeros(n)
    with pytest.raises(ValueError):
        skysegmentor.segmentpointsN(phi, the, Npartitions=1)

def test_segmentpointsN_partition_numbers_increasing():
    n = 30
    phi = np.linspace(0, 3*np.pi, n)
    the = np.linspace(0, np.pi, n)
    weights = np.linspace(1, 5, n)
    partitionID = skysegmentor.segmentpointsN(phi, the, Npartitions=6, weights=weights)
    unique_parts = np.unique(partitionID)
    # Partition IDs should be between 1 and Npartitions
    assert all((unique_parts >= 1) & (unique_parts <= 6))

def test_segmentpointsN_single_point():
    phi = np.array([1.0])
    the = np.array([0.5])
    weights = np.array([1.0])
    partitionID = skysegmentor.segmentpointsN(phi, the, Npartitions=2, weights=weights)
    # Only one point, partition should remain 1
    assert np.array_equal(partitionID, np.array([1]))