import numpy as np
import healpy as hp
import skysegmentor
import pytest


def test_no_cascade():
    labels = np.array([1, 2, 3, 4, 5])  # labels[i] == i+1
    for i in range(1, 6):
        assert skysegmentor._cascade(labels, i) == i


def test_single_step_cascade():
    labels = np.array([2, 2, 3, 4, 5])  # label 1 -> 2, others unchanged
    assert skysegmentor._cascade(labels, 1) == 2
    assert skysegmentor._cascade(labels, 2) == 2
    assert skysegmentor._cascade(labels, 3) == 3


def test_multi_step_cascade():
    # 1->2, 2->4, 4->4 (stop)
    labels = np.array([2, 4, 3, 4, 5])
    assert skysegmentor._cascade(labels, 1) == 4
    assert skysegmentor._cascade(labels, 2) == 4
    assert skysegmentor._cascade(labels, 3) == 3
    assert skysegmentor._cascade(labels, 4) == 4


def test_cascade_index_out_of_range():
    labels = np.array([1, 2, 3])
    with pytest.raises(IndexError):
        skysegmentor._cascade(labels, 4)  # indexin out of range


def test_cascade_edge_case_empty_labels():
    labels = np.array([])
    with pytest.raises(IndexError):
        skysegmentor._cascade(labels, 1)


def test_cascade_all_no_cascade():
    labels = np.array([1, 2, 3, 4, 5])
    expected = labels.copy()
    result = skysegmentor._cascade_all(labels)
    assert np.array_equal(result, expected)


def test_cascade_all_single_step():
    labels = np.array([2, 2, 3, 4, 5])  # 1->2, others unchanged
    expected = np.array([2, 2, 3, 4, 5])
    result = skysegmentor._cascade_all(labels)
    assert np.array_equal(result, expected)


def test_cascade_all_multi_step():
    labels = np.array([2, 4, 3, 4, 5])  # 1->2->4, 2->4, others unchanged
    expected = np.array([4, 4, 3, 4, 5])
    result = skysegmentor._cascade_all(labels)
    assert np.array_equal(result, expected)


def test_cascade_all_all_cascade_to_same():
    labels = np.array([2, 3, 4, 5, 5])  # chain 1->2->3->4->5->5
    expected = np.array([5, 5, 5, 5, 5])
    result = skysegmentor._cascade_all(labels)
    assert np.array_equal(result, expected)


def test_cascade_all_empty():
    labels = np.array([])
    expected = np.array([])
    result = skysegmentor._cascade_all(labels)
    assert np.array_equal(result, expected)


def test_cascade_all_single_element():
    labels = np.array([1])
    expected = np.array([1])
    result = skysegmentor._cascade_all(labels)
    assert np.array_equal(result, expected)


def test_cascade_all_with_invalid_index():
    labels = np.array([2, 3, 10])  # 10 is out of bounds
    with pytest.raises(IndexError):
        skysegmentor._cascade_all(labels)


def test_unionise_no_cascade():
    labels = np.array([1, 2, 3, 4])
    ind1out, ind2out, indout = skysegmentor._unionise(1, 3, labels)
    assert ind1out == 1
    assert ind2out == 3
    assert indout == 1  # smaller of (1,3)


def test_unionise_with_cascade():
    labels = np.array([2, 2, 3, 4])  # 1->2, 2->2, 3->3, 4->4
    ind1out, ind2out, indout = skysegmentor._unionise(1, 2, labels)
    assert ind1out == 2
    assert ind2out == 2
    assert indout == 2


def test_unionise_reverse_order():
    labels = np.array([1, 2, 3, 3])  # 1->1, 2->2, 3->3, 4->3
    ind1out, ind2out, indout = skysegmentor._unionise(3, 2, labels)
    assert ind1out == 3
    assert ind2out == 2
    assert indout == 2


def test_unionise_equal_indices():
    labels = np.array([1, 2, 3])
    ind1out, ind2out, indout = skysegmentor._unionise(2, 2, labels)
    assert ind1out == 2
    assert ind2out == 2
    assert indout == 2


def test_unionise_invalid_index():
    labels = np.array([1, 2, 3])
    with pytest.raises(IndexError):
        skysegmentor._unionise(1, 5, labels)  # 5 is out of bounds


def test_shuffle_down_basic():
    labels = np.array([1, 2, 3, 4])
    result = skysegmentor._shuffle_down(labels)
    expected = np.array([1, 2, 3, 4])
    assert np.array_equal(result, expected)


def test_shuffle_down_with_duplicates():
    labels = np.array([2, 2, 3, 3, 5])
    # labels 1 and 4 are missing, so labels get remapped compactly:
    # unique labels present are 2,3,5 which map to 1,2,3 respectively
    expected = np.array([1, 1, 2, 2, 3])
    result = skysegmentor._shuffle_down(labels)
    assert np.array_equal(result, expected)


def test_shuffle_down_all_same_label():
    labels = np.array([3, 3, 3])
    expected = np.array([1, 1, 1])
    result = skysegmentor._shuffle_down(labels)
    assert np.array_equal(result, expected)


def test_shuffle_down_non_consecutive_labels():
    labels = np.array([3, 4, 3, 4, 5])
    # unique labels present: 3,4,5 → mapped to 1,2,3
    expected = np.array([1, 2, 1, 2, 3])
    result = skysegmentor._shuffle_down(labels)
    assert np.array_equal(result, expected)


def test_shuffle_down_empty():
    labels = np.array([])
    expected = np.array([])
    result = skysegmentor._shuffle_down(labels)
    assert np.array_equal(result, expected)


def test_shuffle_down_single_element():
    labels = np.array([1])
    expected = np.array([1])
    result = skysegmentor._shuffle_down(labels)
    assert np.array_equal(result, expected)


def test_shuffle_down_invalid_labels():
    labels = np.array([0, -1, 3])
    with pytest.raises(AssertionError):
        skysegmentor._shuffle_down(labels)


def test_single_element_returns_same():
    arr = np.array([1, 2, 3])
    result = skysegmentor._if_list_concatenate([np.array([1]), np.array([2, 3])])
    # Should return the input array unchanged (not wrapped in another array)
    assert np.array_equal(result, arr)


def test_multiple_elements_concatenated():
    arr1 = np.array([1, 2])
    arr2 = np.array([3, 4])
    arr3 = np.array([5])
    result = skysegmentor._if_list_concatenate([arr1, arr2, arr3])
    expected = np.array([1, 2, 3, 4, 5])
    assert np.array_equal(result, expected)


def test_empty_list():
    arr = []
    result = skysegmentor._if_list_concatenate(arr)
    assert result == []  # should return the empty list unchanged


def test_single_scalar_element():
    # Although type hint says list of arrays or numbers, what if it's a scalar?
    arr = [np.array(42)]
    result = skysegmentor._if_list_concatenate(arr)
    assert np.array_equal(result, np.array([42]))


def test_mixed_float_and_int_arrays():
    arr1 = np.array([1.0, 2.0])
    arr2 = np.array([3, 4])
    result = skysegmentor._if_list_concatenate([arr1, arr2])
    expected = np.array([1.0, 2.0, 3.0, 4.0])
    assert np.allclose(result, expected)


def test_unionfinder_all_zeros():
    nside = 2
    npix = hp.nside2npix(nside)
    binmap = np.zeros(npix, dtype=int)
    groupID = skysegmentor.unionfinder(binmap)
    assert groupID.shape == binmap.shape
    assert np.all(groupID == 0), "All zeros input should yield all zeros output"


def test_unionfinder_single_active_pixel():
    nside = 2
    npix = hp.nside2npix(nside)
    binmap = np.zeros(npix, dtype=int)
    binmap[5] = 1
    groupID = skysegmentor.unionfinder(binmap)
    assert groupID[5] > 0, "Active pixel should get a positive label"
    assert np.sum(groupID == groupID[5]) == 1, "Only one pixel should have that label"


def test_unionfinder_two_disconnected_active_pixels():
    nside = 2
    npix = hp.nside2npix(nside)
    binmap = np.zeros(npix, dtype=int)
    # Select two pixels far apart (e.g., 0 and half way)
    binmap[0] = 1
    binmap[npix // 2] = 1
    groupID = skysegmentor.unionfinder(binmap)
    assert groupID[0] > 0 and groupID[npix // 2] > 0
    assert groupID[0] != groupID[npix // 2], "Disconnected pixels should have distinct labels"


def test_unionfinder_two_connected_active_pixels():
    nside = 2
    npix = hp.nside2npix(nside)
    binmap = np.zeros(npix, dtype=int)
    # Pick pixel 0 and one of its neighbors
    neighbors = hp.get_all_neighbours(nside, 0)
    neighbors = neighbors[neighbors != -1]
    if len(neighbors) == 0:
        pytest.skip("No neighbors found for pixel 0")
    binmap[0] = 1
    binmap[neighbors[0]] = 1
    groupID = skysegmentor.unionfinder(binmap)
    assert groupID[0] > 0 and groupID[neighbors[0]] > 0
    assert groupID[0] == groupID[neighbors[0]], "Connected pixels should share the same label"


def test_unionfinder_multiple_clusters():
    nside = 3
    npix = hp.nside2npix(nside)
    binmap = np.zeros(npix, dtype=int)

    # Create cluster 1: pixels 0, neighbors[0]
    cluster1 = [0]
    neighbors = hp.get_all_neighbours(nside, 0)
    neighbors = neighbors[neighbors != -1]
    if len(neighbors) > 0:
        cluster1.append(neighbors[0])
    
    # Create cluster 2: pixels npix//2 and one neighbor
    cluster2 = [npix//2]
    neighbors2 = hp.get_all_neighbours(nside, npix//2)
    neighbors2 = neighbors2[neighbors2 != -1]
    if len(neighbors2) > 0:
        cluster2.append(neighbors2[0])
    
    for pix in cluster1 + cluster2:
        binmap[pix] = 1
    
    groupID = skysegmentor.unionfinder(binmap)
    
    # Check labels for cluster 1 pixels (all same)
    label1 = groupID[cluster1[0]]
    assert all(groupID[p] == label1 for p in cluster1), "Cluster 1 pixels should share the same label"
    
    # Check labels for cluster 2 pixels (all same)
    label2 = groupID[cluster2[0]]
    assert all(groupID[p] == label2 for p in cluster2), "Cluster 2 pixels should share the same label"
    
    # Labels of clusters should differ
    assert label1 != label2, "Separate clusters should have different labels"


def test_output_integrity():
    nside = 3
    npix = hp.nside2npix(nside)
    binmap = np.random.choice([0, 1], size=npix)
    groupID = skysegmentor.unionfinder(binmap)
    
    # Shape matches input
    assert groupID.shape == binmap.shape
    
    # Output labels are integers
    assert np.issubdtype(groupID.dtype, np.integer)
    
    # All active pixels are labeled > 0
    assert np.all(groupID[binmap == 1] > 0)
    
    # All inactive pixels are zero
    assert np.all(groupID[binmap == 0] == 0)


def test_labels_are_contiguous():
    # The _shuffle_down function should ensure labels are compact (1..N)
    nside = 3
    npix = hp.nside2npix(nside)
    binmap = np.zeros(npix, dtype=int)
    
    # Randomly activate some pixels
    active_pixels = np.random.choice(npix, size=5, replace=False)
    for p in active_pixels:
        binmap[p] = 1
    
    groupID = skysegmentor.unionfinder(binmap)
    
    unique_labels = np.unique(groupID[groupID > 0])
    # Labels should be contiguous from 1 to number of unique labels
    expected_labels = np.arange(1, len(unique_labels) + 1)
    assert np.array_equal(unique_labels, expected_labels)