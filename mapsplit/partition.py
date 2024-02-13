import numpy as np
import healpy as hp

from . import coords, rotate


def get_partition_IDs(partition):
    """Returns the total weights of each partition.

    Parameters
    ----------
    partition : array
        Partition IDs on a map. Unfilled partitions are assigned partition 0.

    Returns
    -------
    partition_IDs : array
        Unique partition IDs, including zero.
    """
    partition_IDs = np.unique(partition)
    return partition_IDs


def total_partition_weights(partition, weights):
    """Returns the total weights of each partition.

    Parameters
    ----------
    partition : array
        Partition IDs on a map. Unfilled partitions are assigned partition 0.
    weights : array
        A weight assigned to each element of the partition array.

    Returns
    -------
    partition_IDs : array
        Unique partition IDs, including zero.
    partition_weights : array
        The total weight for each partition.
    """
    partition_IDs = get_partition_IDs(partition)
    Npartitions = np.max(partition_IDs)+1
    #partitions_weights = magpie.pixels.bin_pix(partition, Npartition, weights=weights)
    partitions_weights = np.zeros(Npartition)
    np.add.at(partition_weights, partition_IDs, weights)
    return partition_IDs, partition_weights


def remove_val4array(array, val):
    """Removes a given value from an array.

    Parameters
    ----------
    array : array
        Data vector.
    val : float
        Value to be removed from an array.
    """
    return array[array != val]


def fill_map(pixID, nside, val=1.):
    """Fill a Healpix map with a given value val at given pixel locations.

    Parameters
    ----------
    pixID : int array
        Pixel index.
    nside : int
        HEalpix map nside.
    val : float, optional
        Value to fill certain pixels with.
    """
    tmap = np.zeros(hp.nside2npix(nside))
    tmap[pixID] = val
    return tmap


def find_boundary_pix(nside, bnmap):
    """Locate pixels on the boundaries of a binary map.

    Parameters
    ----------
    nside : int
        Healpix map nside.
    bnmap : int array
        Binary healpix map.

    Returns
    -------
    pixboundID : array
        Boundary pixel IDs of the given healpix binary map.
    """
    pixID = np.nonzero(bnmap)[0]
    pixneiID = hp.get_all_neighbours(nside, pixID)
    isbound = [any(bnmap[remove_val4array(pixneiID.T[i], -1)]==0) for i in range(0, len(pixneiID.T))]
    cond = np.where(np.array(isbound))[0]
    pixboundID = pixID[cond]
    return pixboundID


def get_most_dist_points(nside, bnmap):
    """Returns the most distant points on a binary map. Note there is an implicit
    assumption that the map does not span regions larger than a hemisphere.

    Parameters
    ----------
    nside : int
        HEalpix map nside.
    bnmap : int array
        Binary healpix map.

    Returns
    -------
    p1, t1, p2, t2 : float
        Angular coordinates (phi, theta) for the most distant points (1 and 2) on
        the binary map.
    """
    pixboundID = find_boundary_pix(nside, bnmap)
    the_bound, phi_bound = hp.pix2ang(nside, pixboundID)

    pp1, pp2 = np.meshgrid(phi_bound, phi_bound, indexing='ij')
    tt1, tt2 = np.meshgrid(the_bound, the_bound, indexing='ij')
    pp1, pp2, tt1, tt2 = pp1.flatten(), pp2.flatten(), tt1.flatten(), tt2.flatten()

    dist = coords.distusphere(pp1, tt1, pp2, tt2)

    ind = np.argmax(dist)
    p1, p2 = pp1[ind], pp2[ind]
    t1, t2 = tt1[ind], tt2[ind]
    return p1, t1, p2, t2


def weight_dif(phi_split, phi, weights, balance=1):
    """Compute the difference between the weights on either side of phi_split.

    Parameters
    ----------
    phi_split : float
        Longitude split.
    phi : array
        Longitude coordinates.
    weights : array
        Weights corresponding to each longitude coordinates.
    balance : float, optional
        A multiplication factor assigned to weights below phi_split.
    """
    cond = np.where(phi <= phi_split)[0]
    weights1 = balance*np.sum(weights[cond])
    cond = np.where(phi > phi_split)[0]
    weights2 = np.sum(weights[cond])
    return abs(weights1-weights2)


def find_dphi(phi, weights, balance=1):
    """Determines the splitting longitude required for partitioning, either with
    1-to-1 weights on either side or unbalanced weighting if balance != 1.

    Parameters
    ----------
    hi : array
        Longitude coordinates.
    weights : array
        Weights corresponding to each longitude coordinates.

    Returns
    -------
    dphi : float
        Splitting longitude.
    """
    dphis = np.linspace(phi.min(), phi.max(), 100)
    _dphi = dphis[1]-dphis[0]
    weights_dif = np.array([weight_dif(dphi, phi, weights, balance=balance) for dphi in dphis])

    ind = np.argmin(weights_dif)
    dphi = dphis[ind]

    dphis = np.linspace(dphi-2*_dphi, dphi+2*_dphi, 100)
    weights_dif = np.array([weight_dif(dphi, phi, weights, balance=balance) for dphi in dphis])

    ind = np.argmin(weights_dif)
    dphi = dphis[ind]
    return dphi


def split_into_2(weightmap, balance=1, partitionmap=None, partition=None):
    """Splits a map with weights into 2 equal (unequal in balance != 1).

    Parameters
    ----------
    weightmap : array
        Healpix weight map.
    balance : float, optional
        Balance of the weights for the partitioning.
    partitionmap : int array, optional
        Partitioned map IDs.
    partition : int, optional
        A singular partition to be partitioned in two pieces.

    Returns
    -------
    partitionmap : int array
        Partitioned map IDs.
    """
    npix = len(weightmap)
    nside = hp.npix2nside(npix)
    pixID = np.nonzero(weightmap)[0]
    bnmap = np.zeros(npix)
    bnmap[pixID] = 1

    if partitionmap is None:
        partitionmap = np.copy(bnmap)
        maxpartition = 1
        partition = 1
    else:
        maxpartition = int(np.max(partitionmap))

    _pixID = np.where(partitionmap == partition)[0]
    _bnmap = np.zeros(npix)
    _bnmap[_pixID] = 1

    _the, _phi = hp.pix2ang(nside, _pixID)
    _weights = weightmap[_pixID]

    p1, t1, p2, t2 = get_most_dist_points(nside, _bnmap)
    a1, a2, a3 = rotate.rotate2plane([p1, t1], [p2, t2])

    _phi, _the = rotate.forward_rotate(_phi, _the, a1, a2, a3)

    _dphi = find_dphi(_phi, _weights, balance=balance)

    _cond = np.where(_phi > _dphi)[0]
    partitionmap[_pixID[_cond]] = maxpartition + 1

    return partitionmap


def split_into_N(weightmap, Npartitions):
    """Splits a map with weights into equal Npartition sides.

    Parameters
    ----------
    weightmap : array
        Healpix weight map.
    Npartitions : int
        Number of partitioned regions

    Returns
    -------
    partitionmap : int array
        Partitioned map IDs.
    """
    # The number of partitions currently assigned for each partition ID.
    part_Npart = np.zeros(Npartitions)
    part_Npart[0] = Npartitions
    maxpartition = 1

    partitionmap = np.zeros(len(weightmap))
    pixID = np.nonzero(weightmap)
    partitionmap[pixID] = 1.

    while any(part_Npart == 0):
        for i in range(0, len(part_Npart)):
            partition = i+1
            if part_Npart[i] > 1:
                wei1 = int(np.floor(part_Npart[i]/2.))
                wei2 = part_Npart[i] - wei1
                part_Npart[i] = wei1
                part_Npart[maxpartition] = wei2
                maxpartition += 1

                balance = wei2/wei1

                partitionmap = split_into_2(weightmap, balance, partitionmap=partitionmap,
                    partition=partition)

    return partitionmap
