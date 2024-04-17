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


def find_map_barycenter(bnmap, wmap=None):
    """Determines the barycenter of center of mass direction of the input binary map.

    Parameters
    ----------
    bnmap : array
        binary map.
    wmap : array, optional
        The weights.

    Returns
    -------
    phic, thec : float
        The center
    """
    if wmap is None:
        wmap = np.ones(len(bnmap))
    nside = hp.npix2nside(len(bnmap))
    pixID = np.where(bnmap != 0.)[0]
    the, phi = hp.pix2ang(nside, pixID)
    wei = wmap[pixID]
    x, y, z = coords.sphere2cart(np.ones(len(phi)), phi, the, center=[0., 0., 0.])
    xc = np.sum(x*wei)/np.sum(wei)
    yc = np.sum(y*wei)/np.sum(wei)
    zc = np.sum(z*wei)/np.sum(wei)
    _, phic, thec = coords.cart2sphere(xc, yc, zc)
    phir, ther = rotate.rotate_usphere(phi, the, [-phic, -thec, 0.])
    themax = np.max(ther)
    return phic, thec, themax


def find_points_barycenter(phi, the, weights=None):
    """Determines the barycenter of center of mass direction of the input point dataset.

    Parameters
    ----------
    phi, the : array
        Angular coordinates.
    weights : array, optional
        Weights for points.

    Returns
    -------
    phic, thec : float
        The center
    """
    if weights is None:
        weights = np.ones(len(phi))
    x, y, z = coords.sphere2cart(np.ones(len(phi)), phi, the, center=[0., 0., 0.])
    xc = np.sum(x*weights)/np.sum(weights)
    yc = np.sum(y*weights)/np.sum(weights)
    zc = np.sum(z*weights)/np.sum(weights)
    _, phic, thec = coords.cart2sphere(xc, yc, zc)
    phir, ther = rotate.rotate_usphere(phi, the, [-phic, -thec, 0.])
    themax = np.max(ther)
    return phic, thec, themax


def get_map_border(bnmap, wmap=None, res=[200, 100]):
    """Determines the outer border of binary map region.

    Parameters
    ----------
    bnmap : array
        binary map.
    wmap : array, optional
        The weights.
    res : list, optional
        Resolution of spherical cap grid where [phiresolution, thetaresolution]
        to find region border.

    Returns
    -------
    phi_border, the_border : float
        Approximate border region.
    """
    nside = hp.npix2nside(len(bnmap))
    phic, thec, themax = find_map_barycenter(bnmap, wmap=wmap)

    psize = res[0]
    tsize = res[1]

    pedges = np.linspace(0., 2*np.pi, psize + 1)
    pmid = 0.5*(pedges[1:] + pedges[:-1])
    tedges = np.linspace(0., np.max(themax)*1.05, tsize + 1)
    tmid = 0.5*(tedges[1:] + tedges[:-1])

    pcap, tcap = np.meshgrid(pmid, tmid, indexing='ij')
    pshape = np.shape(pcap)

    pcap_rot, tcap_rot = rotate.rotate_usphere(pcap.flatten(), tcap.flatten(), [0., thec, phic])
    pixID = hp.ang2pix(nside, tcap_rot, pcap_rot)
    wcap_rot = bnmap[pixID]

    pcap_rot = pcap_rot.reshape(pshape)
    tcap_rot = tcap_rot.reshape(pshape)
    wcap_rot = wcap_rot.reshape(pshape)

    phi_border, the_border = [], []

    tind = np.arange(len(tmid))

    for i in range(len(wcap_rot)):
        if len(tind[wcap_rot[i] != 0.]) > 0:
            ind = np.max(tind[wcap_rot[i] != 0.])
            phi_border.append(pcap_rot[i,ind])
            the_border.append(tcap_rot[i,ind])

    phi_border = np.array(phi_border)
    the_border = np.array(the_border)

    return phi_border, the_border


def get_points_border(phi, the, weights=None, res=100):
    """Determines the outer border of binary map region.

    Parameters
    ----------
    phi, the : array
        Angular coordinates.
    weights : array, optional
        Weights for points.
    res : int, optional
        Resolution of spherical cap grid for phiresolution to find region border.

    Returns
    -------
    phi_border, the_border : float
        Approximate border region.
    """
    phic, thec, themax = find_points_barycenter(phi, the, weights=weights)

    pedges = np.linspace(0., 2*np.pi, res + 1)

    phi_rot, the_rot = rotate.rotate_usphere(phi, the, [-phic, -thec, 0.])

    phi_border, the_border = [], []

    for i in range(0, len(pedges)-1):
        cond = np.where((phi_rot >= pedges[i]) & (phi_rot <= pedges[i+1]))[0]
        if len(cond) > 0:
            ind = cond[np.argmax(the_rot[cond])]
            phi_border.append(phi[ind])
            the_border.append(the[ind])

    phi_border = np.array(phi_border)
    the_border = np.array(the_border)

    return phi_border, the_border


def get_map_most_dist_points(nside, bnmap, wmap=None, res=[100, 50]):
    """Returns the most distant points on a binary map.

    Parameters
    ----------
    nside : int
        HEalpix map nside.
    bnmap : int array
        Binary healpix map.
    wmap : array, optional
        The weights.
    res : list, optional
        Resolution of spherical cap grid where [phiresolution, thetaresolution]
        to find region border.

    Returns
    -------
    p1, t1, p2, t2 : float
        Angular coordinates (phi, theta) for the most distant points (1 and 2) on
        the binary map.
    """

    phi_border, the_border = get_map_border(bnmap, wmap=wmap, res=res)

    pp1, pp2 = np.meshgrid(phi_border, phi_border, indexing='ij')
    tt1, tt2 = np.meshgrid(the_border, the_border, indexing='ij')
    pp1, pp2, tt1, tt2 = pp1.flatten(), pp2.flatten(), tt1.flatten(), tt2.flatten()

    dist = coords.distusphere(pp1, tt1, pp2, tt2)

    ind = np.argmax(dist)
    p1, p2 = pp1[ind], pp2[ind]
    t1, t2 = tt1[ind], tt2[ind]

    return p1, t1, p2, t2


def get_points_most_dist_points(phi, the, weights=None, res=100):
    """Returns the most distant points from a set of points.

    Parameters
    ----------
    phi, the : array
        Angular coordinates.
    weights : array, optional
        Weights for points.
    res : int, optional
        Resolution of spherical cap grid for phiresolution to find region border.

    Returns
    -------
    p1, t1, p2, t2 : float
        Angular coordinates (phi, theta) for the most distant points (1 and 2) on
        the binary map.
    """

    phi_border, the_border = get_points_border(phi, the, weights=weights, res=res)

    pp1, pp2 = np.meshgrid(phi_border, phi_border, indexing='ij')
    tt1, tt2 = np.meshgrid(the_border, the_border, indexing='ij')
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


def splitmap2(weightmap, balance=1, partitionmap=None, partition=None, res=[100, 50]):
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
    res : list, optional
        Resolution of spherical cap grid where [phiresolution, thetaresolution]
        to find region border.

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

    if np.sum(_bnmap) != len(_bnmap):

        p1, t1, p2, t2 = get_map_most_dist_points(nside, _bnmap, wmap=weightmap, res=res)

        a1, a2, a3 = rotate.rotate2plane([p1, t1], [p2, t2])

        _phi, _the = rotate.forward_rotate(_phi, _the, a1, a2, a3)

    _dphi = find_dphi(_phi, _weights, balance=balance)

    _cond = np.where(_phi > _dphi)[0]
    partitionmap[_pixID[_cond]] = maxpartition + 1

    return partitionmap


def splitpoints2(phi, the, weights=None, balance=1, partitionID=None,
    partition=None, res=100):
    """Splits a map with weights into 2 equal (unequal in balance != 1).

    Parameters
    ----------
    phi, the : array
        Angular positions.
    weights : array, optional
        Angular position weights.
    balance : float, optional
        Balance of the weights for the partitioning.
    partitionID : int array, optional
        Partitioned map IDs.
    partition : int, optional
        A singular partition to be partitioned in two pieces.
    res : float, optional
        Resolution of spherical cap phiresolution to find region border.

    Returns
    -------
    partitionID : int array
        Partitioned map IDs.
    """

    if weights is None:
        weights = np.ones(len(phi))

    if partitionID is None:
        partitionID = np.ones(len(phi))
        maxpartition = 1
        partition = 1
    else:
        maxpartition = int(np.max(partitionID))

    _pixID = np.where(partitionID == partition)[0]
    _bnmap = np.ones(len(_pixID))

    _phi, _the = phi[_pixID], the[_pixID]
    _weights = weights[_pixID]

    p1, t1, p2, t2 = get_points_most_dist_points(_phi, _the, weights=_weights, res=res)

    a1, a2, a3 = rotate.rotate2plane([p1, t1], [p2, t2])

    _phi, _the = rotate.forward_rotate(_phi, _the, a1, a2, a3)

    _dphi = find_dphi(_phi, _weights, balance=balance)

    _cond = np.where(_phi > _dphi)[0]
    partitionID[_pixID[_cond]] = maxpartition + 1

    return partitionID


def splitmapN(weightmap, Npartitions, res=[100, 50]):
    """Splits a map with weights into equal Npartition sides.

    Parameters
    ----------
    weightmap : array
        Healpix weight map.
    Npartitions : int
        Number of partitioned regions
    res : list, optional
        Resolution of spherical cap grid where [phiresolution, thetaresolution]
        to find region border.

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

                partitionmap = splitmap2(weightmap, balance=balance, partitionmap=partitionmap,
                    partition=partition, res=res)

    return partitionmap


def splitpointsN(phi, the, Npartitions, weights=None, res=100):
    """Splits a map with weights into equal Npartition sides.

    Parameters
    ----------
    weightmap : array
        Healpix weight map.
    Npartitions : int
        Number of partitioned regions
    res : list, optional
        Resolution of spherical cap grid where [phiresolution, thetaresolution]
        to find region border.

    Returns
    -------
    partitionmap : int array
        Partitioned map IDs.
    """
    # The number of partitions currently assigned for each partition ID.
    part_Npart = np.zeros(Npartitions)
    part_Npart[0] = Npartitions
    maxpartition = 1

    if weights is None:
        weights = np.ones(len(phi))

    partitionID = np.ones(len(weights))

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

                partitionID = splitpoints2(phi, the, weights=weights, balance=balance, partitionID=partitionID,
                    partition=partition, res=res)

    return partitionID
