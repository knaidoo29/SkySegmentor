import numpy as np
import healpy as hp
from typing import List, Tuple, Union


def _cascade(labels: np.ndarray, indexin: int) -> int:
    """Cascades down a linked list of label reassignments.

    Parameters
    ----------
    labels : int array
        Label array.
    indexin : int
        Label index in.

    Returns
    -------
    indexout : out
        Cascaded label index out.
    """
    indexout = indexin
    while labels[indexout - 1] != indexout:
        indexout = labels[indexout-1]
    return indexout


def _cascade_all(labels: np.ndarray) -> np.ndarray:
    """Cascade label index for an array.

    Parameters
    ----------
    labels : int array
        Label array.

    Returns
    -------
    labelsout : int array
        Cascade all label index in an array.
    """
    labelsout = np.copy(labels)
    for i in range (0, len(labels)):
        labelsout[i] = _cascade(labels, labels[i])
    return labelsout


def _unionise(ind1: int, ind2: int, labels: np.ndarray) -> Tuple[int, int, int]:
    """Finds the union of two label indexes.

    Parameters
    ----------
    ind1, ind2 : int
        Index 1 and 2.
    labels : int array
        Label array.

    Returns
    -------
    ind1out, ind2out : int
        Outputted index 1 and 2.
    indout : int
        Index out.
    """
    ind1out = _cascade(labels, ind1)
    ind2out = _cascade(labels, ind2)
    if ind1out <= ind2out:
        indout = ind1out
    else:
        indout = ind2out
    return ind1out, ind2out, indout


def _shuffle_down(labels: np.ndarray) -> np.ndarray:
    """Shuffle label index for an array.

    Parameters
    ----------
    labels : int array
        Label array.

    Returns
    -------
    labelsout : int array
        Cascade all label index in an array.
    """
    labelsout = np.copy(labels)
    nlabels = np.zeros(len(labels))
    for i in range(0, len(labels)):
        nlabels[labels[i]-1] = nlabels[labels[i]-1] + 1
    maplabel = np.zeros(len(labels)).astype('int')
    j = 0
    for i in range(0, len(labels)):
        if nlabels[i] > 0:
            j = j + 1
            maplabel[i] = j
        else:
            maplabel[i] = 0
    for i in range(0, len(labels)):
        labelsout[i] = maplabel[labels[i]-1]
    return labelsout


def _if_list_concatenate(arr: List[Union[float, int]]) -> np.ndarray:
    """Concatenates list only if length is > 1."""
    if len(arr) > 1:
        arr = np.concatenate(arr)
    return arr


def unionfinder(binmap: np.ndarray) -> np.ndarray:
    """Group or label assignment on a healpix grid using the HoshenKopelman algorithm.

    Parameters
    ----------
    binmap : array
        Binary healpix map.

    Returns
    -------
    groupID : array
        Labelled healpix map.
    """
    """Group or label assignment on a healpix grid using the HoshenKopelman algorithm.

    Parameters
    ----------
    binmap : array
        Binary healpix map.

    Returns
    -------
    groupID : array
        Labelled healpix map.
    """
    nside = hp.npix2nside(len(binmap))
    groupID = np.zeros(hp.nside2npix(nside))
    groupID = groupID.astype('int')
    cond = np.where(binmap == 1)[0]
    groupID[cond] = -1
    groupID_equals = []
    currentID = 0
    for i in range(0, len(binmap)):
        if binmap[i] == 1 and groupID[i] == -1:
            currentID += 1
            groupID[i] = currentID
            groupID_equals.append([])
        if groupID[i] != 0:
            neighs = hp.get_all_neighbours(nside, i)
            neighs = neighs[neighs != -1]
            neighs = neighs[groupID[neighs] != 0]
            if len(neighs) > 0:
                groupID[neighs[groupID[neighs] == -1]] = groupID[i]
                neighs = neighs[groupID[neighs] != groupID[i]]
            if len(neighs) > 0:
                groupID_equals[groupID[i]-1].append(np.unique(groupID[neighs]))
    groupID_equals2 = [np.unique(_if_list_concatenate(_groupID)) for _groupID in groupID_equals]
    groupID_ind = np.arange(len(groupID_equals)) + 1
    groupID_pair1 = []
    groupID_pair2 = []
    for i in range(0, len(groupID_equals2)):
        for j in range(0, len(groupID_equals2[i])):
            groupID_pair1.append(groupID_ind[i])
            groupID_pair2.append(groupID_equals2[i][j])
    groupID_pair1 = np.array(groupID_pair1)
    groupID_pair2 = np.array(groupID_pair2)
    cond = np.where(groupID_pair1 > groupID_pair2)[0]
    temp = groupID_pair2[cond]
    groupID_pair2[cond] = groupID_pair1[cond]
    groupID_pair1[cond] = temp
    for i in range(0, len(groupID_pair1)):
        ind1out, ind2out, indout = _unionise(groupID_pair1[i], groupID_pair2[i], groupID_ind)
        groupID_ind[ind1out-1] = indout
        groupID_ind[ind2out-1] = indout
    groupID_ind = _cascade_all(groupID_ind)
    groupID_ind = _shuffle_down(groupID_ind)
    cond = np.where(groupID != 0)[0]
    groupID[cond] = groupID_ind[groupID[cond]-1]
    return groupID
