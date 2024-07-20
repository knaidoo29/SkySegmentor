# SkySegmentor

|               |                                             |
|---------------|---------------------------------------------|
| Author        | Krishna Naidoo                              |               
| Version       | 0.0.0                                       |
| Repository    | https://github.com/knaidoo29/SkySegmentor   |
| Documentation | TBA                                         |

## Introduction

SkySegmentor is a python 3 package for splitting binary (or weighted) HEALPix maps
or points on the sphere into equally weighted segments. The segmentation uses a
sequential binary space partitioning scheme, a generalisation of the k-d tree
algorithm. By definition all partitions are approximately equal (with errors the
size of the HEALPix pixel scale).

## Dependencies

* `numpy`
* `healpy`

## Installation

SkySegmentor can be installed by first cloning the repository

```
git clone https://github.com/knaidoo29/SkySegmentor.git
cd SkySegmentor
```

and install by either running

```
pip install . [--user]
```

or

```
python setup.py build
python setup.py install
```

You should now be able to import the module:

```python
import skysegmentor
```

## Documentation

Documentation, including tutorials and explanation of API, can be found here TBA.
Alternatively a PDF version of the documentation is located in the ``docs/`` folder
called ``skysegmentor.pdf``. Offline documentation can be generating by running
``make html`` in the ``docs/`` folder which will generate html documentation in
the ``docs/build/html`` folder that can be accessed by opening the ``index.html``
file in a browser.

## Tutorial

### Basic Usage

#### Segmenting Healpix Maps

```python
import healpy
import skysegmentor

# Healpix mask, where zeros are regions outside of the mask and ones inside the
# mask. You can also input a weighted map, where instead of 1s you give weights.
mask = # define mask values

Npartitions = 100 # Number of partitions
partitionmap = skysegmentor.segmentmapN(mask, Npartitions)
```

#### Segmenting Points on the Sphere

```python
import skysegmentor

# Define points on the sphere to be segmented.
phi = # longitude defined in radians from [0, 2*pi]
the = # latitude defined in radians from [0, pi], where 0 = North Pole.

Npartitions = 100 # Number of partitions
partitionIDs = skysegmentor.segmentpointsN(phi, the, Npartitions)
```

if using RA and Dec in degrees you can convert to phi and the using

```python
phi = np.deg2rad(ra)
the = np.deg2rad(90. - dec)
```

if not all points are equal, you can specify a weight

```python
weights = # define point weights
partitionIDs = skysegmentor.segmentpointsN(phi, the, Npartitions, weights=weights)
```

## Citing

You can cite ``SkySegmentor`` using the following BibTex:

```
TBA
```

## Support

If you have any issues with the code or want to suggest ways to improve it please open a new issue ([here](https://github.com/knaidoo29/SkySegmentor/issues))
or (if you don't have a github account) email _krishna.naidoo.11@ucl.ac.uk_.
