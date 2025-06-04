![biglogo](docs/source/_static/SkySegmentor_logo_large_github.jpg)

<p align="center">
    <a href="https://pypi.org/project/skysegmentor/">
    <img src="https://img.shields.io/pypi/v/skysegmentor.svg" alt="PyPI version">
    </a>
    <a href="https://github.com/knaidoo29/SkySegmentor">
    <img src="https://img.shields.io/badge/GitHub-repo-blue?logo=github" alt="GitHub repository">
    </a>
    <a href="https://codecov.io/github/knaidoo29/SkySegmentor" > 
    <img src="https://codecov.io/github/knaidoo29/SkySegmentor/graph/badge.svg?token=C9MXIA22X2"/> 
    </a>
    <a href="https://skysegmentor.readthedocs.io/en/latest/">
    <img src="https://readthedocs.org/projects/skysegmentor/badge/?version=latest" alt="Documentation Status">
    </a>
    <a href="https://opensource.org/licenses/MIT">
    <img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT">
    </a>
    <a href="https://github.com/psf/black">
    <img src="https://img.shields.io/badge/code%20style-black-000000.svg" alt="Code style: black">
    </a>
</p>

<!-- [![PyPI version](https://img.shields.io/pypi/v/skysegmentor.svg)](https://pypi.org/project/skysegmentor/)
[![GitHub](https://img.shields.io/badge/GitHub-repo-blue?logo=github)](https://github.com/knaidoo29/SkySegmentor)
[![Docs](https://readthedocs.org/projects/skysegmentor/badge/?version=latest)](https://skysegmentor.readthedocs.io/en/latest/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black) -->

## Introduction

**SkySegmentor** is a ``python 3`` package for dividing points or maps (in ``HEALPix``
format) on the celestial sphere into equal-sized segments. It employs a sequential 
binary space partitioning scheme -- a generalization of the *k*-d tree algorithm -- 
that supports segmentation of arbitrarily shaped sky regions. By design, all 
partitions are approximately equal in area, with discrepancies no larger than the 
``HEALPix`` pixel scale.

## Dependencies

* `numpy>=1.22,<1.27`
* `healpy>=1.15.0`

## Installation

First clone the repository

```
git clone https://github.com/knaidoo29/SkySegmentor.git
cd SkySegmentor
```

and install by running

```
pip install . [--user]
```

You should now be able to import the module:

```python
import skysegmentor
```

## Documentation

Documentation, including tutorials and explanation of API, can be found here ``https://skysegmentor.readthedocs.io/``. Alternatively a PDF version of the documentation is located in the ``docs/`` folder called ``skysegmentor.pdf``. Offline documentation can be generating by running ``make html`` in the ``docs/`` folder which will generate html documentation in the ``docs/build/html`` folder that can be accessed by opening the ``index.html`` file in a browser.

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

``SkySegmentor`` was developed as part of the Euclid angular power spectra internal covariance pipeline. If you use ``SkySegmentor`` please cite this paper, which can be found here:

**--THESE ARE PLACEHOLDERS TO BE UPDATED LATER--**
* NASA ADS:
* ArXiv: 
* BibTex:
    ```
    @ARTICLE{Naidoo2025,
            author = {{Euclid Collaboration} and {Naidoo}, K. and 
            {Ruiz-Zapatero}, J. and {Tessore}, N. and 
            {Joachimi}, B. and {Loureiro}, A. and
            others ...}
            title = "{Euclid preparation: TBD. Accurate and precise data-driven angular power spectrum covariances}
    }
    ```

and include a link to the SkySegmentor documentation page:

    https://skysegmentor.readthedocs.io/

## Support

If you have any issues with the code or want to suggest ways to improve it please open a new issue ([here](https://github.com/knaidoo29/SkySegmentor/issues))
or (if you don't have a github account) email _krishna.naidoo.11@ucl.ac.uk_.
