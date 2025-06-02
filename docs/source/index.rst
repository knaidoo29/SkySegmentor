
.. image:: _static/SkySegmentor_logo_large_white.jpg
   :align: center
   :class: only-light

.. image:: _static/SkySegmentor_logo_large_black.jpg
   :align: center
   :class: only-dark
  
============
SkySegmentor
============

.. list-table::
    :widths: 20 100
    :header-rows: 0
    :stub-columns: 1

    * - Version
      - 0.0.6
    * - Repository
      - https://github.com/knaidoo29/SkySegmentor
    * - Documentation
      - https://skysegmentor.readthedocs.io/

Contents
========

* `Introduction`_
* `Dependencies`_
* `Installation`_
* `Tutorials`_
* `API`_
* `Contributors`_
* `Support`_
* `Version History`_

Introduction
============

**SkySegmentor** is a ``python 3`` package for dividing points or maps (in ``HEALPix``
format) on the celestial sphere into equal-sized segments. It employs a sequential 
binary space partitioning scheme -- a generalization of the *k*-d tree algorithm -- 
that supports segmentation of arbitrarily shaped sky regions. By design, all 
partitions are approximately equal in area, with discrepancies no larger than the 
``HEALPix`` pixel scale.

Dependencies
============

* `numpy <http://www.numpy.org/>`_
* `healpy <https://healpy.readthedocs.io/en/latest/>`_

Installation
============

SkySegmentor can be installed by first cloning the repository

.. code-block:: bash

    git clone https://github.com/knaidoo29/SkySegmentor.git
    cd SkySegmentor

and install by either running

.. code-block:: bash

    pip install . [--user]

or

.. code-block:: bash

    python setup.py build
    python setup.py install

You should now be able to import the module:

.. code-block:: python

    import skysegmentor

Tutorials
=========

Basic Usage
-----------

Segmenting Healpix Maps
^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    import healpy
    import skysegmentor

    # Healpix mask, where zeros are regions outside of the mask and ones inside the
    # mask. You can also input a weighted map, where instead of 1s you give weights.
    mask = # define mask values

    Npartitions = 100 # Number of partitions
    partitionmap = skysegmentor.segmentmapN(mask, Npartitions)

Segmenting Points on the Sphere
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    import skysegmentor

    # Define points on the sphere to be segmented.
    phi = # longitude defined in radians from [0, 2*pi]
    the = # latitude defined in radians from [0, pi], where 0 = North Pole.

    Npartitions = 100 # Number of partitions
    partitionIDs = skysegmentor.segmentpointsN(phi, the, Npartitions)

if using RA and Dec in degrees you can convert to phi and the using

.. code-block:: python

    phi = np.deg2rad(ra)
    the = np.deg2rad(90. - dec)

if not all points are equal, you can specify a weight

.. code-block:: python

    weights = # define point weights
    partitionIDs = skysegmentor.segmentpointsN(phi, the, Npartitions, weights=weights)

.. toctree::
  :maxdepth: 2

  tutorials

API
===

.. toctree::
  :maxdepth: 2

  api


Contributors
============

If you use ``SkySegmentor`` in a publication please cite::

    TBA

and include a link to the SkySegmentor main page::

    https://github.com/knaidoo29/SkySegmentor

Support
=======

If you have any issues with the code or want to suggest ways to improve it please
open a new issue (`here <https://github.com/knaidoo29/SkySegmentor/issues>`_) or
(if you don't have a github account) email krishna.naidoo.11@ucl.ac.uk.

Version History
===============

* **Version 0.0.0**:
