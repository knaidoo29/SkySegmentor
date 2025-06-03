Basic Usage
===========

Segmenting a Healpix Maps
-------------------------

.. code-block:: python

    import healpy
    import skysegmentor

    # Healpix mask, where zeros are regions outside of the mask and ones inside the
    # mask. You can also input a weighted map, where instead of 1s you give weights.
    mask = # define mask values

    Npartitions = 100 # Number of partitions
    partitionmap = skysegmentor.segmentmapN(mask, Npartitions)

.. image:: _static/partitionmap1.png
    :align: center

Segmenting Points on the Sphere
-------------------------------

.. code-block:: python

    import skysegmentor

    # Define points on the sphere to be segmented.
    phi = # longitude defined in radians from [0, 2*pi]
    the = # latitude defined in radians from [0, pi], where 0 = North Pole.

    Npartitions = 100 # Number of partitions
    partitionIDs = skysegmentor.segmentpointsN(phi, the, Npartitions)

if using RA and Dec in degrees you can convert to logitudinal and latitude coordinates `phi` and `the` using

.. code-block:: python

    phi = np.deg2rad(ra)
    the = np.deg2rad(90. - dec)

.. image:: _static/partitionpoints.png
    :align: center
