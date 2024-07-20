============
SkySegmentor
============

.. list-table::
    :widths: 10 60
    :header-rows: 0
    :stub-columns: 1

    * - Author
      - Krishna Naidoo
    * - Version
      - 0.0.0
    * - Repository
      - https://github.com/knaidoo29/SkySegmentor
    * - Documentation
      - TBA

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

**SkySegmentor** is a ``python 3`` package for splitting binary (or weighted) ``HEALPix``
maps or points on the sphere into equally weighted segments. The segmentation uses
a sequential binary space partitioning scheme, a generalisation of the k-d tree
algorithm. By definition all partitions are approximately equal (with errors the
size of the ``HEALPix`` pixel scale).

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
