
.. image:: _static/SkySegmentor_logo_large_white.jpg
   :align: center
   :class: only-light

.. image:: _static/SkySegmentor_logo_large_black.jpg
   :align: center
   :class: only-dark

.. image:: https://img.shields.io/pypi/v/skysegmentor.svg
   :target: https://pypi.org/project/skysegmentor/
   :alt: PyPI version

.. image:: https://img.shields.io/badge/GitHub-repo-blue?logo=github
   :target: https://github.com/knaidoo29/SkySegmentor
   :alt: GitHub repository

.. image:: https://readthedocs.org/projects/skysegmentor/badge/?version=latest
   :target: https://skysegmentor.readthedocs.io/en/latest/
   :alt: Documentation Status

.. image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: https://opensource.org/licenses/MIT
   :alt: License: MIT

Contents
========

* `Introduction`_
* `Dependencies`_
* `Installation`_
* `Tutorials`_
* `API`_
* `Citing`_
* `Support`_

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

* `numpy <http://www.numpy.org/>`_ -- versions: ``>=1.22,<1.27``
* `healpy <https://healpy.readthedocs.io/en/latest/>`_ -- versions: ``>=1.15.0``

Installation
============

First cloning the repository

.. code-block:: bash

    git clone https://github.com/knaidoo29/SkySegmentor.git
    cd SkySegmentor

and install by running

.. code-block:: bash

    pip install . [--user]


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


Citing
======

.. warning::

  These are placeholders to be replaced upon publication on ArXiv.

If you use ``SkySegmentor`` in a publication please cite:

* NASA ADS:
* ArXiv:

BibTex::

    @ARTICLE{Naidoo2025,
        author = {{Euclid Collaboration} and {Naidoo}, K. and {Ruiz-Zapatero}, J. 
        and {Tessore}, N. and {Joachimi}, B. and {Loureiro}, A. and others ...},
        title = "{Euclid preparation: TBD. Accurate and precise data-driven angular power spectrum covariances}"
     }

and include a link to the SkySegmentor documentation page: https://skysegmentor.readthedocs.io/

Support
=======

If you have any issues with the code or want to suggest ways to improve it please
open a new issue (`here <https://github.com/knaidoo29/SkySegmentor/issues>`_) or
(if you don't have a github account) email krishna.naidoo.11@ucl.ac.uk.
