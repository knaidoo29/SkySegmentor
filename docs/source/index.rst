
.. image:: _static/SkySegmentor_logo_large_white.jpg
   :align: center
   :class: only-light

.. image:: _static/SkySegmentor_logo_large_black.jpg
   :align: center
   :class: only-dark

.. raw:: html
   
   <p align="center">
      <a href="https://github.com/knaidoo29/SkySegmentor/actions/workflows/python-tests.yml">
      <img src="https://github.com/knaidoo29/SkySegmentor/actions/workflows/python-tests.yml/badge.svg" alt="Python Tests">
      </a>
      <a href="https://codecov.io/github/knaidoo29/SkySegmentor" > 
      <img src="https://codecov.io/github/knaidoo29/SkySegmentor/graph/badge.svg?token=C9MXIA22X2"/> 
      </a>
      <a href="https://img.shields.io/badge/Python-3.9%20|%203.10%20|%203.11%20|%203.12-blue">
      <img src="https://img.shields.io/badge/Python-3.9%20|%203.10%20|%203.11%20|%203.12-blue" alt="Python Version Support">
      </a>
      <a href="https://img.shields.io/github/v/release/knaidoo29/skysegmentor">
      <img src="https://img.shields.io/github/v/release/knaidoo29/skysegmentor" alt="Version">
      </a>
      <a href="https://pypi.org/project/skysegmentor/">
      <img src="https://img.shields.io/pypi/v/skysegmentor.svg" alt="PyPI version">
      </a>
      <a href="https://skysegmentor.readthedocs.io/en/latest/">
      <img src="https://readthedocs.org/projects/skysegmentor/badge/?version=latest" alt="Documentation Status">
      </a>
      <a href="https://github.com/knaidoo29/SkySegmentor">
      <img src="https://img.shields.io/badge/GitHub-repo-blue?logo=github" alt="GitHub repository">
      </a>
      <a href="https://opensource.org/licenses/MIT">
      <img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT">
      </a>
      <a href="https://github.com/psf/black">
      <img src="https://img.shields.io/badge/code%20style-black-000000.svg" alt="Code style: black">
      </a>
   </p>

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

**SkySegmentor** is a ``python`` package for dividing points or maps (in ``HEALPix``
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

Pip Installation
----------------

.. code-block:: bash

   pip install skysegmentor

From the source
---------------

Clone the repository

.. code-block:: bash

    git clone https://github.com/knaidoo29/SkySegmentor.git
    cd SkySegmentor

and install by running

.. code-block:: bash

    pip install . [--user]


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
