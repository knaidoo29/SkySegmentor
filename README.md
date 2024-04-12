# MapSplit

Author:         Krishna Naidoo                          
Version:        0.0.0                          
Homepage:       https://github.com/knaidoo29/MapSplit
Documentation:  TBA

## Introduction

MapSplit is a python 3 package for splitting binary (or weighted) HEALPix maps
into equally weighted segments. The segmentation uses a sequential binary space
partitioning scheme, a generalisation of spatial tree algorithms such as the k-d
tree. By definition all partitions are approximately equal (with errors the size
of the HEALPix pixel scale).

## Dependencies

* `numpy`
* `healpix`

## Installation

MapSplit can be installed by first cloning the repository

```
git clone https://github.com/knaidoo29/MapSplit.git
cd MapSplit
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

Similarly, if you would like to work and edit MapSplit you can clone the repository
and install an editable version:

```
git clone https://github.com/knaidoo29/MapSplit.git
cd MapSplit
pip install -e . [--user]
```

You should now be able to import the module:

```python
import mapsplit
```

## Documentation

In depth documentation and tutorials are provided TBA.

## Tutorials

The tutorials in the documentation are supplied as ipython notebooks which can be downloaded from TBA or can be run online using binder TBA.

## Citing

You can cite ``MapSplit`` using the following BibTex:

```
TBA
```

## Support

If you have any issues with the code or want to suggest ways to improve it please open a new issue ([here](https://github.com/knaidoo29/MapSplit/issues))
or (if you don't have a github account) email _krishna.naidoo.11@ucl.ac.uk_.
