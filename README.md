# AbundanceMatching

[![PyPI:AbundanceMatching](https://img.shields.io/pypi/v/AbundanceMatching.svg)](https://pypi.python.org/pypi/AbundanceMatching)
[![ascl:1604.006](https://img.shields.io/badge/ascl-1604.006-blue.svg?colorB=262255)](https://ascl.net/1604.006)

A python module that implements subhalo abundance matching. 
This module can interpolate and extrapolate abundance functions (e.g., stellar mass function, halo mass function) 
and provides Peter Behroozi's fiducial deconvolution implementation ([Behroozi et al. 2010](https://ui.adsabs.harvard.edu/abs/2010ApJ...717..379B/abstract)).

## Installation

```bash
pip install abundancematching
```

## Example

Here's an example to do abundance matching with this code.

```python
"""
Assume you have a numpy structured array `halos`,
which contains a list of halos, with labels of the quantity names.
Assume you also have a luminosity function table `lf`,
whose first column st column is the quantity to match (e.g. magnitude),
and the second column is the abundance (per Mpc^3 per Mag).
"""

import matplotlib.pyplot as plt
from AbundanceMatching import AbundanceFunction, LF_SCATTER_MULT, calc_number_densities

af = AbundanceFunction(lf[:,0], lf[:,1], (-27, -5))

# check the abundance function
plt.semilogy(lf[:,0], lf[:,1])
x = np.linspace(-27, -5, 101)
plt.semilogy(x, af(x))

# deconvolution and check results (it's a good idea to always check this)
scatter = 0.2
remainder = af.deconvolute(scatter*LF_SCATTER_MULT, 20)
x, nd = af.get_number_density_table()
plt.plot(x, remainder/nd);

# get number densities of the halo catalog
nd_halos = calc_number_densities(halos['vpeak'], box_size)

# do abundance matching with no scatter
catalog = af.match(nd_halos)

#do abundance matching with some scatter
catalog_sc = af.match(nd_halos, scatter*LF_SCATTER_MULT)

#if you want to do multiple (100) realizations:
catalog_deconv = af.match(nd_halos, scatter*LF_SCATTER_MULT, False)
for __ in range(100):
    catalog_this = add_scatter(catalog_deconv, scatter*LF_SCATTER_MULT)
    catalog_this = rematch(catalog_this, catalog, af._x_flipped)
    # do something with catalog_this
```
