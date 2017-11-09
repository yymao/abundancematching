# AbundanceMatching

A python module to create (interpolate and extrapolate) abundance functions and also provide fiducial deconvolution (with Peter Behroozi's implementation) and abundance matching.


## Installation

    pip install AbundanceMatching

## Example

Here's an example to do abundance matching with this code.

    #Assume you have a numpy structured array `halos`, 
    #which contains a list of halos, with labels of the quantity names.
    #Assume you also have a luminosity function table `lf`, 
    #whose first column st column is the quantity to match (e.g. magnitude), 
    #and the second column is the abundance (per Mpc^3 per Mag).
    
    import matplotlib.pyploy as plt
    from AbundanceMatching import *
    af = AbundanceFunction(lf[:,0], lf[:,1], (-27, -5))
    
    #check the abundance function
    plt.semilogy(lf[:,0], lf[:,1])
    x = np.linspace(-27, -5, 101)
    plt.semilogy(x, af(x))
    
    #deconvolution and check results (it's a good idea to always check this)
    scatter = 0.2
    remainder = af.deconvolute(scatter*LF_SCATTER_MULT, 20)
    x, nd = af.get_number_density_table()
    plot(x, remainder/nd);
    
    #get number densities of the halo catalog
    nd_halos = calc_number_densities(halos['vpeak'], box_size)
    
    #do abundance matching with no scatter
    catalog = af.match(nd_halos)
    
    #do abundance matching with some scatter
    catalog_sc = af.match(nd_halos, scatter*LF_SCATTER_MULT)
    
    #if you want to do multiple (100) realizations:
    catalog_deconv = af.match(nd_halos, scatter*LF_SCATTER_MULT, False)
    for __ in xrange(100):
        catalog_this = add_scatter(catalog_deconv, scatter*LF_SCATTER_MULT)
        catalog_this = rematch(catalog_this, catalog, af._x_flipped)
        # do something with catalog_this