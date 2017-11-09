import numpy as np

__all__ = ['HaloAbundanceFunction', 'calc_number_densities', 'calc_number_densities_in_bins']


def _to_float(x, default=np.nan, take_log=False):
    try:
        xf = float(x)
    except (ValueError, TypeError):
        return default
    return np.log(xf) if take_log else xf


def calc_number_densities(x, box_size):
    """
    Calculate the number densities for a list of values.
    Number density =  rank / volumn.

    Parameters
    ----------
    x : array_like
        An 1-d array that contains the values of a halo property
        (e.g. vpeak, mpeak).
    box_size : float
        The length of the cubic cosmological box.

    Returns
    -------
    nd : array_like
        Number density for each x.
    """
    N = len(x)
    inv_vol = 1.0/(box_size**3)
    nd = np.empty(N, dtype=np.float64)
    nd[np.argsort(x)] = np.linspace(N*inv_vol, inv_vol, N)
    return nd


def calc_number_densities_in_bins(x, box_size, bins):
    """
    Given a list of values, calculate the number densities of
    at the positions of `bins`.
    Number density =  rank / volumn.

    Parameters
    ----------
    x : array_like
        An 1-d array that contains the values of a halo property
        (e.g. vpeak, mpeak).
    box_size : float
        The length of the cubic cosmological box.
    bins : array_like
        The position (in x) to calculate the number densities for.

    Returns
    -------
    nd : array_like
        Number density for each value in `bins`.
    """
    return (len(x) - np.searchsorted(x, bins, sorter=np.argsort(x))).astype(float)/(box_size**3)


class HaloAbundanceFunction(object):
    def __init__(self, x, box_size, fit_range=(None, None), fit_points=10,\
            nbin=200, log_bins=True):
        """
        Create (and extrapolate) a halo abundance function from
        a halo property (e.g. vpeak, mpeak).

        Parameters
        ----------
        x : array_like or str
            An 1-d array that contains the values of a halo property
            (e.g. vpeak, mpeak).
        box_size : float
            The length of the cubic cosmological box.
        fit_range : tuple
            The range (a, b) to extrapolate the halo abundance function.
            The extrapolation starts below b and fits down to a.
            Should be in the same unit as x.
            If (None, None), do not extrapolate.
        fit_points : int
            Number of bins to do the extrapolation mentioned above.
        nbin : int, optional
            Number of bins to interpolate the halo abundance function.
        log_bins : bool, optional
            Whether to take log of the halo property. Default: True.
        """
        x = np.log(x) if log_bins else np.asarray(x)
        x = x[np.isfinite(x)]

        fit_to, fit_below = fit_range
        fit_to = _to_float(fit_to, x.min(), log_bins)
        fit_below = _to_float(fit_below, x.min(), log_bins)

        bins = np.linspace(min(x.min(), fit_to), x.max(), int(nbin)+1)
        nd = calc_number_densities_in_bins(x, box_size, bins)

        if fit_to < fit_below:
            dlog_nd = np.gradient(np.log(nd))
            dx = (bins[-1]-bins[0])/int(nbin)
            k = np.searchsorted(bins, fit_below)
            s = slice(k, k+fit_points)
            self._slope = dlog_nd[s].mean()/dx
            nd[:k] = np.exp((bins[:k]-bins[k])*self._slope) * nd[k]

        self._log_bins = log_bins
        self._x = bins
        self._nd_log = np.log(nd)

    def __call__(self, x):
        """
        Return the number density at x.

        Parameters
        ----------
        x : array_like
            The halo abundance proxy, e.g. vpeak or mpeak.

        Returns
        -------
        nd : array_like
            Number densities at x.

        """
        x = np.log(x) if self._log_bins else np.asarray(x)
        return np.exp(np.interp(x, self._x, self._nd_log, np.nan, np.nan))

    def get_number_density_table(self):
        """
        Return the inter/extrapolated number density table.

        Returns
        -------
        x : array_like
            The halo abundance proxy.
        nd : array_like
            Number densities at x.
        """
        return np.exp(self._x), np.exp(self._nd_log)
