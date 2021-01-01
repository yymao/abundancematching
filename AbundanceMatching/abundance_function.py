import numpy as np
from scipy.optimize import curve_fit

_error_import_fiducial_deconvolute = None
try:
    from .fiducial_deconv_wrapper import fiducial_deconvolute
except Exception as inst:
    _error_import_fiducial_deconvolute = inst


__all__ = ['AbundanceFunction', 'add_scatter', 'rematch', 'LF_SCATTER_MULT']


LF_SCATTER_MULT = 2.5


def _diff(a):
    return a[1:]-a[:-1]


def _bright_end_func(x, a, b, c, d):
    return -np.exp(a*x+b) + c*x + d


def _convolve_gaussian(y, sigma, truncate=4):
    sd = float(sigma)
    size = int(np.ceil(truncate * sd))
    weights = np.zeros(size*2+1)
    i = np.arange(size+1)
    weights[size:] = np.exp(-(i*i)/(2.0*sd*sd))
    weights[:size] = weights[:size:-1]
    weights /= weights.sum()
    y_full = np.concatenate((np.zeros(size), y, np.ones(size)*y[-1]))
    return np.convolve(y_full, weights, 'valid')


def add_scatter(x, scatter, in_place=False):
    """
    Add a Gaussian scatter to x.

    Parameters
    ----------
    x : array_like
        Values to add scatter to.
    scatter : float
        Standard deviation (sigma) of the Gaussian.
    in_place : bool, optional
        Whether to add the scatter to x in place or return a
        new array.

    Returns
    -------
    x : array_like
        x with the added scatter.
    """
    x = np.asanyarray(x) if in_place else np.array(x)
    r = np.random.randn(*x.shape)
    r *= scatter
    x += r
    return x


def rematch(catalog1, catalog2, greatest_first=True, \
        catalog2_sorted=False):
    """
    Substitute the values in catalog1 with the values in catalog2,
    accroding to the ranks of both arrays. Values of NaN and INF are
    excluded automatically.

    Parameters
    ----------
    catalog1 : array_like
        1-d array in which all the finite values to be substituted by the
        values in catalog2.

    catalog2 : array_like
        1-d array in which the values to be substituted for the values in
        catalog1.

    greatest_first : bool, optional
        If True (default), the assignment starts with the greatest values.

    catalog2_sorted : bool, optional
        If True, do not re-sort catalog2 again.

    Returns
    -------
    catalog : ndarray
        An array that has the same size as catalog1, and all the values are
        substitute by the values in catalog2, according to the ranks.
    """
    arr2 = np.asarray(catalog2)
    if not catalog2_sorted:
        arr2 = arr2[np.isfinite(arr2)]
        arr2.sort()
        if greatest_first:
            arr2 = arr2[::-1]
    arr1 = np.array(catalog1)
    f = np.where(np.isfinite(arr1))[0]
    s = np.argsort(arr1[f])
    if greatest_first:
        s = s[::-1]
    arr1[f[s[:len(arr2)]]] = arr2[:len(s)]
    arr1[f[s[len(arr2):]]] = np.nan
    return arr1


def _to_float(x, default=np.nan):
    try:
        xf = float(x)
    except (ValueError, TypeError):
        return default
    return xf


class AbundanceFunction(object):
    def __init__(self, x, phi, ext_range=(None, None), nbin=1000, \
            faint_end_first=False, faint_end_slope='fit', \
            faint_end_fit_points=3, bright_end_fit_points=-1):
        """
        This class can interpolate and extrapolate an abundance function,
        and also provides fiducial deconvolution and abundance matching.

        Parameters
        ----------
        x : array_like
            The abundance proxy, usually is magnitude or log(stellar mass).
            `log(phi)` should roughly be linear in `x`.
        phi : array_like
            The abundance value, in the unit of x^{-1} vol^{-1}.
            The integrate phi over x should result in number density.
            `x` and `phi` must have the same size.
        ext_range : tuple, optional
            The minimal and maximal value in x to extrapolate abundance
            function.
        nbin : int, optional
            Number of points to interpolate the abundance function.
        faint_end_first : bool, optional
            Whether `x` and `phi` are listed from faint end to bright end.
            If False (default), assumes bright end listed first.
        faint_end_slope : str or float, optional
            If 'fit', fit the faint-end slope from data.
            If a float number, use it as the faint-end slope.
        faint_end_fit_points : int, optional
            Number of points to fit the faint-end slope.
            Only used if `faint_end_slope` is 'fit'.
        bright_end_fit_points : int, optional
            Number of points to fit the bright end.
            If -1 (default), use all data to fit.

        Notes
        -----
        To do abundance matching, see member functions `deconvolute`
        and `match`.
        """
        x = np.ravel(x)
        phi_log = np.log(phi).flatten()

        if len(x) != len(phi_log):
            raise ValueError('`x` and `phi` must have the same size!')

        bright_end_fit_points = min(int(bright_end_fit_points), len(x))
        if bright_end_fit_points < 0:
            bright_end_fit_points = len(x)
        elif bright_end_fit_points < 4:
            raise ValueError('`bright_end_fit_points` must be -1 or larger than 3')

        if faint_end_slope == 'fit':
            faint_end_fit_points = min(int(faint_end_fit_points), len(x))
            if faint_end_fit_points < 2:
                faint_end_fit_points = 0
                faint_end_slope = 0
        else:
            faint_end_slope = float(faint_end_slope)
            faint_end_fit_points = 0

        ext_min, ext_max = ext_range
        ext_min = _to_float(ext_min, x[0])
        ext_max = _to_float(ext_max, x[-1])

        if faint_end_first:
            x = x[::-1]
            phi_log = phi_log[::-1]
            ext_min, ext_max = ext_max, ext_min

        x_new = np.linspace(ext_min, ext_max, num=int(nbin)+1)
        dx = _diff(x)
        if all(dx > 0): #like luminosity
            self._x_flipped = False
            bright_end_flag = (x_new < x[0])
            faint_end_flag = (x_new > x[-1])
        elif all(dx < 0): #like stellar mass
            self._x_flipped = True
            bright_end_flag = (x_new > x[0])
            faint_end_flag = (x_new < x[-1])
        else:
            raise ValueError('x must be a strictly monotonic array.')

        self._s = slice(None, None, -1 if self._x_flipped else None)

        phi_log_new = np.empty_like(x_new)
        flag = ~(bright_end_flag | faint_end_flag)
        phi_log_new[flag] = np.interp(x_new[flag], x[self._s], phi_log[self._s])

        #fit bright end
        a0 = 1.0 if self._x_flipped else -1.0
        s = slice(bright_end_fit_points)
        popt = curve_fit(_bright_end_func, x[s], phi_log[s], [a0, 0, 0, 0], \
                maxfev=100000)[0]
        phi_log_new[bright_end_flag] = \
                _bright_end_func(x_new[bright_end_flag], *popt)

        #fit faint end
        if faint_end_fit_points:
            s = slice(-faint_end_fit_points, None)
            popt = curve_fit(lambda x, a, b: a*x+b, x[s], phi_log[s], [0, 0], \
                    maxfev=100000)[0]
            faint_end_slope = popt[0]
        else:
            faint_end_slope *= (np.log(10.0) if self._x_flipped else -np.log(10.0))
        b = phi_log[-1]-faint_end_slope*x[-1]
        phi_log_new[faint_end_flag] = x_new[faint_end_flag]*faint_end_slope + b

        dx = np.fabs((x_new[-1]-x_new[0])/int(nbin))
        phi_new = np.exp(phi_log_new)
        flag = np.isfinite(phi_new)
        x_new = x_new[flag]
        phi_new = phi_new[flag]

        dphi = _diff(phi_new)
        phi_center = (phi_new[1:]+phi_new[:-1])*0.5
        phi_int = dphi/_diff(phi_log_new)*dx
        flag = (np.fabs(dphi)/phi_center < 1.0e-7)
        if any(flag):
            phi_int[flag] = phi_center[flag]*dx
        phi_int_0 = phi_int[0]*phi_int[0]/phi_int[1]
        phi_int = np.cumsum(np.insert(phi_int, 0, phi_int_0))

        self._x = x_new
        self._phi_log = phi_log_new
        self._nd_log = np.log(phi_int)
        self.nd_bounds = phi_int[0], phi_int[-1]
        self._x_deconv = {}

    def __call__(self, x):
        """
        Return the abundance values at x, i.e. phi(x).

        Parameters
        ----------
        x : array_like
            The abundance proxy, usually is magnitude or log(stellar mass).

        Returns
        -------
        phi : array_like
            The abundance values at x.
        """
        return np.exp(np.interp(x, self._x[self._s], self._phi_log[self._s], \
                np.nan, np.nan))

    def number_density_at(self, x, scatter=0):
        """
        The number density at x, i.e. return nd(x).

        Parameters
        ----------
        x : array_like
            The abundance proxy, usually is magnitude or log(stellar mass).
        scatter : float, optional
            If not zero, it uses an abundance function that has been
            deconvoluted with this amount of scatter.
            Must run `deconvolute` before calling this function.

        Returns
        -------
        nd : array_like
            Number densities.
        """
        scatter = float(scatter)
        if scatter > 0:
            try:
                xp = self._x_deconv[scatter]
            except (KeyError):
                raise ValueError('Please run deconvolute first!')
        else:
            xp = self._x
        return np.exp(np.interp(x, xp[self._s], self._nd_log[self._s], \
                np.nan, np.nan))

    def match(self, nd, scatter=0, do_add_scatter=True, do_rematch=True):
        """
        Abundance matching: match number density to x, i.e. return x(nd).

        Parameters
        ----------
        nd : array_like
            Number densities.
        scatter : float, optional
            If not zero, it uses an abundance function that has been
            deconvoluted with this amount of scatter.
            Must run `deconvolute` before calling this function.
        do_add_scatter : bool, optional
            Add scatter to the final catalog.
        do_rematch : bool, optional
            Rematch the final catalog to the abundance function.

        Returns
        -------
        catalog : array_like
            The abundance proxies (e.g. magnitude or log(stellar mass))
            at the given number densities.
        """
        scatter = float(scatter)
        if scatter > 0:
            try:
                xp = self._x_deconv[scatter]
            except (KeyError):
                raise ValueError('Please run deconvolute first!')
        else:
            xp = self._x
        x = np.interp(np.log(nd), self._nd_log, xp, np.nan, np.nan)

        if scatter > 0:
            if do_add_scatter:
                x = add_scatter(x, scatter, True)
                if do_rematch:
                    x2 = np.interp(np.log(nd), self._nd_log, self._x, np.nan, np.nan)
                    x = rematch(x, x2, self._x_flipped)
        return x

    def deconvolute(self, scatter, repeat=10, sm_step=0.005, return_remainder=True):
        """
        Deconvolute the abundance function with a given scatter (assuming Gaussian)
        This function uses Peter Behroozi's 'fiducial_deconvolute' in c code.
        You must first compile fiducial_deconvolute to use this function.

        Parameters
        ----------
        scatter : float
            Standard deviation (sigma) of the Gaussian, in the unit of x.
        repeat : int, optional
            Number of times to repeat fiducial deconvolute process.
            This value can change the result significantly.
            *Always* check a reasonable value is used.
        sm_step : float, optional
            Some parameter used in fiducial_deconvolute.
            Using 0.01 or 0.005 is fine.
        return_remainder : bool, optional
            If True, calculate the remainder of this deconvolution.
            *Always* check the reminder is reasonable before
            doing abundance matching.

        Returns
        -------
        remainder : array_like
            Returned only if `return_remainder` is True.
        """
        if _error_import_fiducial_deconvolute is not None:
            raise _error_import_fiducial_deconvolute

        af_key = np.empty(len(self._x), float)
        af_val = np.empty_like(af_key)
        af_key[::-1] = self._x
        if not self._x_flipped:
            af_key *= -1.0
        af_val[::-1] = self._phi_log
        af_val /= np.log(10.0)

        smm = np.empty_like(af_key)
        mf = np.empty_like(af_key)
        smm[::-1] = self._x
        mf[::-1] = np.gradient(np.exp(self._nd_log))
        if not self._x_flipped:
            smm *= -1.0
        smm = fiducial_deconvolute(af_key, af_val, smm, mf, scatter, repeat, sm_step)
        if not self._x_flipped:
            smm *= -1.0
        smm = smm[::-1]
        self._x_deconv[float(scatter)] = smm

        if return_remainder:
            nd = np.exp(np.interp(self._x, smm[self._s], self._nd_log[self._s]))
            dx = np.fabs((self._x[-1] - self._x[0])/float(len(self._x)-1))
            nd_conv = _convolve_gaussian(nd, float(scatter)/dx)
            return nd_conv - np.exp(self._nd_log)

    def get_abundance_table(self):
        """
        Return the inter/extrapolated abundance table.

        Returns
        -------
        x : array_like
            Abundance proxy.
        phi : array_like
            Abundance value.
        """
        return self._x, np.exp(self._phi_log)

    def get_number_density_table(self):
        """
        Return the inter/extrapolated number density table.

        Returns
        -------
        x : array_like
            Abundance proxy.
        nd : array_like
            Number density, i.e. int phi(x) dx.
        """
        return self._x, np.exp(self._nd_log)
