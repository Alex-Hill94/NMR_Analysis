"""
Module Documentation here

Home of the createLineshape function that is used to create the 1D 1H intensity arrays for simulated spectra. All simulation types eventually call this function.
Adapted from nmrsim to allow control of the line function and overall peak width.
"""
#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
_copyright_ = ""
_credits_ = ""
_licence_ = ("")
_reference_ = ("")
#=========================================================================================
# Last code modification:
#=========================================================================================
_modifiedBy_ = "$modifiedBy:  $"
_dateModified_ = "$dateModified: 2023-11-24 11:01:07 +0000 (Fri, November, 24, 2023) $"
_version_ = "$Revision: 3.2.0 $"
#=========================================================================================
# Created:
#=========================================================================================
_author_ = "$Author:  $"
_date_ = "$Date: 2023-11-08 14:40:00 +0000 (Wed, November, 08, 2023) $"
#=========================================================================================
# Start of code
#=========================================================================================


import math
import numbers
import numpy as np


def createLineshape(peaklist, points=65536, limits=None, function='lorentzian'):
    """
    Numpy Arrays of the simulated lineshape for a peaklist.
    Parameters
    ----------
    peaklist : [(float, float, float)...]
        A list of (frequency, intensity, width) tuples.
    y_min : float or int
        Minimum intensity for the plot.
    y_max : float or int
        Maximum intensity for the plot.
    points : int
        Number of data points.
    limits : (float, float)
        Frequency limits for the plot.
    function: string
        Plotting function for the peak shape, either lorentzian or gaussian.
    Returns
    -------
    x, y : numpy.array
        Arrays for frequency (x) and intensity (y) for the simulated lineshape.
    """
    peaklist.sort()
    if limits:
        l_limit, r_limit = low_high(limits)
    else:
        l_limit = peaklist[0][0] - 0.5
        r_limit = peaklist[-1][0] + 0.5
    x_inv = np.linspace(l_limit, r_limit, points)
    x = x_inv[::-1]  # reverses the x axis
    if len(peaklist) > 0:
        if function == 'lorentzian':
            y = add_lorentzians(x, peaklist)
        elif function == 'gaussian':
            y = add_gaussians(x, peaklist)
    else:
        y = np.zeros(points)
    return x, y


def is_number(n):
    if isinstance(n, numbers.Real):
        return n
    else:
        raise TypeError('Must be a real number.')


def is_tuple_of_two_numbers(t):
    m, n = t
    return is_number(m), is_number(n)


def low_high(t):
    two_numbers = is_tuple_of_two_numbers(t)
    return min(two_numbers), max(two_numbers)


def lorentz(v, v0, I, w):
    """
    A lorentz function that takes linewidth at half intensity (w) as a
    parameter.
    When `v` = `v0`, and `w` = 0.5 (Hz), the function returns intensity I.
    Arguments
    ---------
    v : float
        The frequency (x coordinate) in Hz at which to evaluate intensity (y
        coordinate).
    v0 : float
        The center of the distribution.
    I : float
        the relative intensity of the signal
    w : float
        the peak width at half maximum intensity
    Returns
    -------
    float
        the intensity (y coordinate) for the Lorentzian distribution
        evaluated at frequency `v`.
    """
    return I * ((0.5 * w)**2 / ((0.5 * w)**2 + (v - v0)**2))


def add_lorentzians(linspace, peaklist):
    """
    Adapted from nmrsim
    Given a numpy linspace, a peaklist of (frequency, intensity, width)
    tuples, and a linewidth, returns an array of y coordinates for the
    total line shape.
    Arguments
    ---------
    linspace : array-like
        Normally a numpy.linspace of x coordinates corresponding to frequency
        in Hz.
    peaklist : [(float, float)...]
        A list of (frequency, intensity) tuples.
    w : float
        Peak width at half maximum intensity.
    Returns
    -------
    [float...]
        an array of y coordinates corresponding to intensity.
    """
    result = lorentz(linspace, peaklist[0][0], peaklist[0][1], peaklist[0][2])
    for v, i, w in peaklist[1:]:
        result += lorentz(linspace, v, i, w)
    return result


def gauss(v, v0, I, w):
    wf = 0.4246609
    return I * (math.e**((-(v - v0)**2) / (2 * ((w * wf)**2))))


def add_gaussians(linspace, peaklist, w):
    result = gauss(linspace, peaklist[0][0], peaklist[0][1], w)
    for v, i, w in peaklist[1:]:
        result += gauss(linspace, v, i, w)
    return result

