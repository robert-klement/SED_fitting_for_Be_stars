import numpy as np
import os

# COLOR CORRECTION


def color_corr(lbd, flux, flag):
    # path to bandpasses transmission files
    dir0 = os.path.join(os.getcwd(), 'bandpasses')
    # either select filter files, or apply default correction
    # for single fluxes
    if flag.lower() == 'akari':
        #       if type(lbd) == float:
        if lbd.size == 1:
            # print('Just one AKARI data point: assuming slope=3')
            K = np.array([1.096, 0.961, 1.017, 1.203, 0.943, 0.990])
            lbdi = np.array([9., 18., 65., 90., 140., 160.])
            flux_corr = flux / K[np.where(np.abs(lbd - lbdi)
                                          == np.min(np.abs(lbd - lbdi)))]
            return flux_corr
        else:
            band = []
            lbd_arr = np.array(['9', '18', '65', '90', '140', '160'])
            fnames = ['irc-s9w.dat', 'irc-l18w.dat', 'fis-n60.dat', 'fis-wide-s.dat', 'fis-wide-l.dat', 'fis-n160.dat']
            for i in range(len(lbd_arr)):
                delt = 1
                if (np.abs(lbd - float(lbd_arr[i])) < delt).any():
                    band.append(os.path.join(dir0, fnames[i]))

    elif flag.lower() == 'iras':
        if lbd.size == 1:
            # print('Just one IRAS data point: assuming slope=3')
            K = np.array([1.25, 1.23, 1.15])  # for alpha = 1.0 in f_nu propto nu^alpha
            lbdi = np.array([12., 25., 60.])
            flux_corr = flux / K[np.where(np.abs(lbd - lbdi)
                                          == np.min(np.abs(lbd - lbdi)))]
            return flux_corr
        else:
            band = []
            lbd_arr = np.array(['12', '25', '60'])
            delt = 9.
            for i in range(len(lbd_arr)):
                if (np.abs(lbd - float(lbd_arr[i])) < delt).any():
                    band.append(os.path.join(
                        dir0, 'IRAS_IRAS.' + lbd_arr[i] + 'mu.dat'))
    elif flag.lower() == 'wise':
        if lbd.size == 1:
            # print('Just one WISE data point: assuming slope=3')
            K = np.array([0.9961, 0.9976, 0.9393, 0.9934])
            lbdi = np.array([3.3526, 4.6028, 11.5608, 22.0883])
            flux_corr = flux / K[np.where(np.abs(lbd - lbdi)
                                          == np.min(np.abs(lbd - lbdi)))]
            return flux_corr
        else:
            band = []
            lbdi = np.array([3.3526, 4.6028, 11.5608, 22.0883])
            delt = 1.
            for i in range(4):
                if (np.abs(lbd - lbdi[i]) < delt).any():
                    band.append(os.path.join(
                        dir0, 'RSR-W' + '{:1d}'.format(i + 1) + '.txt'))
    else:
        print('unknown input flag')
        flux_corr = np.zeros(len(flux))
        return flux_corr

    # iterative color-correction
    nband = len(band)
    K = np.ones(nband)
    K1 = np.zeros(nband)
    delt = 1e-5
    res = 1.
    while res > delt:
        for iband in range(nband):
            tab = np.loadtxt(band[iband])
            lbdi, R = tab[:, 0], tab[:, 1]
            logG = poly_interp(np.log(lbd), np.log((flux / K)
                                                   / (flux[iband] / K[iband])), np.log(lbdi))
            G = np.exp(logG)
            if flag.lower() == 'wise':
                K1[iband] = integral(lbdi, R * G * lbdi) \
                    / integral(lbdi, R * lbdi)
            else:
                K1[iband] = integral(lbdi, R * G) \
                    / integral(lbdi, R * (lbd[iband] / lbdi))

        res = np.sum(np.abs(K - K1) / K1)
        K = K1

    flux_corr = flux / K

    return flux_corr

# INTEGRAL SIMPSON TRAPEZOID


def integral(x, f, cummulative=False):
    '''
    Integration using trapezoidal rule

    Usage:
    integ = integral(x, f, cummulative=False)
    '''
    x = np.array(x).reshape((-1))
    f = np.array(f).reshape((-1))
    if cummulative:
        integ = np.array([np.trapz(f[0:i], x=x[0:i]) for i in range(len(x))])
    else:
        integ = np.trapz(f, x=x)

    return integ

# POLYNOMIAL N-TH ORDER INTERPOLATION


def poly_interp(xi, yi, x):
    '''
    For N pair of points, interpolates a polynomial
    of (N-1) order

    xi, yi: original data points
    x , y : interpolated values (x is an input)

    Usage:
    y = poly_interp(xi, yi, x)
    '''
    # make sure they are numpy arrays
    xi = np.array([xi]).reshape((-1))
    yi = np.array([yi]).reshape((-1))
    x = np.array([x]).reshape((-1))

    # definitions
    n = len(xi)
    y = 0.

    # loop
    for i in range(n):
        num = yi[i]
        den = 1.
        for j in range(n):
            if j != i:
                num = num * (x - xi[j])
                den = den * (xi[i] - xi[j])

        y = y + num / den

    return y


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
