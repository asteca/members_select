
import numpy as np
from scipy.stats import anderson_ksamp, gaussian_kde
import warnings


def runTest(cl_region, fr_region, ad_runs=100):
    """
    """
    pv_plx_cl, pv_pmRA_cl, pv_pmDE_cl = [[] for _ in range(3)]
    for run_num in range(ad_runs):

        data_cl = dataExtract(cl_region)
        data_fr = dataExtract(fr_region)

        # Compare to the field region.
        ad_pv = ADtest(data_cl, data_fr)
        pv_plx_cl.append(ad_pv[0])
        pv_pmRA_cl.append(ad_pv[1])
        pv_pmDE_cl.append(ad_pv[2])

    P1 = KDEProb(pv_plx_cl + pv_pmRA_cl + pv_pmDE_cl)
    # P1 = KDEProb(pv_pmRA_cl + pv_pmDE_cl)

    return P1


def dataExtract(region):
    """
    """
    # Plx + pm_ra + pm_dec
    kins = np.array([region['Plx'], region['pmRA'], region['pmDE']])
    e_kin = np.array([region['e_Plx'], region['e_pmRA'], region['e_pmDE']])
    k_err = []
    for i, k in enumerate(kins):
        # Only process if any star contains at least one not 'nan'
        # data.
        if np.any(~np.isnan(k)):
            k_err.append(normErr(k, e_kin[i]))

    data_all = np.array(k_err)

    return data_all


def normErr(x, e_x):
    # Randomly move mag and color through a Gaussian function.
    return x + np.random.normal(0, 1, len(x)) * e_x


def ADtest(data_x, data_y):
    """
    Obtain Anderson-Darling test for each data dimension.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        ad_vals = []
        # For each dimension
        for i, dd in enumerate(data_x):
            # Store p-value.
            pv0 = anderson_ksamp([dd, data_y[i]])[2]
            ad_vals.append(pv0)

    return ad_vals


def KDEProb(p_vals_cl):
    """
    """
    def kdeLims(xvals):
        xmin, xmax = max(-1., min(xvals)), min(2., max(xvals))
        xrng = (xmax - xmin) * .3
        return xmin - xrng, xmax + xrng

    def regKDE(p_vals):
        xmin, xmax = kdeLims(p_vals)
        x_kde = np.mgrid[xmin:xmax:1000j]
        # Obtain the 1D KDE for this distribution of p-values.
        kernel = gaussian_kde(p_vals)
        #
        prob = abs(kernel.integrate_box_1d(x_kde.min(), .05))
        return prob

    # Ranges for the p-values
    cl_fr_range = np.ptp(p_vals_cl)

    if cl_fr_range > 0.001:
        p_val_cl = regKDE(p_vals_cl)
    else:
        p_val_cl = 1.

    return p_val_cl
