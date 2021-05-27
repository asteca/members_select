
import sys
import numpy as np
from astropy.table import vstack
from scipy.stats import median_abs_deviation as MAD
from sklearn.neighbors import LocalOutlierFactor
from sklearn.ensemble import IsolationForest
from sklearn.covariance import EllipticEnvelope
# https://github.com/erdogant/pca
from pca import pca


def centers(data, pmin=.99, pstep=.01, Nmin_p=0.01, Nstd=5):
    """
    Estimate the center of the cluster in parallax and proper motions, using
    a subset of stars with the largest assigned probabilities.

    pmin: starting probability cut value
    pstep: step used to diminish the probability cut
    Nmin_p: minimum percentage of stars required to estimate the centers
    Nstd: number of STDDEVS used to remove outliers
    """
    def rmOutliers(dd):
        """
        Remove obvious outliers from Plx+PMs, before estimating the centers
        """
        d_plx, d_pmra, d_pmde = MAD(dd['Plx']), MAD(dd['pmRA']),\
            MAD(dd['pmDE'])
        msk = (abs(dd['Plx'] - np.median(dd['Plx'])) < Nstd * d_plx) &\
            (abs(dd['pmRA'] - np.median(dd['pmRA'])) < Nstd * d_pmra) &\
            (abs(dd['pmDE'] - np.median(dd['pmDE'])) < Nstd * d_pmde)
        return dd[msk]

    Nmin = max(10, min(100, int(len(data) * Nmin_p)))

    check_flag = True
    while check_flag:
        msk = data['probs_final'] > pmin
        data_clean = rmOutliers(data[msk])
        if len(data_clean) < Nmin:
            pmin -= pstep
            continue
        cx_c = np.nanmedian(data_clean['_x'])
        cy_c = np.nanmedian(data_clean['_y'])
        pmRA_c = np.nanmedian(data_clean['pmRA'])
        pmDE_c = np.nanmedian(data_clean['pmDE'])
        Plx_c = np.nanmedian(data_clean['Plx'])

        # cx_std = np.std(data_clean['_x'])
        # cy_std = np.std(data_clean['_y'])
        # pmRA_std = np.std(data_clean['pmRA'])
        # pmDE_std = np.std(data_clean['pmDE'])
        # Plx_std = np.std(data_clean['Plx'])

        check_flag = False
        break

    if check_flag:
        raise ValueError("Could not estimate the ra, dec, Plx & PMs center")

    return cx_c, cy_c, Plx_c, pmRA_c, pmDE_c


def manualFilter(outl, prob_range, data, N_memb):
    """
    """
    data.add_column(np.sqrt(
        data['dc_plx'] + data['dc_pmra'] + data['dc_pmde']), name='dd')

    for i, pp in enumerate(prob_range):
        msk = data['probs_final'] >= pp
        if msk.sum() > N_memb:

            data_in_filt, data_out_filt = filterOutliers(outl, data[msk])

            if len(data_in_filt) >= N_memb:
                idx = np.argsort(data_in_filt['dd'])

                memb_d = data_in_filt[idx][:N_memb]
                field_d = vstack([
                    data[~msk], data_out_filt, data_in_filt[idx][N_memb:]])

                min_prob = pp
                break

    return min_prob, memb_d, field_d


def autoFilter(outl, prob_range, data):
    """
    """
    rad_memb_pp = []
    for i, pp in enumerate(prob_range):

        delta, rad_max, N_memb = np.inf, 0, 0

        msk = data['probs_final'] >= pp
        if msk.sum() > 0:
            # Select data by probability cut
            data_in_pp = data[msk]

            # Split the data selected by probability cut, into accepted
            # and outliers
            data_in_filt, data_out_filt = filterOutliers(outl, data_in_pp)

            # Number of stars selected by probability cut, and not
            # rejected as outliers. I.e.: members.
            N_memb = len(data_in_filt)
            # Radius of stars selected as members
            rad_max = data_in_filt['dist_c'].max()
            # Cluster's area
            cl_area = np.pi * rad_max**2
            # Total number of stars in cluster's region
            N_in_rad = (data['dist_c'] <= rad_max).sum()

            # Define ring in the 3-4 MAD region
            rad_3MAD = 3 * MAD(data_in_filt['dist_c'], scale='normal')
            rad_4MAD = 4 * MAD(data_in_filt['dist_c'], scale='normal')
            ring = (data['dist_c'] > rad_3MAD) & (data['dist_c'] < rad_4MAD)
            N_in_ring = ring.sum()
            area_ring = np.pi * (rad_4MAD**2 - rad_3MAD**2)

            # Field density for the ring
            field_dens = N_in_ring / area_ring
            # Estimated number of members in cluster's region
            N_memb_estim = N_in_rad - field_dens * cl_area
            # Members difference (minimized parameter)
            delta = abs(N_memb - N_memb_estim)
            # This approach gives similar but slightly worst results
            # delta = abs(1. - N_memb / N_memb_estim)

        rad_memb_pp.append([delta, rad_max, N_memb, pp])
        # updt(len(prob_range), i + 1)

    return np.array(rad_memb_pp).T


def updt(total, progress, extra=""):
    """
    Displays or updates a console progress bar.

    Original source: https://stackoverflow.com/a/15860757/1391441
    """
    barLength, status = 20, ""
    progress = float(progress) / float(total)
    if progress >= 1.:
        progress, status = 1, "\r\n"
    block = int(round(barLength * progress))
    text = "\r[{}] {:.0f}% {}{}".format(
        "#" * block + "-" * (barLength - block),
        round(progress * 100, 0), extra, status)
    sys.stdout.write(text)
    sys.stdout.flush()


def filterOutliers(outl, data):
    """
    """

    # 3D Ellipsoid method
    if outl == "2DE":
        Ns = 2
        msk = _3dEllip(Ns, data)
        return data[msk], data[~msk]
    if outl == "3DE":
        Ns = 3
        msk = _3dEllip(Ns, data)
        return data[msk], data[~msk]

    # shape: (n_samples, n_features)
    arr = np.array([data['Plx'], data['pmRA'], data['pmDE']]).T
    # Tested, makes very little difference. Probably because the distances in
    # PM and Plx are comparable
    from sklearn.preprocessing import MinMaxScaler
    arr = MinMaxScaler().fit(arr).transform(arr)

    if outl == "LOF":
        # Local Outlier Factor
        msk = LocalOutlierFactor().fit_predict(arr) > 0
    elif outl == "IF":
        # Isolation forest
        msk = IsolationForest().fit_predict(arr) > 0
    elif outl == "EE":
        msk = EllipticEnvelope(store_precision=False).fit_predict(arr) > 0

    # # Tested. The results are awful
    # elif outl == "DBSCAN":
    #     from sklearn.cluster import DBSCAN
    #     from sklearn.neighbors import NearestNeighbors
    #     # from matplotlib import pyplot as plt
    #     neighbors = NearestNeighbors(n_neighbors=20)
    #     neighbors_fit = neighbors.fit(arr)
    #     distances, indices = neighbors_fit.kneighbors(arr)
    #     distances = np.sort(distances, axis=0)
    #     distances = distances[:, 1]
    #     # plt.plot(distances)
    #     eps = rotate(np.array([np.arange(len(distances)), distances]).T)
    #     msk = DBSCAN(eps=eps).fit_predict(arr) == -1

    if outl in ('HT2', 'SPE'):
        # Initialize model.
        # alpha : float, default: 0.05
        #     Alpha to set the threshold to determine the outliers based on
        #     the Hoteling T2 test.
        # n_std : int, default: 2
        #     Number of standard deviations to determine the outliers using
        #     SPE/DmodX method.
        model = pca(alpha=0.05, n_std=2)
        # Fit transform
        out = model.fit_transform(arr, verbose=0)
        if outl == 'HT2':
            # Hotelling T^2
            msk = ~(out['outliers']['y_bool'].values)
        elif outl == 'SPE':
            # SPE/DmodX
            msk = ~(out['outliers']['y_bool_spe'].values)

    return data[msk], data[~msk]


def rotate(data):
    """
    Rotate a 2d vector.

    (Very) Stripped down version of the great 'kneebow' package, by Georg
    Unterholzner. Source:

    https://github.com/georg-un/kneebow

    data   : 2d numpy array. The data that should be rotated.
    return : probability corresponding to the elbow.
    """
    # The angle of rotation in radians.
    theta = np.arctan2(
        data[-1, 1] - data[0, 1], data[-1, 0] - data[0, 0])

    # Rotation matrix
    co, si = np.cos(theta), np.sin(theta)
    rot_matrix = np.array(((co, -si), (si, co)))

    # Rotate data vector
    rot_data = data.dot(rot_matrix)

    # Find elbow index
    elbow_idx = np.where(rot_data == rot_data.min())[0][0]

    # Adding a small percentage to the selected probability improves the
    # results by making the cut slightly stricter.
    prob_cut = data[elbow_idx][1]

    return prob_cut


def _3dEllip(Ns, dd):
    """
    Filter stars outside of the 3D ellipsoid defined in Plx+PMs, centered
    on the center coordinates. The axis lengths of the ellipsoid are given
    by a multiple of the MAD.
    """
    a = Ns * MAD(dd['Plx'], scale='normal')
    b = Ns * MAD(dd['pmRA'], scale='normal')
    c = Ns * MAD(dd['pmDE'], scale='normal')
    msk = dd['dc_plx'] / a**2 + dd['dc_pmra'] / b**2\
        + dd['dc_pmde'] / c**2 < 1
    return msk
