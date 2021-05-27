
import numpy as np
from scipy.spatial.distance import cdist
from astropy.table import vstack
from modules import dataIO, process, plot


def main(pstep=0.005):
    """
    """
    files = dataIO.readFiles()

    for file in files:
        fname = file.split('.')[0]

        # if fname not in ('trumpler_5', 'trumpler_20'):
        #     continue

        out_plot, cx, cy, Plx_c, pmRA_c, pmDE_c, N_memb_input,\
            pmin, outl = dataIO.readParams(fname)

        print("Outlier rejection method:", outl)
        pmin = 0.5 if pmin == 'a' else pmin
        print("Minimum probability used in range:", pmin)
        prob_range = np.arange(1., pmin, -pstep)

        data = dataIO.readData(file)
        # The centers are estimated using a high probability subset
        cx_e, cy_e, Plx_e, pmRA_e, pmDE_e = process.centers(data)

        if cx == 'a':
            cx = cx_e
        if cy == 'a':
            cy = cy_e
        if Plx_c == 'a':
            Plx_c = Plx_e
        if pmRA_c == 'a':
            pmRA_c = pmRA_e
        if pmDE_c == 'a':
            pmDE_c = pmDE_e

        print("cx={:.4f}, cy={:.4f}".format(cx, cy))
        print("Plx_c={:.3f}, pmRA_c={:.3f}, pmDE_c={:.3f}".format(
            Plx_c, pmRA_c, pmDE_c))

        xy = np.array([data['_x'], data['_y']]).T
        data.add_column(cdist([(cx, cy)], xy)[0], name='dist_c')

        # Add distances to Plx+PMs centers. Square here to avoid repeating the
        # operation in filterOutliers()
        if outl in ('2DE', '3DE'):
            data.add_column((data['Plx'] - Plx_c)**2, name='dc_plx')
            data.add_column((data['pmRA'] - pmRA_c)**2, name='dc_pmra')
            data.add_column((data['pmDE'] - pmDE_c)**2, name='dc_pmde')

        rad_memb_pp, idx, min_prob, rad_cl, memb_d, field_d = dataProcess(
            N_memb_input, outl, prob_range, data)

        print("*** {} {:.2f} {:.2f} {:.0f}".format(
            fname, min_prob, rad_cl, len(memb_d)))

        # # Store split data
        # dataIO.writeData(outl, file, memb_d, field_d)

        if out_plot:
            plot.make(
                outl, fname, data, cx, cy, Plx_c, pmRA_c, pmDE_c,
                rad_memb_pp, idx, min_prob, rad_cl, memb_d, field_d)

        print("Finished")


def dataProcess(N_memb_input, outl, prob_range, data):
    """
    """
    if N_memb_input != 'a':
        # TODO finish
        N_memb = int(N_memb_input)
        min_prob, memb_d, field_d = process.manualFilter(
            outl, prob_range, data, N_memb)

        rad_cl = memb_d['dist_c'].max()

        rad_memb_pp = np.array([[0, rad_cl, N_memb, min_prob]]).T
        idx = 0

    else:
        rad_memb_pp = process.autoFilter(outl, prob_range, data)

        # Extract final optimal values
        idx = np.argmin(rad_memb_pp[0])
        rad_cl, min_prob = rad_memb_pp[1][idx], rad_memb_pp[-1][idx]
        # Split into members and field
        msk = data['probs_final'] >= min_prob

        memb_d, data_rjct = process.filterOutliers(outl, data[msk])
        field_d = vstack([data_rjct, data[~msk]])

        msk_in_rad = memb_d['dist_c'] <= rad_cl
        memb_dd = memb_d[msk_in_rad]
        field_dd = vstack([field_d, memb_d[~msk_in_rad]])

    return rad_memb_pp, idx, min_prob, rad_cl, memb_dd, field_dd


if __name__ == '__main__':
    main()
