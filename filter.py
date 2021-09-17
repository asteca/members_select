
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

        out_field, out_plot, cx, cy, Plx_c, pmRA_c, pmDE_c, rad_input,\
            N_memb_input, pmin, outl = dataIO.readParams(fname)

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
        if rad_input != 'a':
            print("Manual radius: {}".format(rad_input))
        if N_memb_input != 'a':
            print("Manual number of members: {}".format(N_memb_input))

        print("cx={:.4f}, cy={:.4f}".format(cx, cy))
        print("Plx_c={:.3f}, pmRA_c={:.3f}, pmDE_c={:.3f}".format(
            Plx_c, pmRA_c, pmDE_c))

        xy = np.array([data['_x'], data['_y']]).T
        data.add_column(cdist([(cx, cy)], xy)[0], name='dist_c')

        # Add distances to Plx+PMs centers. Square here to avoid repeating the
        # operation in filterOutliers()
        if outl in ('2DE', '3DE') or N_memb_input != 'a':
            data.add_column((data['Plx'] - Plx_c)**2, name='dc_plx')
            data.add_column((data['pmRA'] - pmRA_c)**2, name='dc_pmra')
            data.add_column((data['pmDE'] - pmDE_c)**2, name='dc_pmde')

        rad_memb_pp, idx, min_prob, rad_cl, memb_d, field_d = dataProcess(
            rad_input, N_memb_input, outl, prob_range, data, (cx, cy))

        # Store split data
        dataIO.writeData(outl, N_memb_input, file, memb_d, field_d, out_field)

        if out_plot:
            plot.make(
                outl, fname, data, cx, cy, Plx_c, pmRA_c, pmDE_c,
                rad_memb_pp, idx, min_prob, rad_cl, memb_d, field_d)

        print("Finished")


def dataProcess(rad_input, N_memb_input, outl, prob_range, data, center):
    """
    """
    if N_memb_input == 'a' and rad_input == 'a':
        # TODO
        print("FINISH THIS")
        print("Using AD filter")
        # rad_memb_pp = process.autoFilter(outl, prob_range, data)
        rad_memb_pp = process.ADFilter(outl, prob_range, data, center)

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

    else:

        if N_memb_input != 'a' and rad_input == 'a':
            # TODO finish
            N_memb = int(N_memb_input)
            min_prob, memb_dd, field_dd = process.manualFilterNmemb(
                outl, prob_range, data, N_memb)
            rad_cl = memb_dd['dist_c'].max()

        if N_memb_input == 'a' and rad_input != 'a':
            rad_m = float(rad_input)
            min_prob, memb_dd, field_dd = process.manualFilterRad(
                outl, prob_range, data, rad_m)
            rad_cl = rad_m

        elif N_memb_input != 'a' and rad_input != 'a':
            rad_m = float(rad_input)
            N_memb = int(N_memb_input)
            min_prob, memb_dd, field_dd = process.manualFilterNmembRad(
                outl, prob_range, data, N_memb, rad_m)
            rad_cl = rad_m

        rad_memb_pp = np.array([[0, rad_cl, N_memb, min_prob]]).T
        idx = 0

    return rad_memb_pp, idx, min_prob, rad_cl, memb_dd, field_dd


if __name__ == '__main__':
    main()
