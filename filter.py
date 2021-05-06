
import os
import numpy as np
from scipy.stats import median_abs_deviation as MAD
from scipy.spatial.distance import cdist
# from sklearn.preprocessing import MinMaxScaler
from astropy.table import Table
# from astropy.io import ascii
import matplotlib.pyplot as plt

# # THESE PARAMETERS ARE VERY IMPORTANT AND RATHER ARBITRARY
# #
# # Filter pyUPMASK members
# # Filter coordinates by maximum radius (in deg)
# xyRad = {'haf14': 0.1, 'rup41': 0.08, 'rup42': 0.08, 'rup44': 0.1,
#          'rup152': 0.05}
# # Number of STDDEVS to filter in PMs and parallax space
# Nstd = 3

# files = ("haf14",) #"rup41", "rup42", "rup44", "rup152")

# in_file_folder = "../2_pipeline/2_pyUPMASK/out/"
# out_file_folder = '../2_pipeline/3_members_filter/out/'
# out_fig_folder = '../2_pipeline/3_members_filter/tmp/'
# cols_keep = 'EDR3Name _x _y RA_ICRS DE_ICRS Plx e_Plx pmRA e_pmRA pmDE ' +\
#     'e_pmDE Gmag e_Gmag BP-RP e_BP-RP V eV BV eBV UB eUB VI eVI probs_final'
# cols_keep = cols_keep.split()
# # Remove GMM probability columns
# data = data[cols_keep]

rad_memb = {'haf14': (-.01, -.01, 0.083, 200), 'rup41': (0., 0., 0.5, 400),
            'rup42': (0.0165, -0.014, 0.1, 510),
            'rup44': (0, 0.014, 0.05, 100),
            'rup152': (0.01642, -0.02629, 0.05, 150),
            'NGC2516': (0., 0., 1., 1200)}
N_std = 3
print("Standard deviations used in filter", N_std)

in_folder = "./input"
files = os.listdir(in_folder)


def main():
    """
    """
    for file in files:
        print("\n{}".format(file))
        data = Table.read(in_folder + "/" + file, format='ascii')
        print("Total number of stars: {}".format(len(data)))

        fname = file.split('.')[0]

        cx, cy, rad, N_memb = rad_memb[fname]
        print("Target:", N_memb)

        xy = np.array([data['_x'], data['_y']]).T
        c_dist = cdist([(cx, cy)], xy)[0]
        data_rad = data[c_dist < rad]
        print("Rad filter:", len(data_rad))

        Plx_c, pmRA_c, pmDE_c = centers(data_rad)
        print("Plx_c={:.3f}, pmRA_c={:.3f}, pmDE_c={:.3f}".format(
            Plx_c, pmRA_c, pmDE_c))

        prob_range = np.arange(0.5, 1., 0.005)
        rad_memb_pp = []
        for pp in prob_range:
            msk = data_rad['probs_final'] > pp
            data_filt = filterData(data_rad[msk], Plx_c, pmRA_c, pmDE_c)
            N_memb_m = len(data_filt)
            rad_memb_pp.append([abs(N_memb - N_memb_m), N_memb_m])
        rad_memb_pp = np.array(rad_memb_pp).T

        idx = np.argmin(rad_memb_pp[0])
        msk = data_rad['probs_final'] > prob_range[idx]
        data_msk = data_rad[msk]
        data_final = filterData(data_msk, Plx_c, pmRA_c, pmDE_c)

        # msk = data_filt['probs_final'] > prob_range[idx]
        print("Prob cut:", prob_range[idx])

        print("Final:", len(data_final))
        # data_final = data_filt[msk]

        xp, yp = 'pmRA', 'pmDE'  # '_x', '_y'

        plt.subplot(221)
        plt.title("N_memb={}".format(len(data_rad)))
        plt.scatter(data[xp], data[yp], c='grey', alpha=.3)
        plt.scatter(data_rad[xp], data_rad[yp])

        # plt.subplot(222)
        # plt.title("N_memb={}".format(len(data_filt)))
        # plt.scatter(data[xp], data[yp], c='grey', alpha=.3)
        # plt.scatter(data_filt[xp], data_filt[yp])

        plt.subplot(222)
        plt.scatter(prob_range, rad_memb_pp[1])
        plt.axvline(prob_range[idx], c='r')

        plt.subplot(223)
        plt.title("N_memb={}".format(len(data_final)))
        plt.scatter(data[xp], data[yp], c='grey', alpha=.3)
        plt.scatter(data_final[xp], data_final[yp])

        ax = plt.subplot(224)
        plt.scatter(data['BP-RP'], data['Gmag'], c='grey', alpha=.3)
        plt.scatter(data_final['BP-RP'], data_final['Gmag'])
        ax.invert_yaxis()

        plt.show()
        continue

        # msk = data_msk['probs_final'] > .9
        # data_msk = data_msk[msk]

        # print(len(data), len(data_msk))
        # plt.subplot(221)
        # xx = np.linspace(1, 100, 100)
        # yy = np.percentile(data_msk['probs_final'], xx)
        # plt.plot(xx / 100., yy)
        # plt.subplot(222)
        # plt.scatter(data['_x'], data['_y'], c='grey', alpha=.3)
        # plt.scatter(data_msk['_x'], data_msk['_y'])
        # plt.subplot(223)
        # plt.scatter(data['pmRA'], data['pmDE'], c='grey', alpha=.3)
        # plt.scatter(data_msk['pmRA'], data_msk['pmDE'])
        # plt.subplot(224)
        # plt.hist(data['Plx'], 50, color='grey', density=True)
        # plt.hist(data_msk['Plx'], 50, density=True)
        # plt.show()

    # # Store cleaned members
    # fout = out_file_folder + fname + "_" + method + '.dat'
    # ascii.write(memb_d_clean, fout, format='csv', overwrite=True, comment='#')
    # print("Clean members saved to file")

    # makePlot(
    #     method, fname, data, memb_d, field_d, memb_d_clean, pm_Rad, Plx_rad)
    # print("Plot generated")


def centers(data, pmin=.99):
    """
    """
    msk = data['probs_final'] > pmin
    # x_c, y_c = np.median(data[msk]['_x']), np.median(data[msk]['_y'])
    pmRA_c = np.median(data['pmRA'][msk])
    pmDE_c = np.median(data['pmDE'][msk])
    Plx_c = np.median(data['Plx'][msk])

    return Plx_c, pmRA_c, pmDE_c


def filterData(data, Plx_c, pmRA_c, pmDE_c, N_std=N_std):
    """
    """
    d_plx = abs(Plx_c - data['Plx'])
    d_pmRA = abs(pmRA_c - data['pmRA'])
    d_pmDE = abs(pmDE_c - data['pmDE'])

    msk1 = d_plx < N_std * data['Plx'].std()
    msk2 = d_pmRA < N_std * data['pmRA'].std()
    msk3 = d_pmDE < N_std * data['pmDE'].std()
    msk = msk1 & msk2 & msk3

    return data[msk]


def makePlot(
    method, fname, data, memb_d, field_d, memb_d_clean, pm_Rad,
        Plx_rad):
    """
    """
    # Plot
    gs_unit = 5
    gs_x, gs_y = 4, 4
    fig, axs = plt.subplots(gs_y, gs_x, figsize=(
        gs_unit * gs_x, gs_unit * gs_y))

    data_plot = (
        (data['_x'], data['_y'], data['pmRA'], data['pmDE'], data['Plx'],
         data['Gmag'], data['BP-RP']),
        (memb_d['_x'], memb_d['_y'], memb_d['pmRA'], memb_d['pmDE'],
         memb_d['Plx'], memb_d['Gmag'], memb_d['BP-RP']),
        (field_d['_x'], field_d['_y'], field_d['pmRA'], field_d['pmDE'],
         field_d['Plx'], field_d['Gmag'], field_d['BP-RP']),
        (memb_d_clean['_x'], memb_d_clean['_y'], memb_d_clean['pmRA'],
         memb_d_clean['pmDE'], memb_d_clean['Plx'], memb_d_clean['Gmag'],
         memb_d_clean['BP-RP']))

    xmin, xmax = data['_x'].min(), data['_x'].max()
    ymin, ymax = data['_y'].min(), data['_y'].max()
    for ax_y in range(gs_y):
        x, y, pmRA, pmDE, Plx, G, BPRP = data_plot[ax_y]

        color, mrk, sz = None, None, None
        if ax_y in (0, 2):
            color, mrk, sz = 'grey', '.', 4

        ax = axs[ax_y, 0]
        if ax_y in (0, 2):
            ax.set_title("N={}".format(len(x)))
        elif ax_y == 1:
            ax.set_title("N={} | P>={:.2f}".format(len(x), min_prob[fname]))
        elif ax_y == 3:
            ax.set_title("N={} | ({:.4f}, {:.4f})".format(
                len(x), np.median(x), np.median(y)))
        ax.scatter(x, y, alpha=.5, marker=mrk, s=sz, color=color)
        # axs[ax_y, 0].invert_xaxis()
        ax.set_xlim(xmax, xmin)
        ax.set_ylim(ymin, ymax)

        #
        if ax_y in (0, 2):
            cent, std = np.median([pmRA, pmDE], 1),\
                np.mean(MAD([pmRA, pmDE], 1))
            msk = (pmRA < cent[0] + 3 * std) & (pmRA > cent[0] - 3 * std) &\
                (pmDE < cent[1] + 3 * std) & (pmDE > cent[1] - 3 * std)
        else:
            msk = np.array([True for _ in x])
        axs[ax_y, 1].scatter(
            pmRA[msk], pmDE[msk], c=color, alpha=.5, marker=mrk, s=sz)
        axs[ax_y, 1].invert_xaxis()
        if ax_y == 3:
            axs[ax_y, 1].set_title("N={}, PM_rad={:.2f}".format(Nstd, pm_Rad))

        #
        if ax_y == 3:
            axs[ax_y, 2].set_title("N={}, Plx_rad={:.2f}".format(
                Nstd, Plx_rad))
        msk = (Plx > -1) & (Plx < 5)
        axs[ax_y, 2].hist(Plx[msk], 50, color=color)

        #
        axs[ax_y, 3].scatter(
            BPRP, G, c=color, alpha=.5, marker=mrk, s=sz)
        axs[ax_y, 3].grid(ls=':', lw=.7, zorder=-1)
        axs[ax_y, 3].invert_yaxis()

    fig.tight_layout()
    fout = out_fig_folder + fname + "_" + method + '.png'
    plt.savefig(fout, dpi=150, bbox_inches='tight')


if __name__ == '__main__':
    main()
