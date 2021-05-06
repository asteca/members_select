
import os
import numpy as np
from scipy.stats import median_abs_deviation as MAD
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator
from astropy.table import Table


# Names of required columns in the input data file
col_names = {
    'ra': '_x', 'dec': '_y', 'pmra': 'pmRA', 'pmde': 'pmDE',
    'plx': 'Plx'}

# This parameter determines where the "probability cut" is performed, i.e,
# which stars are considered to be real members
prob_cut = 0.99


def main():
    """
    """
    in_folder = "./input"
    files = os.listdir(in_folder)

    for file in files:
        print("\n{}".format(file))
        data = Table.read(in_folder + "/" + file, format='ascii')
        print("Total number of stars: {}".format(len(data)))

        fname, fext = file.split('.')
        probs_mean = data['probs_final']

        makePlot(data, probs_mean, fname)


def makePlot(data_all, probs_mean, fname):
    """
    Make plots using the final probabilities in the "_probs.dat" files.
    """
    # prob_cut = min(prob_cut, probs_mean.max())
    msk_memb = probs_mean >= min(prob_cut, probs_mean.max())

    fig = plt.figure(figsize=(15, 10))
    G = gridspec.GridSpec(4, 3)
    ax1 = plt.subplot(G[0:2, 0])
    ax21 = plt.subplot(G[0:1, 1])
    ax22 = plt.subplot(G[1:2, 1])
    ax3 = plt.subplot(G[0:2, 2])
    ax4 = plt.subplot(G[2:4, 0])
    ax5 = plt.subplot(G[2:4, 1])
    ax6 = plt.subplot(G[2:4, 2])

    ax1.hist(probs_mean, 50, alpha=.5)
    ax1.axvline(prob_cut, c='g', zorder=6)
    ax12 = ax1.twinx()
    ax12.yaxis.set_minor_locator(AutoMinorLocator(4))
    xx, yy = np.linspace(probs_mean.min(), probs_mean.max(), 100), []
    for pr in xx:
        yy.append((probs_mean >= pr).sum())
    ax12.plot(xx, yy, c='r', lw=3, label=r"$P_{{cut}}$={:.2f}".format(
        prob_cut))
    ax12.set_ylabel(r'$N_{stars}$')
    # ax12.set_yscale('log')
    # ax12.set_yticks([])
    ax1.set_xlabel('Probabilities')
    ax1.set_ylabel(r'$\log (N)$')
    ax1.set_yscale('log')
    ax1.set_xlim(0, 1.01)
    plt.legend()

    plx_prob, pmRA_prob, pmDE_prob = [], [], []
    for pp in np.arange(.5, .99, .01):
        msk = probs_mean >= pp
        plx_prob.append([
            pp, data_all[col_names['plx']][msk].mean(),
            data_all[col_names['plx']][msk].std()])
        pmRA_prob.append([
            pp, data_all[col_names['pmra']][msk].mean(),
            data_all[col_names['pmra']][msk].std()])
        pmDE_prob.append([
            pp, data_all[col_names['pmde']][msk].mean(),
            data_all[col_names['pmde']][msk].std()])

    plx_prob, pmRA_prob, pmDE_prob = [
        np.array(_).T for _ in (plx_prob, pmRA_prob, pmDE_prob)]

    ax21.plot(pmDE_prob[0], pmDE_prob[1], c='b', label='pmDE')
    ax21.fill_between(
        pmDE_prob[0], pmDE_prob[1] - pmDE_prob[2], pmDE_prob[1] + pmDE_prob[2],
        color='blue', alpha=0.1)
    ax21.axvline(prob_cut, c='g', zorder=6)
    ax21.set_xlabel(r'$P_{cut}$')
    ax21.set_ylabel(r'$\mu_{\delta}$ [mas/yr]')

    ax22.plot(pmRA_prob[0], pmRA_prob[1], c='r', label='pmRA')
    ax22.fill_between(
        pmRA_prob[0], pmRA_prob[1] - pmRA_prob[2], pmRA_prob[1] + pmRA_prob[2],
        color='red', alpha=0.1)
    ax22.axvline(prob_cut, c='g', zorder=6)
    ax22.set_ylabel(r'$\mu_{\alpha} \cos \delta$ [mas/yr]')
    ax22.set_xlabel(r'$P_{cut}$')

    ax3.plot(plx_prob[0], plx_prob[1])
    ax3.fill_between(
        plx_prob[0], plx_prob[1] - plx_prob[2], plx_prob[1] + plx_prob[2],
        color='k', alpha=0.1)
    ax3.axvline(prob_cut, c='g', zorder=6)
    ax3.set_ylabel(r'$Plx$')
    ax3.set_xlabel(r'$P_{cut}$')

    ax4.set_title("N={}".format(len(data_all[col_names['ra']][msk_memb])))
    ax4.scatter(
        data_all[col_names['ra']][msk_memb],
        data_all[col_names['dec']][msk_memb],
        marker='o', edgecolor='w', lw=.3, zorder=5)
    ax4.plot(
        data_all[col_names['ra']][~msk_memb],
        data_all[col_names['dec']][~msk_memb], 'k.', alpha=0.3, zorder=1)
    ax4.set_xlabel('RA')
    ax4.set_ylabel('DE')
    ax4.set_xlim(
        max(data_all[col_names['ra']]), min(data_all[col_names['ra']]))
    ax4.set_ylim(
        min(data_all[col_names['dec']]), max(data_all[col_names['dec']]))

    # PMs plot
    pmRA_mean = (
        data_all[col_names['pmra']][msk_memb] /
        np.cos(np.deg2rad(data_all[col_names['pmde']][msk_memb]))).mean()
    pmDE_mean = data_all[col_names['pmde']][msk_memb].mean()
    ax5.set_title(r"$(\mu_{\alpha}, \mu_{\delta})=$" +
                  "({:.3f}, {:.3f})".format(pmRA_mean, pmDE_mean))
    ax5.scatter(
        data_all[col_names['pmra']][msk_memb],
        data_all[col_names['pmde']][msk_memb], marker='.',
        edgecolor='w', lw=.1, alpha=.7, zorder=5)
    ax5.plot(
        data_all[col_names['pmra']][~msk_memb],
        data_all[col_names['pmde']][~msk_memb],
        'k+', alpha=0.3, zorder=1)
    ax5.set_xlabel(r'$\mu_{\alpha} \cos \delta$ [mas/yr]')
    ax5.set_ylabel(r'$\mu_{\delta}$ [mas/yr]')
    xmin, xmax = np.percentile(
        data_all[col_names['pmra']][~msk_memb], (25, 75))
    ymin, ymax = np.percentile(
        data_all[col_names['pmde']][~msk_memb], (25, 75))
    pmra_MAD = MAD(data_all[col_names['pmra']][msk_memb], scale='normal')
    pmde_MAD = MAD(data_all[col_names['pmde']][msk_memb], scale='normal')
    xclmin = data_all[col_names['pmra']][msk_memb].mean() - 5. + pmra_MAD
    xclmax = data_all[col_names['pmra']][msk_memb].mean() + 5. + pmra_MAD
    yclmin = data_all[col_names['pmde']][msk_memb].mean() - 5. + pmde_MAD
    yclmax = data_all[col_names['pmde']][msk_memb].mean() + 5. + pmde_MAD
    ax5.set_xlim(max(xmax, xclmax), min(xmin, xclmin))
    ax5.set_ylim(min(ymin, yclmin), max(ymax, yclmax))

    # Plxs plot
    xmin = np.percentile(data_all[col_names['plx']][msk_memb], 1) -\
        3. * data_all[col_names['plx']][msk_memb].std()
    xmax = np.percentile(data_all[col_names['plx']][msk_memb], 95) +\
        3. * data_all[col_names['plx']][msk_memb].std()
    msk1 = np.logical_and.reduce([
        (data_all[col_names['plx']] > -2), (data_all[col_names['plx']] < 4),
        (msk_memb)])
    msk2 = np.logical_and.reduce([
        (data_all[col_names['plx']] > -2), (data_all[col_names['plx']] < 4),
        (~msk_memb)])
    ax6.hist(
        data_all[col_names['plx']][msk1], 20, density=True, alpha=.7, zorder=5)
    ax6.hist(
        data_all[col_names['plx']][msk2], 75, color='k', alpha=0.3,
        density=True, zorder=1)
    plx_mean = data_all[col_names['plx']][msk_memb].mean()
    plx_16, plx_84 = np.percentile(
        data_all[col_names['plx']][msk_memb], (16, 84))
    ax6.axvline(
        plx_mean, c='r',
        label=r"$Plx_{{mean}}={:.3f}_{{{:.3f}}}^{{{:.3f}}}$".format(
            plx_mean, plx_16, plx_84), zorder=6)
    ax6.axvline(plx_16, c='orange', ls='--', zorder=6)
    ax6.axvline(plx_84, c='orange', ls='--', zorder=6)
    ax6.set_xlabel('Plx')
    ax6.set_xlim(xmin, xmax)
    ax6.legend()

    file_out = 'output/' + fname + '.png'
    fig.tight_layout()
    plt.savefig(file_out, dpi=150, bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    main()
