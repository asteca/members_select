
import numpy as np
from scipy.stats import median_abs_deviation as MAD
import matplotlib.pyplot as plt


def make(
    outl, fname, data, cx, cy, Plx_c, pmRA_c, pmDE_c, rad_memb_pp, idx,
        min_prob, rad_cl, memb_d, field_d):
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
         field_d['Plx'], field_d['Gmag'], field_d['BP-RP']))

    ax = axs[0, 0]
    msk = data['probs_final'] > .5
    ax.hist(data['probs_final'][msk], 50)
    ax.axvline(min_prob, c='r', ls=':', lw=2)
    ax.set_yscale('log')
    ax.set_xlabel("Probability")
    ax.set_ylabel(r'$\log (N)$')

    prob_range = rad_memb_pp[-1]
    #
    ax = axs[0, 1]
    ax.scatter(prob_range, rad_memb_pp[0])
    ax.axvline(min_prob, c='r', ls=':', lw=2)
    ax.set_title("Prob cut={:.3f}".format(min_prob))
    ax.set_xlabel("Probability")
    ax.set_ylabel("Delta")
    #
    ax = axs[0, 2]
    ax.scatter(prob_range, rad_memb_pp[1] * 60)
    ax.axvline(min_prob, c='r', ls=':', lw=2)
    ax.set_title("Radius={:.1f}".format(rad_memb_pp[1][idx] * 60))
    ax.set_ylabel("Radius [arcmin]")
    ax.set_xlabel("Probability")
    #
    ax = axs[0, 3]
    ax.scatter(prob_range, rad_memb_pp[2])
    ax.axvline(min_prob, c='r', ls=':', lw=2)
    ax.set_title("N_members={:.0f}".format(rad_memb_pp[2][idx]))
    ax.set_ylabel("N members")
    ax.set_xlabel("Probability")

    xmin, xmax = data['_x'].min(), data['_x'].max()
    ymin, ymax = data['_y'].min(), data['_y'].max()
    for ax_y in range(1, 4):
        x, y, pmRA, pmDE, Plx, G, BPRP = data_plot[ax_y - 1]

        color, mrk, sz = None, None, None
        if ax_y in (1, 3):
            color, mrk, sz = 'grey', '.', 4

        ax = axs[ax_y, 0]
        if ax_y in (1, 3):
            ax.set_title("N={}".format(len(x)))
        elif ax_y == 2:
            ax.set_title("P>={:.3f} | r={:.1f} | N={}".format(
                min_prob, rad_memb_pp[1][idx] * 60, len(x)))
        # elif ax_y == 3:
        #     ax.set_title("N={} | ({:.4f}, {:.4f})".format(
        #         len(x), np.median(x), np.median(y)))
        ax.scatter(x, y, alpha=.5, marker=mrk, s=sz, color=color)
        ax.scatter(cx, cy, marker='x', s=15, color='r')
        circle = plt.Circle((cx, cy), rad_cl, color='green', fill=False)
        ax.add_artist(circle)
        ax.set_xlim(xmax, xmin)
        ax.set_ylim(ymin, ymax)
        ax.set_xlabel("ra")
        ax.set_ylabel("de")

        #
        if ax_y in (1, 3):
            cent, std = np.median([pmRA, pmDE], 1),\
                np.mean(MAD([pmRA, pmDE], 1))
            msk = (pmRA < cent[0] + 4 * std) & (pmRA > cent[0] - 4 * std) &\
                (pmDE < cent[1] + 4 * std) & (pmDE > cent[1] - 4 * std)
        else:
            msk = np.array([True for _ in x])
        ax = axs[ax_y, 1]
        ax.scatter(
            pmRA[msk], pmDE[msk], c=color, alpha=.5, marker=mrk, s=sz)
        ax.scatter(pmRA_c, pmDE_c, marker='x', c='r', s=15)
        ax.invert_xaxis()
        if ax_y == 2:
            ax.set_title("{}".format(outl))
        ax.set_xlabel("pmRA")
        ax.set_ylabel("pmDE")

        #
        ax = axs[ax_y, 2]
        msk = (Plx > -1) & (Plx < 5)
        ax.hist(Plx[msk], 50, color=color)
        ax.axvline(Plx_c, c='r', ls=':', lw=2)
        ax.set_xlabel("Plx")

        #
        ax = axs[ax_y, 3]
        ax.scatter(
            BPRP, G, c=color, alpha=.5, marker=mrk, s=sz)
        ax.grid(ls=':', lw=.7, zorder=-1)
        ax.invert_yaxis()
        ax.set_xlabel("BP-RP")
        ax.set_ylabel("G")

    fig.tight_layout()
    fout = 'output/' + fname + '.png'
    plt.savefig(fout, dpi=150, bbox_inches='tight')
    # Close to release memory.
    plt.clf()
    plt.close("all")
