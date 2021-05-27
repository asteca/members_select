
import os
import configparser
from astropy.io import ascii
from astropy.table import Table, vstack

in_folder = "./input"
out_folder = "./output"


def readParams(fname):
    """
    """
    print("\n{}".format(fname))

    in_params = configparser.ConfigParser()
    in_params.read('params.ini')

    gen_pars = in_params['General parameters']
    outl, out_plot = gen_pars.get('outl_method'),\
        gen_pars.getboolean('out_plot')

    data_columns = in_params['Clusters data']

    cx, cy, cPlx, cpmRA, cpmDE, N_memb, pmin = ['a'] * 7
    clpars = data_columns.get(fname)
    if clpars is not None:
        pars = []
        for pp in clpars.split():
            try:
                pars.append(float(pp))
            except ValueError:
                pars.append(pp)
        cx, cy, cPlx, cpmRA, cpmDE, N_memb, pmin = pars
    else:
        print("Cluster not found in 'params.ini' file")

    if outl not in ('2DE', '3DE', 'LOF', 'IF', 'EE', 'HT2', 'SPE'):
        raise ValueError("Outlier method not recognized: {}".format(outl))

    return out_plot, cx, cy, cPlx, cpmRA, cpmDE, N_memb, pmin, outl


def readFiles():
    """
    """
    files = os.listdir(in_folder)
    files.remove("dont_read")

    # Return a case-insensitive sorted list
    return sorted(files, key=str.casefold)


def readData(file):
    """
    """
    data = Table.read(in_folder + "/" + file, format='ascii')
    print("Total number of stars: {}".format(len(data)))
    return data


def writeData(outl, file, memb_d, field_d):
    """
    """
    # Remove added columns
    memb_d.remove_columns(['dist_c'])
    field_d.remove_columns(['dist_c'])
    if outl in ('2DE', '3DE'):
        memb_d.remove_columns(['dc_plx', 'dc_pmra', 'dc_pmde'])
        field_d.remove_columns(['dc_plx', 'dc_pmra', 'dc_pmde'])
    # Identify members and field stars
    memb_d.add_column(1, name='membs_select')
    field_d.add_column(0, name='membs_select')

    fout = out_folder + '/membs_' + file
    ascii.write(vstack([memb_d, field_d]), fout, format='csv', overwrite=True,
                comment='#')

    # if out_field:
    #     field_d.remove_columns(['dist_c'])
    #     if outl in ('2DE', '3DE'):
    #         field_d.remove_columns(['dc_plx', 'dc_pmra', 'dc_pmde'])
    #     field_d.add_column(0, name='membs_select')

    #     fout = out_folder + '/field_' + file
    #     ascii.write(field_d, fout, format='csv', overwrite=True, comment='#')
