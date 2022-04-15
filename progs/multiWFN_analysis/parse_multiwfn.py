#!/usr/bin/env python3

################
## User Input ##
################

# this dict converts orca atom indices to sensible names (required, if orca1.out not present)
# fill here, if automatic determination does not work
orca2name = { 
#    16: 'Fe1',
#    17: 'Fe2',
#    18: 'Fe3',
#    19: 'Fe4',
#    20: 'Fe5',
#    21: 'Fe6',
#    22: 'Fe7',
#    15: 'Mo',
}

import numpy as np
import pandas as pd
import re

from pathlib import Path
import argparse
import sys

try:
    from seaborn import heatmap
    from matplotlib import rcParams
    from matplotlib.pylab import subplots
    rcParams['figure.dpi'] = 200
    rich_output = True
except ModuleNotFoundError:
    rich_output = False
    print('NOTICE: Seaborn/matplotlib not found: only creating plain text files')


#################
## Fuzzy atoms ##
#################
def get_fuzzy_spin(output):
    '''get values of a fuzzy volume integration

    return values of block given below as df (here: spin)
    ATTENTION: MultiWFN counts from 0, returned df counts from 0 

   Atomic space        Value                % of sum            % of sum abs
     1(C )            0.00001069             0.000356             0.000054
     2(H )           -0.00000212            -0.000071            -0.000011
     3(H )            0.00000017             0.000006             0.000001
     4(C )            0.00009353             0.003118             0.000476
     5(H )            0.00002562             0.000854             0.000130
     6(H )            0.00000051             0.000017             0.000003
     7(C )            0.00059285             0.019762             0.003018'''

    df = pd.DataFrame() # empty dataframe
    parse = False # parsing switch

    with open(output) as f:
        for line in f:
            l = line.split()

            if 'Atomic space Value % of sum % of sum abs' in ' '.join(l): # start parsing
                parse = True
                continue
            elif 'Summing up above values:' in line: # stop parsing 
                parse = False
                continue

            if parse:
                i = int(l[0].split('(')[0]) - 1 # counting from 0
                c = float(l[-3])
                df.loc['Fuzzy spin', i] = c

    return df

def get_hirsh_charge(output):
    '''return CM5 and Hirshfeld charges

    parses the block given below and return dataframe with Hirshfeld and CM5 charges
    ATTENTION: MultiWFN counts from 0, returned df counts from 0 

                         ======= Summary of CM5 charges =======
     Atom:    1C   CM5 charge:   -0.227214  Hirshfeld charge:   -0.077632
     Atom:    2H   CM5 charge:    0.079908  Hirshfeld charge:    0.027221
     Atom:    3H   CM5 charge:    0.092847  Hirshfeld charge:    0.041910
     Atom:    4C   CM5 charge:   -0.155276  Hirshfeld charge:   -0.048841'''

    df = pd.DataFrame() # initiate dataframe
    parse = False # parsing switch

    with open(output) as f:
        for line in f:

            if '======= Summary of CM5 charges =======' in line: # start parsing
                parse = True
                continue
            elif 'Summing up all CM5 charges:' in line: # stop parsing
                parse = False
                continue

            if parse:
                l = line.split()

                i = int(re.findall('^\d+', l[1])[0]) - 1 # counting from 0
                hirsh = float(l[-1])
                cm5 = float(l[-4])

                df.loc['Hirshfeld charge', i] = hirsh
                df.loc['CM5 charge', i] = cm5

    return df

##############################
## LIDI.txt and orbcomp.txt ##
##############################

def read_lidi(output):
    '''returns contents of LIDI.txt as dataframes
    
    return
    df_da       delocalization indices for alpha electrons
    df_db       delocalization indices for beta electrons
    df_dt       delocalization indices for all electrons
    df_l        localization indices for alpha, beta, and all electrons'''

    df_da, df_db, df_dt = [ pd.DataFrame(dtype=float) for _ in range(3) ]
    df_la = pd.DataFrame(dtype=float, columns=['a'])
    df_lb = pd.DataFrame(dtype=float, columns=['b'])
    df_lt = pd.DataFrame(dtype=float, columns=['tot'])

    deloc, loc = False, False 
    with open(output) as f:
        for line in f:

            # decide on which block we are in and set df accordingly 
            if 'delocalization' in line.lower():
                deloc = True
                if 'alpha' in line:
                    df = df_da
                elif 'beta' in line:
                    df = df_db
                elif 'Total' in line:
                    df = df_dt
                continue
            elif 'localization' in line.lower():
                loc = True
                if 'alpha' in line:
                    df = df_la
                elif 'beta' in line:
                    df = df_lb
                else:
                    df = df_lt
                continue
            elif len(line.split()) == 0:
                deloc, loc = False, False 
                continue

            # now parse
            else:
                if deloc:
                    l = line.split()

                    try:
                        x = l[1]
                    except IndexError:
                        x = ''

                    if '.' not in x: # if not float, it is column name
                        i_col_range = [ int(i) for i in l ]

                    elif '.' in x: # if float, it is a contribution
                        i_row = int(l[0])
                        c_range = [ float(c) for c in l[1:] ]

                        for i_col, c in zip(i_col_range, c_range): # fill dataframe
                            df.loc[i_row - 1, i_col - 1] = c # start counting from 0

                elif loc:
                    line_ = re.sub('[a-zA-Z\(\):]', '', line)
                    l = line_.split()

                    for l1, l2 in zip(l[::2], l[1::2]):
                        i = int(l1)
                        c = float(l2)

                        df.loc[ i - 1,:] = c # start counting from 0

        df_l = pd.concat([df_la, df_lb, df_lt], axis=1)
        
        return df_da, df_db, df_dt, df_l


def read_orbcomp(output, orb_i=0, orb_f=np.inf, orb_0=0):
    '''creates pandas dataframe from orbcomp.txt, counting from 0

    optional
    orb_i       first orbital index (counting from 0)
    orb_f       last orbital index
    orb_0       index for start counting (default: 0), required for beta'''

    df = pd.DataFrame() # empty DataFrame

    with open(output) as f:
        for line in f: 
            
            if 'Orbital' in line: # get orbital number from line
                l = line.split()
                n = int(l[1]) - 1 # adjust to ORCA counting
                continue

            # dont parse if out of range
            if orb_i > n:
                continue
            elif orb_f < n:
                break

            l = line.split()
            i = int(l[0]) - 1 # adjust to ORCA counting
            c = float(l[1]) / 100 # convert to 0 < c < 1

            df.loc[n - orb_0, i] = c

    return df

def get_locorbs(multi, orbcomp, spin, orca2name, minerg, minsum, thresh):
    '''creates dataframe with orbitals that are localized on certain atoms

    requires
    multi       multiwfn output with "-1 # Print basic information of all orbitals" called
    orbcomp     location of 'ORBCOMP.txt' file
    spin        electron spin. 0 for alpha, 1 for beta
    orca2name   dictionary connection orca index with user defined name
    minerg      ignore orbitals with energies lower than this
    minsum      ignore orbitals with summed contributions (according to orca2name) lower than this
    thresh      remove single atomic contributions lower than this (to enhance readability)'''

    # first: get orbital numbers from multiwfn output
    # determine orbital ranges in multiwfn numbering
    df = get_occ_orbitals(multi)

    df_subset = df.loc[ (df.loc[:,'spin'] == spin) & (df.loc[:,'erg'] > minerg) ]

    orb_i, orb_f = df_subset.index.min(), df_subset.index.max()
    orb_0 = df.loc[ df.loc[:,'spin'] == 1 ].index.min() if spin == 1 else 0

    # second: only get specific orbitals from orbcomp.txt
    df = read_orbcomp(orbcomp, orb_i, orb_f, orb_0) # new dataframe df!
    
    # third: only keep orbitals with certain contributions
    df = rename_columns(df, orca2name)
    
    # delete uninteresting orbitals from dataframe
    df.loc[:, 'sum'] = df.loc[:, : ].sum(axis=1) # create sum column
    df = df.loc[ df['sum'] > minsum, : ]

    # remove small contributions
    df = df.apply( lambda x: [y if y > thresh else np.nan for y in x])
    
    # sort values into blocks with > 0.5
    sort_mask = df.where(df > 0.5).sort_values(by=list( df.columns ), ascending=False) # sort mask
    sort_idx = sort_mask.index # order of row for sorting
    df = df.loc[sort_idx, :] # apply

    # append 'a' or 'b' to number to distinguish spin
    df = df.rename(lambda x: '{}{}'.format(x, 'a' if spin == 0 else 'b'))

    return df

def get_lidi(lidi, orca2name, thresh):
    '''creates dataframe with valences and bond orders / localization and delocalization indices
    for certain atoms

    requires
    lidi        location of 'LIDI.txt' file
    orca2name   dictionary connection orca index with user defined name
    thresh      remove single atomic contributions lower than this (to enhance readability)'''

    df_da, df_db, df_dt, df_l = read_lidi(lidi)

    # delocalization indices
    df_da = df_da.loc[orca2name.keys(),orca2name.keys()]
    df_db = df_db.loc[orca2name.keys(),orca2name.keys()]
    df_dt = df_dt.loc[orca2name.keys(),orca2name.keys()]
    df_da.rename(index=orca2name, columns=orca2name, inplace=True)
    df_db.rename(index=orca2name, columns=orca2name, inplace=True)
    df_dt.rename(index=orca2name, columns=orca2name, inplace=True)
    df_da = df_da.apply( lambda x: [y if y > thresh else np.nan for y in x])
    df_db = df_db.apply( lambda x: [y if y > thresh else np.nan for y in x])
    df_dt = df_dt.apply( lambda x: [y if y > thresh else np.nan for y in x])
    
    # for total: keep only upper triange
    df_dt = df_dt.where(np.triu(np.ones(df_dt.shape)).astype(np.bool), other=np.nan)
    # merge alpha and beta: upper tringle alhpa, lower beta
    df_dab = df_da.where(np.triu(np.ones(df_da.shape)).astype(np.bool), other=-df_db)
    for i in df_dab.columns:
        df_dab.loc[i, i] = np.nan


    # localization indices

    df_l = df_l.loc[orca2name.keys(),:]
    df_l.rename(index=orca2name, inplace=True)
    df_l = df_l.apply( lambda x: [y if y > thresh else np.nan for y in x] )
    
    return df_dt, df_dab, df_l



########################
## General processing ##
########################

def get_orca2name(output):
    '''find Mo or V in orca output file and define names accordingly

    requires
    orcaout     orca output file

    find index of either Mo or V, then defines Fe1, Fe2 etc. as following directly Mo or V
    returns dict with { orca index: name (e.g. 'M', 'Fe1') }'''

    with open(output) as f:

        parse = False
        for line in f:

            if r'CARTESIAN COORDINATES (A.U.)' in line:
                parse = True
                next(f)
                next(f)
                continue

            if parse:
                l = line.split()
                if l[1] == 'Mo' or l[1] == 'V':
                    idx = int(l[0])
                    break
    
    d = { idx + i: 'Fe{}'.format(i) for i in range(1, 8)}
    d[idx] = 'M'

    return d


def get_occ_orbitals(multiout):
    '''parse multiwfn output to get information on occupied alpha and beta orbitals
    ATTENTION: MultiWFN counts from 1, returned values count from 0 (df.index)
    returns df with columns=['spin', 'occ', 'erg']

    requires
    mulitout    output file of multiwfn with "Pirint basic information of all orbitals" called'''

    df = pd.DataFrame() # initiate dataframe

    spindict = { # translate strings in multiwfn to integers
        'Alpha':        0,
        'Beta':         1,
        'Alpha&Beta':   2,}

    with open(multiout) as f:
        for line in f:
            l = line.replace(':', ' ').split()

            if len(l) == 8 and l[0] == 'Orbital' and l[-2] == 'Type':
                    n = int(l[1]) - 1 # adjust to counting from 0
                    e = float(l[3])
                    s = l[-1]
                    o = float(l[-3])

                    if o != 0.0: # only consider occupied orbitals
                        df.loc[n,'spin'] = spindict[s]
                        df.loc[n, 'occ'] = o
                        df.loc[n, 'erg'] = e

    return df


def get_qmatoms(qmatoms):
    '''
    returns a list of indices appearing in the same order of the orca xyz 
    reads the single line file qmatoms
    '''
    with open(qmatoms) as f:
        # get rid of anything but numbers

        line = f.readline().strip()
        line = line.lstrip('set qmatoms {')
        line = line.rstrip('}')

        # make list of chemshell indices
        indices = [ int(i) for i in line.split() ]
        indices.sort()

        return indices

def get_atom_psf(psf, aName, resName=None, resID=None, segName=None, fast=True):    
    '''find atom index in psf file based on atom information

    get atom serial number from psf file (XPLOR format)
    atom number as appearing in psf file (numbering starts with 1)

    first match for atom name [optional: + residue name, +residue id, + segName]

    required: 
    aName           atom name e.g. HE1
    optional: 
    resName         residue name such as MET
    resID           residue ID such as 356
    segName         segment name such as ENZ1
    fast            bool: parsing algorithm: False for regex, True for format string (faster)'''

    match = re.compile( # regex for slow matching
        '^\s*'
        '(?P<i>\d+)'
        '\s*'
        '(?P<sn>[\w\d]+)'
        '\s*'
        '(?P<ri>\d+)'
        '\s*'
        '(?P<rn>[\w\d]+)'
        '\s*'
        '(?P<an>[\w\d]+)'
        '\s*'
        '(?P<at>[\w\d]+)'
        '\s*'
        '(?P<c>-?\d+\.\d+)'
        '\s*'
        '(?P<w>\d+\.\d+)'
        '\s*'
        '(?P<n>\d+)'
        '\s*$'
        )


    with open(psf) as f:

        parse = False
        for line in f:
            if '!NATOM' in line: # start here
                parse = True
                continue
            elif len(line.split()) == 0: # stop on emtpy line
                parse = False
            
            if parse:
                if fast: # fast
                    an = line[38:47].strip() 
                    if an == aName: # only continue if atomname matches
                        i = int(line[0:10].strip())
                        sn = line[11:20].strip()
                        ri = int(line[20:29].strip())
                        rn = line[29:38].strip()
                        if ( rn == resName or resName == None ) and ( ri == resID or resID == None ) and ( sn == segName or segName == None ):
                            return i # terminate once found
                else: # regex
                    if aName in line:
                        d = match.search(line).groupdict()
                        i = int(d['i'])
                        sn = d['sn']
                        ri = int(d['ri'])
                        rn = d['rn']
                        an = d['an']
                        
                        if an == aName and ( rn == resName or resName == None ) and ( ri == resID or resID == None ) and ( sn == segName or segName == None ):
                            return i

        return None # if not found


def rename_columns(df, d):
    '''remove and rename columns in a dataframe

    all coloumns in df that are not keys in d are removed
    all columns in df that are keys in d are renamed according to d
    returns modified dataframe

    df      dataframe
    d       dict with oldcol: newcol'''

    df = df.loc[:,d.keys()]
    df.rename(columns=d, inplace=True)

    return df


#########################
## Plotting and Saving ##
#########################

def write_plaintext(df, path):
    ' write df as plaintext to path '
    with open(path,'w') as f:
        df.to_string(f,columns=df.columns)

def save_fig(fig, path):
    if path:
        fig.savefig(path, transparent=False)

def plt_charge_spin(df_charge, df_spin):
    '''plot charge and spin dataframe in one figure

    required
    df_charge   dataframe with charges
    df_spin     dataframe with spin populations
    path        file name to save plot

    returns nothing'''

    x = len(df_charge.columns) 
    fig, axarr = subplots(nrows=2, figsize=(3*x, 8) )

    kw_args = { # settings for both subplots
        'linewidth': 0.5,
        'yticklabels': True,
        'annot': True, 
        'fmt': '.3f', 
        'annot_kws': { 'fontfamily': 'monospace'}, }

    # plot charge to first axis
    ax = axarr[0]
    ax.set_title('Charge')
    heatmap( df_charge.iloc[0:1,:], ax=ax, cmap='viridis_r', **kw_args )

    # plot spin to second axis
    ax = axarr[1]
    ax.set_title('Spin')
    heatmap( df_spin, ax=ax, cmap='seismic_r', **kw_args )

    for ax in axarr:
        ax.tick_params(axis='both', labelrotation=0)

    fig.tight_layout()

    return fig

def plt_orb(df, cmap='coolwarm'):
    '''create plot for orbital composition
    takes dataframe with orbital compositions and creates a plot at path

    requires
    df          dataframe
    path        path to save figure

    returns nothing'''
    
    
    y, x = .225 * len(df.index), .55 * len(df.columns) # autogenerate figure size
    fig, ax = subplots(figsize=(x, y)) # create figure and axis

    heatmap( # plot seaborn heatmap
        data=df, ax=ax,
        cmap=cmap, vmin=-1, vmax=1, cbar=False,
        xticklabels=True, yticklabels=True, linewidth=0.5,
        annot=True, fmt='.2f', annot_kws={ 'fontfamily': 'monospace'},)

    ax.tick_params(axis='y', which='major')
    ax.tick_params(axis='both', labelrotation=.5)
    ax.tick_params(top=True, labeltop=True)

    fig.tight_layout()

    return fig

def run(orca2name):
    'run from command line'

    # command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('multi_out', help='MultiWFN output file')
    parser.add_argument('-o', '--orca-out', metavar='OUT', help='Orca output file', default='orca1.out')
    parser.add_argument('-e', '--orb_erg', metavar='ERG', default=-10000, type=float, 
        help='Only show orbitals with energy (in Ha) higher than ERG. Default: -10000')
    parser.add_argument('-s', '--orb_sum', metavar='SUM', default=.6, type=float, 
        help='Only show orbitals that have a summed contribution (for atoms defined in orca2name) higher than SUM. Default: 0.6')
    parser.add_argument('-t', '--orb_thresh', metavar='THRESH', default=.02, type=float, 
        help='Remove orbital contributions on a single atom below THRESH. Default: 0.02')
    args = parser.parse_args()

    # define input files
    multi = Path(args.multi_out)
    lidi = multi.parent / 'LIDI.txt' 
    orbcomp = multi.parent / 'orbcomp.txt'
    orca = Path(args.orca_out)

    # define thresholds
    minsum = args.orb_sum
    minerg = args.orb_erg
    thresh = args.orb_thresh

    # start processing
    if multi.is_file():
        print('Now processing {} ...'.format(multi.name))
        # get orca indices, if orca2name dict is empty and orca output exists
        if orca.is_file() and not orca2name:
            print('    ... found {}'.format(orca.name))
            orca2name = get_orca2name(orca)

        # charge
        print('    ... getting charges'.format(multi.name))
        df_charge = get_hirsh_charge(multi.name)
        df_charge = rename_columns(df_charge, orca2name) 

        #spin
        print('    ... getting spin populations'.format(multi.name))
        df_spin = get_fuzzy_spin(multi.name)
        df_spin = rename_columns(df_spin, orca2name) 

        # write output
        df_charge_spin = pd.concat([df_charge, df_spin]) # merge charge and spin


        if rich_output:
            charge_spin_sheet = Path('charge_spin.xlsx')
            charge_spin_fig = Path('charge_spin.png')
            print('    ... writing charges and spin population:\n        {}\n        {}'.format(charge_spin_sheet, charge_spin_fig))
            df_charge_spin.to_excel(charge_spin_sheet)
            fig_charge_spin = plt_charge_spin( df_charge, df_spin)
            save_fig(fig_charge_spin, path=charge_spin_fig)
        else:
            charge_spin_sheet = Path('charge_spin.txt')
            print('    ... writing charges and spin population:\n        {}'.format(charge_spin_sheet))
            write_plaintext(df_charge_spin, charge_spin_sheet)

        print('    ... done!')

        # orbital compositions
        if orbcomp.is_file():
            print('Now processing {} ...'.format(orbcomp.name))
            print('    ... reading {}'.format(orbcomp.name))
            df_orba = get_locorbs(multi, orbcomp, spin=0, orca2name=orca2name, minerg=minerg, minsum=minsum, thresh=thresh)
            df_orbb = get_locorbs(multi, orbcomp, spin=1, orca2name=orca2name, minerg=minerg, minsum=minsum, thresh=thresh)

            # excel sheet and figure
            df_orbab = pd.concat([ df_orba, -df_orbb ]) # concatenate a and b, make beta negative

            # sort values into blocks with > 0.5
            sort_index = df_orbab.where(np.abs(df_orbab) > .5).sort_values(by=list(df_orbab.columns), ascending=False).index
            df_orbab = df_orbab.loc[sort_index, :]
            if rich_output:
                orb_sheet = Path('orbital_composition.xlsx')
                orb_fig = Path('orbital_composition.png')
                print('    ... writing orbital compositions:\n        {}\n        {}'.format(orb_sheet, orb_fig))
                df_orbab.to_excel(orb_sheet)
                fig_orb = plt_orb(df=df_orbab)
                save_fig(fig_orb, path=orb_fig)
            else:
                orb_sheet = Path('orbital_composition.txt')
                print('    ... writing orbital compositions:\n        {}'.format(orb_sheet))
                write_plaintext(df_orbab.fillna(''), orb_sheet)

            print('    ... done!')
        else:
            print('{} NOT found. Orbitals will not be analyzed'.format(orbcomp.name))

        # valences bond orders / localization and delocalization indices
        if lidi.is_file():
            print('Now processing {} ...'.format(lidi.name))
            print('    ... reading {}'.format(lidi.name))
            df_dt, df_dab, df_l = get_lidi(lidi=lidi, orca2name=orca2name, thresh=thresh)
            
            l_sheet = Path('localization_indices.txt')
            print('    ... writing localization indices / atomic valence:\n        {}'.format(l_sheet))
            write_plaintext(df_l.fillna(''), l_sheet)

            dab_sheet = Path('delocalization_indices_ab.txt')
            print('    ... writing alpha and beta delocalization indices / bond orders:\n        {}'.format(dab_sheet))
            write_plaintext(df_dab.fillna(''), dab_sheet)

            dt_sheet = Path('delocalization_indices_tot.txt')
            print('    ... writing total delocalization indices / bond orders:\n        {}'.format(dt_sheet))
            write_plaintext(df_dt.fillna(''), dt_sheet)

        else:
            print('{} NOT found. Localization and delocalization indices will not be analyzed'.format(lidi.name))


    else: 
        print('{} NOT found. Cannot continue and will exit now'.format(multi.name))
        sys.exit()

if __name__ == '__main__':
    run(orca2name)
