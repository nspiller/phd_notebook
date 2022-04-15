#!/usr/bin/env python3
import sys

import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pylab as plt
import argparse
from matplotlib import rcParams
rcParams.update({
    'figure.dpi': 200,
    'font.family': 'serif',
    })


def get_data(output):
    '''get LOEWDIN REDUCED ORBITAL POPULATIONS PER MO data from output file
    
    returns list of blocks
    if unrestricted, list[0] contains spin up, list[1] contains spin down
    if restricted, list[0] contains restricted orbitals
    
    
    example unrestricted:
    
    ------------------------------------------                                      
    LOEWDIN REDUCED ORBITAL POPULATIONS PER MO                                      
    -------------------------------------------                                     
    THRESHOLD FOR PRINTING IS 0.1%                                                  
    SPIN UP                                                                         
                          0         1         2         3         4         5       
                     -144.18793 -18.63385 -14.84937 -12.05575 -12.05213 -12.05213   
                       1.00000   1.00000   1.00000   1.00000   1.00000   1.00000    
                      --------  --------  --------  --------  --------  --------    
     0 Ca s             100.0       0.0      99.7       0.0       0.0       0.0     
     0 Ca pz              0.0       0.0       0.0     100.0       0.0       0.0     
     
     
     example restricted:
     
    ------------------------------------------                                    
    LOEWDIN REDUCED ORBITAL POPULATIONS PER MO                                    
    -------------------------------------------                                   
    THRESHOLD FOR PRINTING IS 0.1%                                                
                          0         1         2         3         4         5     
                     -144.18789 -18.63406 -14.84933 -12.05572 -12.05209 -12.05209 
                       2.00000   2.00000   2.00000   2.00000   2.00000   2.00000  
                      --------  --------  --------  --------  --------  --------  
     0 Ca s             100.0       0.0      99.7       0.0       0.0       0.0   
     0 Ca pz              0.0       0.0       0.0     100.0       0.0       0.0   
     0 Ca px              0.0       0.0       0.0       0.0       0.0     100.0  
    '''

    empty_lines = 0 # counter for empty lines
    match = 'LOEWDIN REDUCED ORBITAL POPULATIONS PER MO' # start parsing here
    match_casscf = 'LOEWDIN ORBITAL-COMPOSITIONS' # CASSCF prints different string
    break_casscf = '----------------------------' # break early for CASSCF
    parse = False # switch parsing
    spin = 2 # 0: spin up, 1: spin down, 2: restricted
    data = [ [], [], [] ]

    with open(output) as f:
        for line in f:

            # turn on parsing
            if not parse and ( match in line or match_casscf in line ): 
                parse = True
                next(f) # skip next two lines
                next(f)
                continue
                
            # stop reading
            elif parse and ( not len(line.split()) or break_casscf in line ):
                empty_lines += 1
                # stop if two consecutive empty lines found in spin = 1 or 2 blocks
                if empty_lines == 2: 
                    if spin:
                        break
                    else:
                        continue
                        
            else:
                empty_lines = 0

            # read data
            if parse: 
                if 'SPIN UP' in line:
                    spin = 0
                    continue
                elif 'SPIN DOWN' in line:
                    spin = 1
                    continue
                    
                # writes to: data[0] for alpha, data[1] for beta, data[2] for restricted
                data[spin].append(line)

    # remove empty blocks
    data = [ i for i in data if i ]

    return data


def get_dataframe(block, nmax=None, only_occ=False):
    '''convert plain text data into pandas dataframes
    returns two dataframes, the first containing the contributions, the second orbital information
    
    takes a list containing the plain text data from ORCA output file
    
    optional arguments
        nmax: integer or None, index of highest MO to consider
        only_occ: True or False, stop at first occupation 0.0
    '''        
    
    # get reduced ao basis for df: 
    # dict of dict of set: {a_no1: {l1: {ml11, ml12 ,...}, l2: {ml21, ml22}, }, ... }
    basis = dict()
    # dict with atoms to index list: {'Fe': [0, 1], 'S': [2, 3]}
    atom2index = dict()
    
    for line in block:
        
        l = line[:12].split() # look at fixed width column only
        
        if len(l): # only lines with data in first couple of characters
            a_no = int(l[0])
            a = l[1]
            ao_l = l[-1][0]
            ao_ml = l[-1][1:]

            # fill nested dict, or create subdict 
            if a_no in basis: 
                
                if ao_l in basis[a_no]:
                    basis[a_no][ao_l].add(ao_ml)
                else:
                    basis[a_no][ao_l] = { ao_ml }

            else:
                basis[a_no] = {ao_l: { ao_ml } }
                
            # fill atom2index
            if a in atom2index:
                atom2index[a].add(a_no)
            else:
                atom2index[a] = { a_no }
    
    # ... dict with basis is complete now
    
    # convert to list of tuples for df index          
    columns1 = []
    for k1, v1 in basis.items():
        for k2, v2 in v1.items():
            for v in v2:
                columns1.append((k1, k2, v))
    
    # determine # of mos
    empty_lines = [ i for i, line in enumerate(block) if not len(line.split())]
    last_empty = empty_lines[-2]
    last_mo = int( block[last_empty + 1].split()[-1] )

    # basis for pandas dataframe
    columns1 = pd.MultiIndex.from_tuples(columns1)
    index = range( last_mo + 1)

    # empty numpy array for contributions
    arr1 = np.zeros((len(index), len(columns1)) )
    # flatten multiindex: simple dict with mapping of each index to integer 
    columns2col_idx = { '{}{}{}'.format(*c): i for i, c in enumerate(columns1) }
    
    # empty numpy array for orbital information
    arr2 = np.empty((len(index), 3)) # no, erg, occ
    arr2[:] = np.nan
        
    # start parsing
    parse_rows = False
    header = True
    parse_finished = False
    
    block_iterator = iter(block)
    for line in block_iterator:
        l = line.split()
        
        # check for new header (necessary to capture other than first)
        if len(l) == 0:
            header = True
            continue

        # parse header lines: every 4 lines after empy line
        if header:
            # stop, if premature end
            if parse_finished:
                last_mo = row_idx[-1]
                break
            
            row_idx = [int(i) for i in l] 

            l = next(block_iterator).split()

            es = [ i for i in l ]
            l = next(block_iterator).split()

            occs = [ float(i) for i in l ]
            next(block_iterator)

            # check if reached nmax
            if nmax and nmax in row_idx:
                final_idx = row_idx.index(nmax)
                parse_finished = True

            # check if still occ orbitals
            elif only_occ and 0.0 in occs:
                final_idx = occs.index(0.0)
                parse_finished = True

            if parse_finished:
                row_idx = row_idx[:final_idx+1]
                es = es[:final_idx+1]
                occs = occs[:final_idx+1]
            
            arr2[row_idx, 0] = row_idx
            arr2[row_idx, 1] = es
            arr2[row_idx, 2] = occs

            header = False
            continue

        # parse data lines
        contr = l[3:len(row_idx)+3]
        a_no = l[0]
        ele = l[1]
        ao = l[2]
        ao_l = ao[0]
        ao_ml = ao[1:]
        
        # fill numpy array 
        col_idx = columns2col_idx['{}{}{}'.format(a_no, ao_l, ao_ml)]
        arr1[row_idx, col_idx] = contr
        
            
    # df1 is a multiindex dataframe
    df1 = pd.DataFrame(data=arr1, index=index, columns=columns1)
    df1 = df1.sort_index(axis=1, level=0, sort_remaining=False) # sort columns by atom numbers
    df1 = df1.loc[:last_mo,]
    
    # df2 is a regular 2D dataframe
    columns2 = [ 'mo', 'erg', 'occ' ]
    df2 = pd.DataFrame(data=arr2, index=index, columns=columns2)
    df2.loc[:, 'mo'] = pd.to_numeric(df2.loc[:, 'mo'], downcast='integer', errors='coerce')
    df2 = df2.loc[:last_mo,]

    return df1, df2, atom2index


def select_contr(df, atom2index, mos=[], atoms=[], subshells=[], collapse=0):
    '''Select only part of the dataframe and add atom names
    mos     list of orbital indice (default: all)
    atoms   list of atom indices to 
    return smaller df'''
    
    if not atoms:
        atoms = slice(None)
        print('INFO selecting atoms: all')
    else:
        try:
            atoms = [ atom2index[i] if i in atom2index else int(i) for i in atoms ] # substitute atom names by indices
        except ValueError:
            print('ERROR: One of {} a valid atom'.format(' '.join(atoms)))
            sys.exit(1)

        atoms = list(set(list(pd.core.common.flatten(atoms)))) # flatten possible sets and make sure only unique
        print('INFO selecting atoms: {}'.format(' '.join([ str(i) for i in atoms])))
        
    if not subshells:
        subshells = slice(None)
        print('INFO selecting l: all')
    else:
        print('INFO selecting l: {}'.format(' '.join(subshells)))
        
    if mos:
        mos = range(mos[0], mos[1] + 1)
    else:
        mos = slice(None)
    
    
    df = df.loc[mos,  ( atoms, subshells, )  ]
    
    if collapse == 1:
        df = df.sum(level=(0, 1,), axis=1)
        print('INFO collapsing data: ml')
    elif collapse == 2:
        df = df.sum(level=(0, ), axis=1)
        print('INFO collapsing data: l and ml')
        
    # rename df columns to include atom names
    renamedict = dict()
    for k, v in atom2index.items():
        for i in v:
            renamedict[i] = '{}{}'.format(k, i)
    df = df.rename(columns=renamedict)
        
    return df

def plot(dfarr, path=None, annot=True, s=1.0):
    '''plot list of dataframes as rows
    annot    turn annotations on or off (and cmap off or on)
    s        scaling factor for figure'''
    
    n = len(dfarr)

    y, x = [ i/2/s for i in dfarr[0].shape ]
    
    fig, axarr = plt.subplots(nrows=n, figsize=(x, n*y))
    try:
        len(axarr)
    except TypeError:
        axarr = [ axarr ]

    for df, ax in zip(dfarr, axarr):

        sns.heatmap(data=df, ax=ax, 
                annot=annot, square=True, 
                linewidth=0.01, vmin=0, vmax=100, 
                cbar=not annot, fmt='g', cmap='magma')
        ax.set_ylabel('MO #')
        ax.set_xlabel('AO contr (Atom #-l-ml)')

    fig.tight_layout()
    
    if path:
        fig.savefig(path, transparent=False)

def suggest_active(df, actorbs):
    'compare dataframe with active orbitals and suggest rotation of largest component'
    
    act_i, act_f = actorbs 
    act_n = act_f - act_i + 1

    maxidx = []
    arr = np.array(df)
    for _ in range(act_n):
        row = np.unravel_index(np.argmax(arr, axis=None), arr.shape)[0]
        maxidx.append(row)
        arr[row, :] = -1

    sact = set(range(act_i, act_f + 1))
    smax = set(df.iloc[maxidx, :].index)
    sact, smax = sact - smax, smax - sact

    if sact:
        print('')
        print('rotate')
        for iact, imax in zip(sact, smax):
            print('    {{ {}, {}, 90 }}'.format(imax, iact))
        print('    end')
        print('')
    else:
        print('    ... cannot improve active space')

    
def run():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('orca_output', help='Name of the output file (requires ! normalprint keyword)')
    parser.add_argument('-o', '--output', help='Name of the plot file. Accepts file endigs supported by matplotlib (default orca_output.png)')
    parser.add_argument('-r', '--range', metavar='N', nargs=2, type=int, help='Orbital index range, counting starts with 0. default: all occupied')
    parser.add_argument('-c', '--collapse', metavar='N', type=int, default=0, 
                        help='Sum up contributions by collapsing subshells l (1) and ml (2). default: collapse nothing')
    parser.add_argument('-a', '--atoms', metavar='A', nargs='+',
                        help='Atom names (such as C Fe) or indices (such as 0 1 2). default: all atoms')
    parser.add_argument('-s', '--subshells', metavar='L', nargs='+', type=str,
                        help='Subshells (i.e. orbital angular momentum) to display. default: all (s, p, d, ...)')
    parser.add_argument('--scaling', metavar='S', type=float, default=1.0,
                        help='Scaling factor for text in plot. default: 1')
    parser.add_argument('--numbersoff', action='store_false', 
                        help='Do not plot numbers within plot, but plot colorbar instead')
    parser.add_argument('--suggest', metavar='N', nargs=2, type=int, help='First and last active orbital. Suggest orbital rotations based on largest orbital contributions in selection (only alpha/restricted) for given active orbital range') 

    args = parser.parse_args()
    print('INFO reading data from {}'.format(args.orca_output))
    
    blocks = get_data(args.orca_output)

    if args.range:
        only_occ = False
        nmax = args.range[1]
        print('INFO selecting orbitals: {} to {}'.format(*args.range))
    else:
        only_occ = True
        nmax = None
        print('INFO selecting all orbitals: occupied')
    
    dfarr = []
    for b in blocks:
        df, _, atom2index = get_dataframe(b, nmax=nmax, only_occ=only_occ)
        dfarr.append(df)
        
    # reduce dataframes
    for i, df in enumerate(dfarr):
        print('INFO processing block: {}'.format(i))
        dfarr[i] = select_contr(df, atom2index, mos=args.range, atoms=args.atoms, subshells=args.subshells, collapse=args.collapse)


    # optional: active space suggestion
    if args.suggest:
        print('INFO suggesting orbital rotations for active space ...')
        suggest_active(dfarr[0], args.suggest)
        

    # plot
    output = args.output if args.output else '{}.png'.format(args.orca_output)
    print('INFO plotting: {}'.format(output))
    plot(dfarr, path=output, annot=args.numbersoff, s=args.scaling)
    plt.show()
    


if __name__ == '__main__':
    run()
