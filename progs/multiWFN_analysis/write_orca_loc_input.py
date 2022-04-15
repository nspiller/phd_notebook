#!/usr/bin/env python3

import pandas as pd
import argparse
import re

def get_orbitals(output):
    '''return orbital information as pandas dataframe

    takes one argument: orca output file
    returns pandas dataframe of the last occuring ORBITAL ENERGY section'''

    parse, alpha = False, False # initial switches

    with open(output) as f:
        for line in f:
            l = line.split()

            if 'SPIN UP ORBITALS' in line: # (re)initial dataframe here 
                df = pd.DataFrame(columns=['n', 'occ', 'e_ha', 'e_ev', 'spin'])
                parse, alpha = True, True # parse alpha orbitals
                continue
            elif 'SPIN DOWN ORBITALS' in line:
                parse, alpha = True, False #swtich to beta orbitals
                continue
            elif len(l) == 0: # termitate parsing on emtpy line
                parse = False
                continue
            elif l[0] == 'NO': # skip header line
                continue

            if parse:
                n = int(l[0])         # orbital numer
                occ = float(l[1])     # occupation
                ha = float(l[2])      # energy Ha
                ev = float(l[3])      # energy eV
                s = 0 if alpha else 1 # spin
                
                row = len(df.index)   # determine row number
                df.loc[row,:] = [ n, occ, ha, ev, s ] 


    return df

def write_input_file(orb_i, orb_f, spin, fname, gbwin, gbwout):
    ' write orca_loc input file with initial, final orbital and spin'

    # template with placeholders
    template = '''\
{gbwin:10}  # input gbw
{gbwout:10}  # output gbw
{orb_i:10}  # orbital window: first orbital to be localized 
{orb_f:10}  # orbital window: last orbital to be localized 
2           # localization method 1=PIPEK-MEZEY,2=FOSTER-BOYS,3=IAO-IBO,4=IAO-BOYS, 5=NEW-BOYS 6=FB AUGHESS
{spin:10}  # operator:0 for alpha, 1 for beta
1200        # maximum number of iterations
1e-6        # convergence tolerance of the localization functional value
0.0         # relative convergence tolerance of the localization functional value
100         # printing thresh to call an orbital strongly localized: all should be printed as deloc
100         # printing thresh to call an orbital bond-like: all should be printed as deloc
2           # print level
1           # use Cholesky Decomposition (0=false, 1=true, default is true,optional)
0           # Randomize seed for localization(optional)
'''

    # dict to fill template
    context = {
        'gbwin' : str(gbwin),
        'gbwout': str(gbwout),
        'orb_i' : str(orb_i),
        'orb_f' : str(orb_f),
        'spin'  : str(spin), }

    # write file
    with open(fname, 'w') as f:
        f.write(template.format(**context))


def get_orbital_range(df, spin, minerg):
    '''return minimal and maximal orbital index 

    returns tuple with the number of the highest occupied orbital and lowest with energy of at least minerg

    requires
    df       dataframe created by get_orbitals()
    spin     alpha: 0, beta: 1
    minerg   minimal orbital energy in Ha, default: 1'''
    
    # select only occupied & of certain spin & with at least minerg in ha
    df = df.loc[ ( df.loc[:,'occ'] != 0 ) & ( df.loc[:,'spin'] == spin ) & (df.loc[:,'e_ha'] > minerg)]
    
    # get min and max orbital index
    nmin = df.loc[:,'n'].min()
    nmax = df.loc[:,'n'].max()

    return (nmin, nmax)


def run():

    # setup command line behavior
    parser = argparse.ArgumentParser(description='determine orbital ranges from orca output and write orca_loc input files')
    parser.add_argument('output', nargs='?', help='give output file of orca calulation', default='orca1.out')
    args = parser.parse_args()

    # files
    orcaout      = args.output
    gbw          = re.sub(r'(mpi\d*\.)?out$', 'gbw', orcaout)
    alphainput   = 'orca_loc_a.input'
    betainput    = 'orca_loc_b.input'
    alphagbw     = 'loc_a.gbw'
    alphabetagbw = 'loc.gbw'
    
    # get orbital raneg of occupied orbitals with energy > -1 Ha
    print('Reading {}'.format(orcaout))
    df = get_orbitals(orcaout)
    print('    ... getting alpha orbital range')
    a_min, a_max = get_orbital_range(df, spin=0, minerg=-1)
    print('    ... getting beta orbital range')
    b_min, b_max = get_orbital_range(df, spin=1, minerg=-1)
    
    # write files
    print('Writing files for orca_loc utility')
    print('    ... writing {}'.format(alphainput))
    write_input_file(a_min, a_max, 0, alphainput, gbw, alphagbw)
    print('    ... writing {}'.format(betainput))
    write_input_file(b_min, b_max, 1, betainput, alphagbw, alphabetagbw)

    print('    ... done!')
if __name__ == '__main__':
    run()
