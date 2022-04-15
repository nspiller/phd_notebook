#!/usr/bin/env python3
import numpy as np
from math import factorial
from fractions import Fraction
from pandas import DataFrame
import sys


def bincoef(n, k):
    """
    return the binomial coefficient for two numbers:
    n over k = n! / (k! (n-k)!)
    """
    bc = factorial(n) / (factorial(k) * factorial(n - k))
    return bc

def calculate_states(e, o):
    '''
    Return two list of lists for distributing  
    e electrons among
    o spatial orbitals
    
    first list: possible Ms states >= 0

    second list: number of states for each spin 

    output read by root_table()
    ''' 
    if e > 2 * o or e < 0:
        raise ValueError('Number of electrons must be bewteen 0 and 2 * o')

    # α and β for high spin state
    nα_hs = o if e > o else e
    nβ_hs = e - o if e > o else 0
    
    ms_data = []
    s_data = []
    for i in range(nα_hs): # flip α-electrons one at a time
        nα = nα_hs - i
        nβ = nβ_hs + i 
        Ms = (nα - nβ) / 2
        if Ms < 0:
            break

        ms_states = int(bincoef(o, nα) * bincoef(o, nβ)) # possible states with ms
        ms_data.append([Ms, nα, nβ, ms_states])

        s_roots = int((1 -  (o - nα) * nβ / ( (nα + 1) * (o - nβ + 1) ) ) * bincoef(o, nα) * bincoef(o, nβ)) # possible states with s
        mult = int(2 * Ms + 1)
        s_data.append([mult, s_roots])

        
    return [ms_data, s_data]

    
def root_table(data):
    '''
    Print a table for list of lists:
    '''
  
    ms_data, s_data = data
    table = [[
        'Ms>=0',
        'nα',
        'nβ',
        'microstates',
        '|',
        'S',
        'mult',
        'states/roots'
        ]]

    total_ms_states = 0
    total_s_states = 0
    for d_ms, d_s in zip(ms_data, s_data):
        ms, nα, nβ, ms_states = d_ms 
        total_ms_states += ms_states

        mult, s_states = d_s
        s = (mult - 1) / 2
        total_s_states += s_states

        table.append([str(Fraction(ms)), nα, nβ, ms_states, '|', str(Fraction(s)), mult, s_states])

    table.append(['total',  '', '', str(total_ms_states), '|', 'total', '', str(total_s_states)])


    return table 

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'e', 
        help='number of electrons', 
        type=int)
    parser.add_argument(
        'o',
        help='number of orbitals',
        type=int)
    args = parser.parse_args()

    df = DataFrame( root_table(calculate_states(args.e,args.o)) ) 
    print(df.to_string(index=False, header=False))

