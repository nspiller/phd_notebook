#!/usr/bin/env python3

import numpy as np
from scipy.stats import linregress
import glob
import re
import matplotlib.pylab as plt

compund_dict = {
    '01': ( 0.9  , 'fe2cl4'),
    '02': ( -0.02, 'fe2cn6'),
    '03': ( 1.34 , 'fe2f6'),
    '04': ( 1.05 , 'fe2poroac'),
    '05': ( 0.56 , 'fe2sr3'),
    '06': ( 0.29 , 'fe3az'),
    '07': ( 0.19 , 'fe3cl4'),
#    '08': ( ?    , 'fe3cl6'),
    '09': ( -0.13, 'fe3cn6'),
    '10': ( 0.48 , 'fe3f6'),
    '11': ( 0.51 , 'fe3h2o6'),
    '12': ( 0.15 , 'fe3mac'),
    '13': ( 0.2  , 'fe3oeppy'),
    '14': ( 0.67 , 'fe3poro2'),
    '15': ( 0.4  , 'fe3poroac'),
    '16': ( -0.02, 'fe4mac'),
    '17': ( 0.08 , 'fe4poro'),
    '18': ( 0.17 , 'fe4tmco'),
    '19': ( -0.87, 'fe6o4'),
    '20': ( 0    , 'feco5'),
    '21': ( 0.04 , 'feno6'),
    '22': ( 0.33 , 'feno7'),
    '23': ( 0.34 , 'feph3'),
    '24': ( 0.44 , 'fesme'),
}


def get_moess(outputfile):
    '''
    return rho 0 from output file for single iron center

    matches 
    -----------------------------------------
    ELECTRIC AND MAGNETIC HYPERFINE STRUCTURE
    -----------------------------------------
    […]
     Delta-EQ=(1/2{e**2qQ}*sqrt(1+1/3*eta**2) =   -7.887 MHz =   -0.680 mm/s
     RHO(0)=  14923.321525954 a.u.**-3
     '''

    match_block = re.compile('ELECTRIC AND MAGNETIC HYPERFINE STRUCTURE')
    match_rho = re.compile('RHO\(0\)=\s*(?P<rho>\d+\.\d+) a.u.\*\*-3')

    with open(outputfile) as f:
        switch_block = False
        for line in f:
            if match_block.search(line):
                switch_block = True
            if switch_block == True:
                if match_rho.search(line):
                    rho = match_rho.search(line)
                    rho = rho.groupdict()['rho']
                    rho = np.float(rho)

        return rho

def terminated_normally(output):
    '''check for termition string: ****ORCA TERMINATED NORMALLY****'''
    with open(output) as f:
        for line in f:
            if 'ORCA TERMINATED NORMALLY' in line:
                return True
        return False

rho, delta = [], []
for i in compund_dict:
    output = glob.glob('cal/' + i + '/*.out')[0] 
    if terminated_normally(output):
        r = get_moess(output)
        d = compund_dict[i][0]
        rho.append(r)
        delta.append(d)

rho, delta = np.array(rho), np.array(delta)

# fit
a, b, r_value, p_value, std_err = linregress(rho, delta)

# write to file
with open('./fit_values.dat', 'w') as f:
    f.write('δ   = a * ρ0 + b\n')
    f.write('a   = {:7.10f}\n'.format(a))
    f.write('b   = {:7.10f}\n'.format(b))
    f.write('R   = {:7.10f}\n'.format(r_value))
    f.write('R^2 = {:7.10f}\n'.format(r_value**2))
    f.write('p   = {:7.10f}\n'.format(p_value))
    f.write('σ   = {:7.10f}\n'.format(std_err))
    f.write('(according to python package scipy.stats.linregress)')

# create figure
fig, ax = plt.subplots(figsize=(10, 10))

# scatter plot
ax.scatter(rho, delta, marker='x') 

# linear fit 
xmin, xmax = ax.get_xlim()
x = np.linspace(xmin, xmax, 10) 
fitcurve = lambda x, a, b : a * x + b
ax.plot(x, fitcurve(x, a, b))

# fitted values
ax.text(.9, .9, 
        r'$\delta = \alpha * ( \rho_0 - C ) + \beta  = ax + b $' + '\n' + \
        r'$ a = {:1.5f}$'.format(a) + '\n' + \
        r'$ b = {:1.5f}$'.format(b) + '\n' + \
        r'$ R^2 = {:1.3}$'.format(r_value**2), 
        ha='right', va='top', transform=ax.transAxes)

# formatting
ax.set_xlabel(r'DFT: $\rho (0)$')
ax.set_ylabel(r'exp: $\delta [mm/s]$')
# ax.ticklabel_format(style='sci', axis='x', scilimits=(4,4))
fig.tight_layout()

# save
fig.savefig('./fit_plot.pdf')
