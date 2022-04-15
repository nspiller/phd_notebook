#!/usr/bin/env python3
import pandas as pd
import numpy as np
from numpy import sin, cos, sqrt 
import re
from scipy.optimize import minimize
import argparse

#############
## USER INPUT

# each element in this list is one ligand: L'-L-M 
# define through position in xyz file
# e.g. [2, 1, 0] indicates H, O, Fe in xyz file
#   Fe 0 0 0
#   O  1 0 0
#   H  1 1 0
connectivity = [
    [5, 1, 0],
    [6, 2, 0],
    [7, 3, 0],
    [8, 4, 0]
]

# define AOM parameters as list of tuples
# (parameter, the starting value, upper constraint, lower constraint)
# valid parameters: σ, π, πs, πc, σπc, σπs
aom_params = [
    ('σ', 5e3, (None, None)),
    ('πs', 1e3, (None, None)),
    ('πc', 0, (None, None)),
    ('σπc', 0, (None, None))
]

# additional constraints passed to scipy's minimize
constraints = ( # constraints defined as functions based on aom_params
             {'type': 'eq', 'fun': lambda x:  x[0] * x[2] - x[3]**2 }, # σ * πx = σπx**2
    )

# use either CASSCF or NEVPT2 AILFT elements
correlation = 'NEVPT2'

############
## FUNCTIONS

def get_xyz(output):
    'get xyz block from orca output file and return numpy array'

    with open(output) as f:
        parse = False
        for line in f:

            if 'CARTESIAN COORDINATES (ANGSTROEM)' in line:
                parse = True
                xyz_lines = []
                next(f) # skip line with dashes
                continue

            if parse:
                if len(line.split()) == 0:
                    break
                else:
                    xyz_lines.append(line.split())

    xyz_coords = [ [i[1], i[2], i[3]] for i in xyz_lines ]
    xyz = np.array(xyz_coords)
    xyz = xyz.astype(float)

    return xyz

def get_normalized(v):
    'return normalized numpy array'
    vn = v / np.linalg.norm(v)
    return vn

def measure_dihedral(xyz, index_list):
    'return dihedral angle that is formed by index_list with z axis'
    i1, i2, i3 = index_list

    vec1 = xyz[int(i1)]
    vec2 = xyz[int(i2)]
    vec3 = xyz[int(i3)]
    vec4 = np.array([0, 0, 1])

    vec21_norm = get_normalized(vec1 - vec2)
    vec32_norm = get_normalized(vec2 - vec3)
    vec34_norm = get_normalized(vec4 - vec3)

    # build orthonormal coordninate system:
    # z axis is along vec32
    # y axis is orthogonal to vec12 and vec23
    # x axis is projection of vec12 in xy plane
    ez = vec32_norm
    ey = np.cross(vec21_norm, ez)
    ex = np.cross(ez, ey)

    vec34_normx = np.dot(ex, vec34_norm)
    vec34_normy = np.dot(ey, vec34_norm)
    
    dihedral = np.arctan2(vec34_normy, vec34_normx)
    
    return dihedral


def get_ligand_angles(xyz, index_list):
    '''returns the tuple of angles θ, ϕ, ψ in radians for a ligand
    
    takes two arguments:
    1: numpy array with xyz data
    2: index list: [ X, L, O ] to describe the connectivity
    L is atom index of the ligand
    X the index of the atom connected to L,
    O is the index of the origin

    θ: angle between orig->L and z-axis
    ϕ: angle between projection of orig->L to xy-plane and x-axis
    ψ: dihedral angle spanned by X->L and z axis 
    '''
    
    X, L, O = index_list # unpack X, L, O

    # unit vector Z
    EZ = np.array([0, 0, 1]) - xyz[O]
    EZ = get_normalized(EZ)

    vOL = xyz[L] - xyz[O] # origin -> L
    vOL = get_normalized(vOL) # normalize

    # θ: angle between orig->L and z-axis
    θ = np.arccos( np.dot(vOL, EZ) )

    # ϕ: angle between projection of orig->L to xy-plane and x-axis
    ϕ = np.arctan2(vOL[1], vOL[0]) # np.arctan2 gives correct sign for polar coordinates

    #ψ: dihedral angle 
    ψ = measure_dihedral(xyz, index_list) 

    # write list of [θ, ϕ, ψ]
    angles = np.array( [ θ, ϕ, ψ ], dtype=float)
    
    return angles

def load_ligand_angles(xyz, con):
    '''
    returns the ligand angles for a given connectivity and xyz file
    takes two arguments
        xyz             xyz file with coordinates, indexing starts with 0
        con_list        list of connectivity lists 
                        [ [ X, L, O ], [ X2, L2, O], … ]
    '''
    lig_ang = []
    for c in con:
        la = get_ligand_angles(xyz, c)
        lig_ang.append(la)

    return lig_ang

def traceless(numpy_array):
    'subtract from diagonal elements to make matrix traceless'

    dim = numpy_array.shape[0]
    uni = np.identity(dim)
    sub = numpy_array.trace() / dim * np.identity(dim)
    numpy_array_traceless = numpy_array - sub
    return numpy_array_traceless

def ailft_matele(outputfile, correlation):
    'returns an np.array with the traceless dorb energy matrix'

    # get parameter block from file
    params_block = []
    with open(outputfile) as f:
        switch_params_block = False
        if correlation == 'CASSCF':
            match_params_block = 'AILFT MATRIX ELEMENTS (CASSCF)'
        elif correlation == 'NEVPT2':
            match_params_block = 'AILFT MATRIX ELEMENTS (NEVPT2)'
        for line in f:
            if match_params_block in line:
                switch_params_block = True
                next(f)
                next(f)
                continue
            if re.match('^\s*$', line):
                switch_params_block = False
                continue
            if switch_params_block == True:
                params_block.append(line)

    # assemble dict for matrix creation
    d = {}
    for line in params_block:
        l = line.split('=')
        l[0] = ''.join(l[0].split())
        l[3] = l[3].strip()
        l[3] = l[3].split()[0]
        d[ l[0] ] = float( l[3] )

    d_ = {} # only keep matrix elements starting with H
    for key in d:
        if re.match('^H', key):
            d_[key] = d[key]
    d = d_

    # build symmetric matrix
    dorb_matrix = np.array([
    [ d['H(dx2-y2,dx2-y2)'], d['H(dx2-y2,dz2)'], d['H(dx2-y2,dxy)'], d['H(dx2-y2,dxz)'], d['H(dx2-y2,dyz)'] ],
    [ d['H(dx2-y2,dz2)']   , d['H(dz2,dz2)'],    d['H(dz2,dxy)'],    d['H(dxz,dz2)'],    d['H(dz2,dyz)'] ],    
    [ d['H(dx2-y2,dxy)']   , d['H(dz2,dxy)'],    d['H(dxy,dxy)'],    d['H(dxz,dxy)'],    d['H(dyz,dxy)'] ],
    [ d['H(dx2-y2,dxz)']   , d['H(dxz,dz2)'],    d['H(dxz,dxy)'],    d['H(dxz,dxz)'],    d['H(dxz,dyz)'] ],
    [ d['H(dx2-y2,dyz)']   , d['H(dz2,dyz)'],    d['H(dyz,dxy)'],    d['H(dxz,dyz)'],    d['H(dyz,dyz)'] ]
    ], dtype='float_')
    dorb_matrix = traceless(dorb_matrix) # make traceless
    
    # vector from matrix
    basis = [
    (0,0),
    (1,0), (1,1), 
    (2,0), (2,1), (2,2),
    (3,0), (3,1), (3,2), (3,3),
    (4,0), (4,1), (4,2), (4,3), (4,4)
    ]
    vec = np.array([dorb_matrix[i] for i in basis], dtype='float_')

    return vec

###############
## AOM factors
def x2y2(angles):
    θ, ϕ, ψ = angles    
    Fσ = sqrt(3) / 4 * cos(2 * ϕ) * (1 - cos(2 * θ))
    Fπs = -sin(2 * ϕ) * sin(θ) * cos(ψ) - 0.5 * cos(2 * ϕ) * sin(2 * θ) * sin(ψ)
    Fπc = -sin(2 * ϕ) * sin(θ) * sin(ψ) + 0.5 * cos(2 * ϕ) * sin(2 * θ) * cos(ψ)
    return Fσ, Fπs, Fπc

def z2(angles):
    θ, ϕ, ψ = angles     
    Fσ = (1 + 3 * cos(2 * θ)) / 4
    Fπs = sqrt(3) / 2 * sin(2 * θ) * sin(ψ)
    Fπc = -sqrt(3) / 2 * sin(2 * θ) * cos(ψ)
    return Fσ, Fπs, Fπc

def xy(angles):
    θ, ϕ, ψ = angles     
    Fσ = sqrt(3) / 4 * sin(2 * ϕ) * (1 - cos(2 * θ))
    Fπs = cos(2 * ϕ) * sin(θ) * cos(ψ) - 0.5 * sin(2 * ϕ) * sin(2 * θ) * sin(ψ)
    Fπc = cos(2 * ϕ) * sin(θ) * sin(ψ) + 0.5 * sin(2 * ϕ) * sin(2 * θ) * cos(ψ)
    return Fσ, Fπs, Fπc

def xz(angles):
    θ, ϕ, ψ = angles     
    Fσ = sqrt(3) / 2 * cos(ϕ) * sin(2 * θ)
    Fπs = -sin(ϕ) * cos(θ) * cos(ψ) - cos(ϕ) * cos(2 * θ) * sin(ψ)
    Fπc = -sin(ϕ) * cos(θ) * sin(ψ) + cos(ϕ) * cos(2 * θ) * cos(ψ)
    return Fσ, Fπs, Fπc

def yz(angles):
    θ, ϕ, ψ = angles    
    Fσ = sqrt(3) / 2 * sin(ϕ) * sin(2 * θ) 
    Fπs = cos(ϕ) * cos(θ) * cos(ψ) - sin(ϕ) * cos(2 * θ) * sin(ψ)
    Fπc = cos(ϕ) * cos(θ) * sin(ψ) + sin(ϕ) * cos(2 * θ) * cos(ψ)
    return Fσ, Fπs, Fπc

def aom_vec(angles, param_type):
    '''
    takes np.array with three angles for ligand: θ, ϕ, ψ
    returns np.array with shape (15,1) with aom factors according to param_type:
    σ
    πs
    πc
    
    same structure as ailft vector (also "traceless")
    '''

    # assing aom factors
    x2y2_σ, x2y2_πs, x2y2_πc = x2y2(angles)
    z2_σ, z2_πs, z2_πc = z2(angles)
    xy_σ, xy_πs, xy_πc = xy(angles)
    xz_σ, xz_πs, xz_πc = xz(angles)
    yz_σ, yz_πs, yz_πc = yz(angles)

    if param_type == 'σ': # σ parameters
        vec = np.array([
            x2y2_σ * x2y2_σ,
            z2_σ * x2y2_σ, z2_σ * z2_σ,
            xy_σ * x2y2_σ, xy_σ * z2_σ, xy_σ * xy_σ,
            xz_σ * x2y2_σ, xz_σ * z2_σ, xz_σ * xy_σ, xz_σ * xz_σ,
            yz_σ * x2y2_σ, yz_σ * z2_σ, yz_σ * xy_σ, yz_σ * xz_σ, yz_σ * yz_σ
        ])
        
    elif param_type == 'πs':
        vec = np.array([
            x2y2_πs * x2y2_πs,
            z2_πs * x2y2_πs, z2_πs * z2_πs,
            xy_πs * x2y2_πs, xy_πs * z2_πs, xy_πs * xy_πs,
            xz_πs * x2y2_πs, xz_πs * z2_πs, xz_πs * xy_πs, xz_πs * xz_πs,
            yz_πs * x2y2_πs, yz_πs * z2_πs, yz_πs * xy_πs, yz_πs * xz_πs, yz_πs * yz_πs
        ])
        
    elif param_type == 'πc':
        vec = np.array([
            x2y2_πc * x2y2_πc,
            z2_πc * x2y2_πc, z2_πc * z2_πc,
            xy_πc * x2y2_πc, xy_πc * z2_πc, xy_πc * xy_πc,
            xz_πc * x2y2_πc, xz_πc * z2_πc, xz_πc * xy_πc, xz_πc * xz_πc,
            yz_πc * x2y2_πc, yz_πc * z2_πc, yz_πc * xy_πc, yz_πc * xz_πc, yz_πc * yz_πc
        ])
        
    elif param_type == 'π':
        vec_πc = np.array([
            x2y2_πc * x2y2_πc,
            z2_πc * x2y2_πc, z2_πc * z2_πc,
            xy_πc * x2y2_πc, xy_πc * z2_πc, xy_πc * xy_πc,
            xz_πc * x2y2_πc, xz_πc * z2_πc, xz_πc * xy_πc, xz_πc * xz_πc,
            yz_πc * x2y2_πc, yz_πc * z2_πc, yz_πc * xy_πc, yz_πc * xz_πc, yz_πc * yz_πc
        ])

        vec_πs = np.array([
            x2y2_πs * x2y2_πs,
            z2_πs * x2y2_πs, z2_πs * z2_πs,
            xy_πs * x2y2_πs, xy_πs * z2_πs, xy_πs * xy_πs,
            xz_πs * x2y2_πs, xz_πs * z2_πs, xz_πs * xy_πs, xz_πs * xz_πs,
            yz_πs * x2y2_πs, yz_πs * z2_πs, yz_πs * xy_πs, yz_πs * xz_πs, yz_πs * yz_πs
        ])

        vec = vec_πs + vec_πc

    elif param_type == 'σπs':
        
        # set ψ = 0, because ψ has been taken into account in πs
        angles_ = np.array([ angles[0], angles[1], 0 ]) 
        
        # assing factors again
        x2y2_σ, x2y2_πs, x2y2_πc = x2y2(angles)
        z2_σ, z2_πs, z2_πc = z2(angles)
        xy_σ, xy_πs, xy_πc = xy(angles)
        xz_σ, xz_πs, xz_πc = xz(angles)
        yz_σ, yz_πs, yz_πc = yz(angles)
        
        vec = np.array([ 
            x2y2_σ * x2y2_πs * 4,
            0, z2_σ * z2_πs * 4,
            0, 0, xy_σ * xy_πs * 4,
            0, 0, 0, xz_σ * xz_πs * 4,
            0, 0, 0, 0, yz_σ * yz_πs * 4
        ])
        
    elif param_type == 'σπc':
        
        # set ψ = 0, because ψ has been taken into account in πc
        angles_ = np.array([ angles[0], angles[1], 0 ]) 
        
        # assing factors again
        x2y2_σ, x2y2_πs, x2y2_πc = x2y2(angles)
        z2_σ, z2_πs, z2_πc = z2(angles)
        xy_σ, xy_πs, xy_πc = xy(angles)
        xz_σ, xz_πs, xz_πc = xz(angles)
        yz_σ, yz_πs, yz_πc = yz(angles)
        
        vec = np.array([ 
            x2y2_σ * x2y2_πc * 4,
            0, z2_σ * z2_πc * 4,
            0, 0, xy_σ * xy_πc * 4,
            0, 0, 0, xz_σ * xz_πc * 4,
            0, 0, 0, 0, yz_σ * yz_πc * 4
        ])

    vec = traceless_aom_param(vec) # make traceless as with AILFT vec

    vec = vec.reshape( (15,1) ) # to make matrix multiplication possible

    return vec

def aom_matrix(lig_ang, params):
    '''
    create matrix containing the AOM factors in the order of ailft_matele (see OrcaAILFT)
    arguments
        lig_ang     ligand angles (list of np.array([θ, ϕ, ψ))
        params      list of parameters 
                    possible: ['σ', 'π', 'πs', 'πc', 'σπs', 'σπc']

    returns hstack of vectors as matrix in order given by params
    '''

    param_dict = { p: np.zeros( (15,1) ) for p in params }
    
    for ang in lig_ang:
        for p in param_dict:
            param_dict[p] += aom_vec(ang, param_type=p)

    M = np.hstack( [ param_dict[p] for p in param_dict ] )

    return M


def traceless_aom_param(vec, diag=[0,2,5,9,14]):
    '''
    make vec with shape (15,) containing the AOM factors so that it 
    corresponds to a traceless matrix select the positions of the diagonal 
    elements according to diag (default for ailft calculation)
    returns vec with shape (15,)
    '''
    
    trace = sum(vec[diag]) # determine trace
    sub = trace / len(diag) # devide by number of diagonal elements
    vec[diag] = vec[diag] - sub # subtract from diagonal elements
    
    return vec

def matrify(a):
    '''make 5x5 matrix out of array shaped like:
        x2y2 * x2y2,
        z2 * x2y2, z2 * z2,
        xy * x2y2, xy * z2, xy * xy,
        xz * x2y2, xz * z2, xz * xy, xz * xz,
        yz * x2y2, yz * z2, yz * xy, yz * xz, yz * yz
    '''
    a = a.reshape( (15,) )
    m = pd.DataFrame([
        [ a[0],  a[1],  a[3],  a[6],  a[10] ],
        [ a[1],  a[2],  a[4],  a[7],  a[11] ],
        [ a[3],  a[4],  a[5],  a[8],  a[12] ],
        [ a[6],  a[7],  a[8],  a[9],  a[13] ],
        [ a[10], a[11], a[12], a[13], a[14] ],
    ])
    return m

def f_rms(x, M, a):
    'return RMS for the linear equations given by M*x = a'
    val = np.sqrt( ( (np.dot(M, x) - a)**2 ).mean() )
    return val

#############
## RUN SCRIPT

def run(out, con, par, cstr, cor):

    # get angles
    xyz = get_xyz(out)
    ang = load_ligand_angles(xyz, con)
    print('Ligand angles (in rad):')
    for a in ang:
        print('θ: {: .4f},  ϕ: {: .4f}, ψ:{: .4f}'.format(a[0], a[1], a[2]))
    print()

    # get AILFT matrix
    ailft_vec = ailft_matele(out, cor)
    print('AILFT maitrx ({}, traceless, in cm-1, < 10 set to 0):'.format(cor))
    ailft_vec[np.abs(ailft_vec) < 10] = 0
    print(matrify(ailft_vec))
    print()

    # construct AOM matrix
    aom_mat = aom_matrix(ang, [p[0] for p in par])
    for i, p in enumerate(par):
        print('AOM matrix for {} (traceless, < 1e-10 set to 0):'.format(p[0]))
        m = aom_mat[:,i]
        m[np.abs(m) < 1e-10] = 0
        print(matrify(m))
        print()

    # fit 
    par0 = [ p[1] for p in aom_params ]
    bnds = [ p[2] for p in aom_params ]

    min_res = minimize(
        f_rms, par0, method='SLSQP', jac=False, args=(aom_mat, ailft_vec), 
        bounds=bnds, constraints=cstr, 
        options={'maxiter': 10000, 'ftol': 1e-8, 'eps': 1e-8, })


    # results
    print('Results of the AOM fit (in cm-1)' )
    fmt = '    {:4s}:  {:1.2f}'
    for i, p in enumerate(par):
        print(fmt.format(p[0], min_res.x[i]))
    print(fmt.format('RMS', min_res.fun))

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('orca_output', help='Name of the output file containing AILFT calculation')
    args = parser.parse_args()
    
    run(args.orca_output, connectivity, aom_params, constraints, correlation)
