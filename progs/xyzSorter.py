#!/usr/bin/env python3
import numpy as np
from copy import deepcopy
import sys

class xyz_file:
    '''
    class to read, manipulate and write xyz files
    
    create with path to xyz file

    file must not contain empty lines
    '''
    def __init__(self, infile):
        self.infile = infile
        self.__read_xyz()

    def __read_xyz(self):
        '''from xyz extract:
            - 2 header lines: self.header
            - list of np.array([x, y ,z]) for each atom: xyz_list
            - list of 'atom1', 'atom2', … : atom_list
        '''
        with open(self.infile) as f:
        
            # prase first two lines as headers
            header1 = f.readline().rstrip() 
            header2 = f.readline().rstrip()
        
            # create and fill xyz and atom list
            xyz_list, atom_list = [], [] 
            for line in f:

                ls = line.split()  
                atom, x, y, z = ls[:5]

                xyz_list.append(np.array([x, y, z], dtype='float_'))
                atom_list.append(str(atom))

        self.header1 = header1
        self.header2 = header2
        self.xyz_list = xyz_list
        self.atom_list = atom_list

    def switch(self, i, j):
        '''switch the positions (lines) of atom i and j (numbering starts with 0)'''
        self.xyz_list[j], self.xyz_list[i] = self.xyz_list[i], self.xyz_list[j]
        self.atom_list[j], self.atom_list[i] = self.atom_list[i], self.atom_list[j]


    def __make_xyz(self):
        '''
        make list containing the xyz file based on current:
            - self.header
            - self.xyz_list
            - self.atom_list
        '''
        contents = []
        contents.append(self.header1 + '\n')
        contents.append(self.header2 + '\n')
        for i, j in zip(self.atom_list, self.xyz_list):
            atom = i
            x, y, z = j[0], j[1], j[2]
            line = '{:4} {: 15.10f} {: 15.10f} {: 15.10f}{}'.format(atom, x, y, z, '\n')
            contents.append(line)
        
        return contents
        
    def write_xyz(self, outfile):
        '''
        write outfile as xyz based on current:
            - self.header
            - self.xyz_list
            - self.atom_list
        '''
        self.outfile = outfile
        contents = self.__make_xyz()
        with open(self.outfile, 'w+') as f:
            for line in contents:
                f.write(line)

    def reorder(self, template):
        '''
        reorder indices by minimizing (individual) distancen to atoms
        in template

        requires path to template xyz file

        very inefficient implementation, only works for small molecules and small displacements
        '''
        # xyz structure to use as template
        tem = xyz_file(template)

        # list to contain the new order
        l = []

        # cycle through each coordinate of templeta and compare distance
        for atom_tem, xyz_tem in zip(tem.atom_list, tem.xyz_list):
            d_comp = 10000 # start with large distance
            for i, (atom_self, xyz_self) in enumerate(zip(self.atom_list, self.xyz_list)):
                if atom_self == atom_tem: # only compare same atom type
                    d = np.linalg.norm(xyz_tem - xyz_self)
                    if d < d_comp: # compare distance
                        d_comp = d # declare new shortest distance and corresponding index
                        li = i
            l.append(li) # append that index with shortest distance

        if len(l) != len(set(l)): # check if l contains only unique indices
            print('ERROR: unambiguous assignment not possible')
            exit() 

        # reorder xyz and atomlist according to l
        self_ = deepcopy(self)
        for i, j in enumerate(l):
            self.xyz_list[i] = self_.xyz_list[j]
            self.atom_list[i] = self_.atom_list[j]


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('xyz', help='xyz file that needs sorting',)
    parser.add_argument('xyz_template', help='xyz file that serves as template')
    parser.add_argument('xyz_output', help='new file to be written')

    args = parser.parse_args()

    inp = args.xyz
    tem = args.xyz_template
    out = args.xyz_output

    geom = xyz_file(inp)
    geom.reorder(tem)
    geom.write_xyz(out)
