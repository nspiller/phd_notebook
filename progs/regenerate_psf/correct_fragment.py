#!/usr/bin/env python3

def extract_elements_pdb(pdb):
    '''
    extract element names as appearing in pdb file
    return list of strings
    '''
    
    el = []
    with open(pdb) as f:
        for line in f:
            if line.startswith('ATOM'):
                el.append(line[76:78].strip())
    return el


def change_elements_fragment(frag1, frag2, el_list):
    '''
    change elements in frag1 to el_list and write to frag2
    '''

    with open(frag1, 'r') as fin:
        with open(frag2, 'w') as fout:
            i_i, i_f = 3, len(el_list) + 3
            for i, line in enumerate(fin):
                i_ = i - 4
                if  4 <= i <= len(el_list) + 3:
                    l = line.split()
                    l[0] = '{:2s}'.format(el_list[ i_ ])
                    line_out = ' '.join(l) + '\n'
                else:
                    line_out = line
                fout.write(line_out)



if __name__ == '__main__':
    import argparse, os

    parser = argparse.ArgumentParser(description='Change elements in fragment file definitions from a parent PDB')
    parser.add_argument('frag', metavar='fragment-file',
        help='Fragment file, e.g. system.c generated from a PDB file from chemshell')
    parser.add_argument('pdb', metavar='pdb-file',
        help='PDB file with proper element names according to PDB format')
    parser.add_argument('-b', metavar='BACKUPFILE',
        help='Name of the backup file. Default: old-*.c')
    args = parser.parse_args()

    frag = args.frag
    pdb = args.pdb
    bak = args.b if args.b else 'old-' + frag

    el_list = extract_elements_pdb(pdb) # extract element list from pdb

    os.rename(frag, bak) # move frag to backup file
    change_elements_fragment(bak, frag, el_list) # read in backup and write frag

