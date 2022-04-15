#!/usr/bin/env python3

def updateIndices(f1, f2, n, x, i0=1):
    '''
    read in file f1 and output to file f2 atom indices above n adjusted by x
    f1          read file
    f2          written file
    n           index at which the change occurs
    x           number of atoms removed or added after n
    i0          default 1 (chemshell counting, vmd starts with 0)
    '''

    with open(f1, 'r') as fin:
        with open(f2, 'w') as fout: 
            for line in fin:
                print('Modifying indices in {}. Old indices in {}.'.format(f2, f1),
                    'Counting indices starts with {}.'.format(i0)) 
                for c in '}{':
                    line = line.replace(c, '')
                l = line.split()
                name, ix = l[0:2], l[2:]
                
                qm = True if name[1] == 'qmatoms' else False
                
                old_ix = [ int(i) for i in ix ]
                chn_ix = [ i + n for i in range(1, x+2) ]

                if x > 0:
                    upd_old_ix = [ i + x if i > n else i for i in old_ix ]
                    new_ix = upd_old_ix if qm == True else upd_old_ix + chn_ix
                
                elif x < 0: 
                    rmv_old_ix = [ i for i in old_ix if i not in chn_ix ]
                    upd_old_ix = [ i + x if i > n else i for i in rmv_old_ix ]
                    new_ix = upd_old_ix 

                name = [str(i) for i in name]
                new_ix = [str(i) for i in new_ix]
                line_new = ' '.join(name) + ' { ' + ' '.join(new_ix) + ' } ' + '\n'
                
                fout.write(line_new)

if __name__ == "__main__":
    import argparse, os
    
    # parse arguments
    parser = argparse.ArgumentParser(description='Update indices in a file containing one line like "set atoms {i1 i2 i3 …}')
    parser.add_argument('region', metavar='region-file',
        help='file containing one tcl list of atom indices')
    parser.add_argument('N', type=int,
        help='atom index after which change occured (see index counting below)')
    parser.add_argument('X', type=int,
        help='number of atoms added (> 0) or deleted (< 0) after N')
    parser.add_argument('-b', metavar='BACKUPFILE', 
        help='Name of backup file. Default: old-*')
    parser.add_argument('-c', metavar='I0', type=int, default=1,
        help='''Number with which index counting starts.
        Chemshell starts with 1 (default). VMD starts with 0''')
    args = parser.parse_args()

    # define variables
    region = args.region
    bak = args.b if args.b else 'old-' + region
    x = args.X
    n = args.N
    i0 = args.c

    # move chm to backup file
    os.rename(region, bak)

    # read in backup file and write changes to new file
    updateIndices(bak, region, n, x, i0=1)

