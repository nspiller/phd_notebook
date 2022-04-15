#!/bin/bash

# carry out multiwfn analysis of orca output files
#
# User needs to set the path for orca, multiwfn and the parser folder
# expects orca1.out and orca1.gbw to be present 


#############################
## Environmental Variables ##
ORCAPATH= # path to orca
MPILIB= # oath to openmpi/lib folder (only for shared orca)
MPIPATH= # path to openmpi/bin folder (only for shared orca)

MULTIWFNPATH= # path to multiwfn executable

PARSERPATH= # path to multiWFN_analysis folder
MULTI_INP=$PARSERPATH/multi_hirsh.inp 
WRITE_LOC=$PARSERPATH/write_orca_loc_input.py
PARSE_MULTI=$PARSERPATH/parse_multiwfn.py

##########
## Orca ##
export PATH="$MPIPATH:$ORCAPATH:$PATH"
export LD_LIBRARY_PATH="$MPILIB:$ORCAPATH:$LD_LIBRARY_PATH"
ORCA_LOC=$(which orca_loc)
ORCA_2MKL=$(which orca_2mkl)

$WRITE_LOC

$ORCA_LOC orca_loc_a.input && $ORCA_LOC orca_loc_b.input 
$ORCA_2MKL loc -molden

##############
## MultiWFN ##
ulimit -u 8191                                                                  
ulimit -c 0                                                                     
ulimit -s unlimited                                                             
export KMP_STACKSIZE=64000000                                                   
export PATH="$MULTIWFNPATH:$PATH"
MULTI=$(which Multiwfn)
MULTI_OUT=$( basename $MULTI_INP inp)out

$MULTI loc.molden.input < $MULTI_INP | tee $MULTI_OUT

#############
## Parsing ##
$PARSE_MULTI
