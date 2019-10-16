#!/bin/bash
#
# Aidan Fike
# June 27, 2019
#
# Program to obtain the trajectory of each omega angle of the cyclic peptide. 

structNum=$1
numRes=$2
gromacs=$3

if [ $gromacs -eq 4 ]
then
    editconf_mpi -f s${structNum}_int_noh.pdb -o curr.gro &> edit.log
else
    gmx_mpi editconf -f s${structNum}_int_noh.pdb -o curr.gro &> edit.log
fi

sed -e "s#NUMRES#${numRes}#g" code/VMD_GenOmegaIndex_template.sh > VMD_GenOmegaIndex.sh

vmd -e VMD_GenOmegaIndex.sh &> vmd.log #Compute the index file to indicate the omega dihedrals

if [ $gromacs -eq 4 ]
then
    g_angle_mpi -f curr.gro -n index_omega.ndx -ov struct_omega.xvg -all -type dihedral &> cis.log
else
    gmx_mpi angle -f curr.gro -n index_omega.ndx -ov struct_omega.xvg -all -type dihedral &> cis.log
fi
    

python code/Py_searchForCis.py > cisOut${structNum}.txt
