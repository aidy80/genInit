#!/usr/bin/env bash
# the next line restart using vmd \
exec vmd "-dispdev text -e $0"

################################################################################
#                                                                              #
# Generate the index file for calculating the Omega angles of the residue      #
# by using Gromacs g_dihedral programs.                                        #
#                                                                              #
# This script is designed for cyclic peptide, since the g_rama fails to work   #
# for the calculation of N-terminal Phi angle and C-terminal Psi angle.        #
#                                                                              #
################################################################################


proc getAtomID {molid resid name} {
  set sel [atomselect $molid "resid $resid and name $name"]
  set atomid [$sel get serial]
  return $atomid
}


proc genPhiPsiAtomList {molid resid0 resid1} {
  set atomidList [list]

  lappend atomidList [getAtomID $molid $resid0 "CA"]
  lappend atomidList [getAtomID $molid $resid0 "C"]
  lappend atomidList [getAtomID $molid $resid1 "N"]
  lappend atomidList [getAtomID $molid $resid1 "CA"]
  return $atomidList
}


proc writeAtomList {ofp atomidList} {
  lassign $atomidList a0 a1 a2 a3
  puts $ofp "$a0 $a1 $a2 $a3"
}



#############################################
set numRes NUMRES
set GRO "curr.gro"

set molid [ mol new $GRO type gro ]

# Set up residue list based on CA atoms
set calpha [atomselect $molid "name CA"]
set resids [$calpha get resid]


set ofp [open "index_omega.ndx" w]
puts $ofp {[ Omega ]}

foreach resi $resids {
  if { $resi == 1} {
    set atomidList [ genPhiPsiAtomList $molid 1 2 ]
  } elseif { $resi == $numRes} {
    set atomidList [ genPhiPsiAtomList $molid $numRes 1 ]
  } else {
    set atomidList [ genPhiPsiAtomList $molid $resi [expr $resi+1]]
  }

  writeAtomList $ofp $atomidList

}

close $ofp

exit 0
