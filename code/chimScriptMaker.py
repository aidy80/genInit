#Aidan Fike
#July 21, 2019
#
#File to create a chimera script which creates a random pdb structure

import os
import numpy as np

#File to create a chimera script to create a structure with amino acids in "seq", and with the
#passed phi/psi angles
def createChimScript(seq, phi, psi, structNum):
    scriptName = "chimScript.py"
    chim = open(scriptName, "w") 
    chim.write("import chimera\n")
    chim.write("from chimera import runCommand\n\n")

    chim.write("runCommand('open code/gly.pdb')\n\n")

    chim.write("runCommand('rotation 1 :1@C :1@CA')\n")
    chim.write("runCommand('rotation 1 %s')\n" % phi[0])
    chim.write("runCommand('~select :1@C :1@CA')\n")

    for index in range(1,len(seq)):
        chim.write("runCommand('addaa gly,%s,%s,%s %s')\n" % (str(index+1), phi[index],\
                                                                psi[index],":"+str(index)+".a"))

    chim.write("runCommand('delete :%s@OXT')\n" % len(seq))
    chim.write("runCommand('bond :%s@C :1@N')\n" % len(seq))

    chim.write("runCommand('minimize nogui True, nsteps 1000')\n")

    for index, amino in enumerate(seq):
        chim.write("runCommand('swapaa %s #:%s.a')\n" % (oneToThree(amino), index+1))

    for i in range(len(seq)):
        chim.write("runCommand('chirality :%d.A@CA')\n" % (i+1)) 

    chim.write("runCommand('minimize nogui True, nsteps 1000')\n")

    chim.write("runCommand('select element.H')\n")
    chim.write("runCommand('~select :1.A@H')\n")
    chim.write("runCommand('delete sel')\n")

    chim.write("runCommand('write #0 temp%d.pdb ')\n" % structNum) 

    chim.close()

    return scriptName

def checkForDAminos(structNum, newSeq):
    chimLog = open("s%d.chimOut" % structNum, "r")

    chirality = []
    
    for line in chimLog:
        words = line.split()
        for i in range(len(newSeq)):
            if len(words) > 0 and words[0] == "#0:%d.A@CA" % (i+1):
                chirality.append(words[2])
    chimLog.close()

    scriptName = "dAminoScript.py"
    damino = open(scriptName, "w")

    damino.write("import chimera\n")
    damino.write("from chimera import runCommand\n\n")

    damino.write("runCommand('open temp%d.pdb')\n" % structNum) 

    dFound = False
    for i in range(len(chirality)):
        if (chirality[i] == 'R' and newSeq[i] != 'C') \
                    or (chirality[i] == 'S' and newSeq[i] == 'C'):
            print "D amino acid found at %d%s in struct %d. Now inverting it to an L amino acid" % (i+1,newSeq[i],structNum)
            damino.write("runCommand('invert :%d.A@CA')\n" % (i+1))
            dFound = True

    if dFound:
        damino.write("runCommand('write #0 temp%d.pdb ')\n" % structNum)
    #    damino.write("runCommand('minimize nogui True, nsteps 1000')\n")

#        damino.write("runCommand('select element.H')\n")
#        damino.write("runCommand('~select :1.A@H')\n")
#        damino.write("runCommand('delete sel')\n")
    else:
        print "No D amino acids found in struct %d" % structNum

    damino.close()

    if dFound:
        os.system("chimera  --script " + scriptName + " --nogui &> s%d.dchimOut" % structNum)

    os.remove(scriptName)


#Insert into a pdb file the dihedral connections which are missed by chimera because of the 
#weird cyclic nature
def insertMissingDihed(structNum, numAminos):
    oldPdb = open("temp%d.pdb" % structNum, "r")
    newPdb = open("s%s_int_noh.pdb" % structNum, "w")

    N1 = -1
    CA1 = -1 
    H1 = -1
    CAx = -1
    Cx = -1
    Ox = -1

    #Find relevant atoms in the current pdb and write the new pdb
    for line in oldPdb:
        words = line.split()
        if words[0] != "END":
            if words[5] == "1":
                if words[2] == "N":
                    N1 = int(words[1])
                if words[2] == "CA":
                    CA1 = int(words[1])
                if words[2] == "H":
                    H1 = int(words[1])
            elif words[5] == str(numAminos):
                if words[2] == "C":
                    Cx = int(words[1])
                if words[2] == "CA":
                    CAx = int(words[1])
                if words[2] == "O":
                    Ox = int(words[1])

            newPdb.write(line) 
        else:
            newPdb.write("CONECT %4d %4d %4d %4d\n" % (N1, CA1, Cx, H1))
            newPdb.write("CONECT %4d %4d %4d %4d\n" % (Cx, N1, CAx, Ox))
            newPdb.write("END\n")

    oldPdb.close()
    newPdb.close()

    os.remove("temp%d.pdb" % structNum)
    print "Successfully created s%d_int_noh.pdb" % structNum
    
#Code to convert one letter amino acid codes to three-letter amino acid codes
def oneToThree(one):
    if one == "A":
        return "ala"
    if one == "C":
        return "cys"
    if one == "D":
        return "asp"
    if one == "E":
        return "glu"
    if one == "F":
        return "phe"
    if one == "G":
        return "gly"
    if one == "H":
        return "his"
    if one == "I":
        return "ile"
    if one == "K":
        return "lys"
    if one == "L":
        return "leu"
    if one == "M":
        return "met"
    if one == "N":
        return "asn"
    if one == "P":
        return "pro"
    if one == "Q":
        return "gln"
    if one == "R":
        return "arg"
    if one == "S":
        return "ser"
    if one == "T":
        return "thr"
    if one == "V":
        return "val"
    if one == "W":
        return "trp"
    if one == "Y":
        return "tyr"
