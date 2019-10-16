#Aidan Fike
#July 21, 2019

#Main function. Able to create initial sequences for a given sequence passed by the user.
#See README for more information

import sys

sys.path.insert(1, 'code/')

import chooseStruct 
import parseCommandLine 
import dihedDir 
import chimScriptMaker 
import os
import shutil

maxTries = 20 #The maximum number of structures that are created and evaluated before the program
              #exits (to prevent infinite looping)

def main():
    #Parse command line arguments/obtain the sequence whose structures will be made
    newSeq, thresh, gro4 = parseCommandLine.getArgvInfo()

    #Mk directories
    chooseStruct.mkNoCisDiffs()

    foundPair = False

    #For numberInitStructs, create a structure with random phi psi angles, then only keep it if is
    #does not have any cis peptide bonds
    structNum = 0
    while not foundPair:
        structNum += 1
        
        if structNum > maxTries:
            print "%d structures have been created and evaluated for rmsd/cis bonds but no pair of\n\
structures with a backbone rmsd greater than %fnm were found. Try rerunning the program \n\
with a lower backbone rmsd threshold value for better success (Remember, the rmsd should \n\
be given in [nm], not angstroms).\n" % (maxTries, thresh)
            os.system("rm *.log") 
            os.system("rm *.out") 
            os.system("rm *.chimOut") 
            os.system("rm *.di") 
            shutil.rmtree("noCis")
            shutil.rmtree("diffs")
            sys.exit()

        phi, psi = dihedDir.createRandomPhiPsi(1, len(newSeq))
        dihedDir.writeToDihedFile(newSeq, phi[0], psi[0], structNum)
        scriptName = chimScriptMaker.createChimScript(newSeq, phi[0], psi[0], structNum)
        os.system("chimera  --script " + scriptName + " --nogui &> s%d.chimOut" % structNum)
        os.remove(scriptName)
        chimScriptMaker.checkForDAminos(structNum, newSeq)
        chimScriptMaker.insertMissingDihed(structNum, len(newSeq))

        cisFound = chooseStruct.removeCisStructs(structNum, newSeq, gro4)

        if not cisFound and chooseStruct.areAboveThresh(thresh, gro4, structNum): 
            foundPair = True

    #Calculate the rmsd of all non-cis structures and only keep the
    #most-different two structures
    chooseStruct.findMostDiff(newSeq, gro4)

    shutil.rmtree("noCis")
    shutil.rmtree("diffs")

    #Move the final structures and log files to a new directory named "newSeq"
    chooseStruct.moveBest(newSeq)

    #os.chdir("..")
    #prefix = input("What do you want as your numeral prefix")
    #os.mkdir(prefix + "_" + newSeq)

main()
