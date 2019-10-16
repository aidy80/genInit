#Aidan Fike
#July 21, 2019

#File with functions which measure omega angles and rmsd.

import os
import sys
import dihedDir
import shutil
import chimScriptMaker


#Function to measure omega angles of a given structure. Then, remove the structure 
#if it has cis bond(s) or move the structure if it doesn't have cis bond(s) to the noCis/ directory
#
#Params: structNum - int: The structure number being evaluated
#        newSeq - string: The sequence of the structure being evaluated 
def removeCisStructs(structNum, newSeq, gro4):
        if gro4:
            os.system("bash code/Sh_calc_omega.sh %d %d %d" % (structNum, len(newSeq), 4))
        else:
            os.system("bash code/Sh_calc_omega.sh %d %d %d" % (structNum, len(newSeq), 5))
        os.remove("VMD_GenOmegaIndex.sh")
        os.remove("struct_omega.xvg")
        os.remove("index_omega.ndx")
        os.remove("curr.gro")
        os.remove("angdist.xvg")
        cisFile = open("cisOut%d.txt" % structNum)
        cisFound = True
        for line in cisFile:
            if line.split()[0] == "No":
                cisFound = False
        cisFile.close()
        os.remove("cisOut%d.txt" % structNum)
        if cisFound:
            print "s%d was found to have a cis bond\n" % structNum
            os.remove("s%d_int_noh.pdb" % structNum)
        else:
            print "s%d is cis-bond free\n" % structNum
            os.rename("s%d_int_noh.pdb" % structNum, "noCis/s%d_int_noh.pdb" % structNum)
        return cisFound

#Function to calulate the rmsd between all non-cis sequences
#
#Return: Dictionary with the rmsd between each structure
def calcRmsd(gro4, shouldPrint):
    diffs = []
    di = "noCis/"
    for name1 in os.listdir("noCis"):
        for name2 in os.listdir("noCis"):
            if name1 != name2 and ({name1,name2} not in diffs):
                diffs.append({name1, name2})
                struct1 = name1.split("_")[0]
                struct2 = name2.split("_")[0]
                if not os.path.exists("diffs/" + struct1+struct2 + ".xvg"):
                    if gro4:
                        os.system("echo 4 4 | g_rms_mpi -what rmsd -fit rot+trans -s %s -f %s -o %s &> rms%s.out" \
                            % (di+ name1, di + name2, "diffs/" + struct1+struct2, struct1+struct2))
                    else:
                        os.system("echo 4 4 | gmx_mpi rms -what rmsd -fit rot+trans -s %s -f %s -o %s &> rms%s.out" \
                            % (di+ name1, di + name2, "diffs/" + struct1+struct2, struct1+struct2))

    dists = {}
    for name in os.listdir("diffs"):
        rms = open("diffs/%s" % name)
        structComb = name.split(".")[0]
        for line in rms:
            words = line.split()
            if words[0][0] != "#" and words[0][0] != "@":
                dists[structComb] = float(words[1])
                if shouldPrint:
                    print "\n%s and %s have a backbone rmsd of %f" % ("s" + structComb.split("s")[1], \
                                        "s" + structComb.split("s")[2], float(words[1]))

    return dists
     
#Find the largest rmsd between the created structues
#
#Params: dists: A dictionary with the rmsd between all of the structures
#
#Return: maxrmsd - int: The maximum rmsd found. If -1.0, there was < 2 structures with
#                       no cis bonds
#        maxpair - string: The structures with the max rmsd, i.e. s2s4
def findLargestRmsd(dists):
    maxrmsd = -1.0
    maxpair = ""
    for pair in dists.keys():
        if dists[pair] > maxrmsd:
            maxrmsd = dists[pair]
            maxpair = pair

    return maxrmsd, maxpair

#Function to move relevant log, rmsd, and dihedral files to the current directory
#
#Params: 
#        maxrmsd - float: the largest rmsd between sequences. -1.0 if < 2 structures with no cis
#                       bonds
#        maxpair - string: The sequences with the largest rmsd. I.e. s2s4
def moveRelevantFiles(maxrmsd, maxpair):
    if maxrmsd != -1.0:
        struct1 = "s" + maxpair.split("s")[1] 
        struct2 = "s" + maxpair.split("s")[2] 
        shutil.copyfile("%s.chimOut" % struct1, "s1Chim.log")
        shutil.copyfile("%s.chimOut" % struct2, "s2Chim.log")
        if (os.path.exists("%s.dchimOut" % struct1)):
            shutil.copyfile("%s.dchimOut" % struct1, "ds1Chim.log")
        if (os.path.exists("%s.dchimOut" % struct2)):
            shutil.copyfile("%s.dchimOut" % struct2, "ds2Chim.log")
        shutil.copyfile("%s.di" % struct1, "s1.dihed")
        shutil.copyfile("%s.di" % struct2, "s2.dihed")
        files = os.listdir(".")
        for name in files:
            if len(name) == 11 and name[0:3] == "rms":
                if name.split(".")[0][3:] == maxpair:
                    shutil.copyfile(name, name[0:3] + "s1s2" + ".log")

        os.rename("diffs/%s%s.xvg" % (struct1, struct2), "rmss1s2.xvg")
           
    os.system("rm *.chimOut")
    os.system("rm *.dchimOut")
    os.remove("edit.log")
    os.remove("cis.log")
    os.remove("vmd.log")
    os.system("rm *.di")

def areAboveThresh(thresh, gro4, structNum):
    #Calc rmsd between all non-cis structures
    dists = calcRmsd(gro4, True) 

    #Find largest rmsd
    maxrmsd, maxpair = findLargestRmsd(dists)


    if maxrmsd > thresh:
        return True
    else:
        if structNum >= 2:
            print "\nA pair of non-cis bonded structures with a backbone rmsd above the threshold is not\n\
yet found. Now creating a new structure to try out.\n"
        return False

#Function to calculate the rmsd between all non-cis structues and choose the two
#with the largest rmsd between them.
#
#Params: newSeq - string: The name sequence whose structures were made
def findMostDiff(newSeq, gro4):
    #Calc rmsd between all non-cis structures
    dists = calcRmsd(gro4, False) 

    #Find largest rmsd
    maxrmsd, maxpair = findLargestRmsd(dists)
    
    #Move log, dihed, and rmsd files around
    moveRelevantFiles(maxrmsd, maxpair)

    #Tell the user about the result of this function, and move the relevant pdb files
    if maxrmsd == -1.0:
        print "\nThere were one or less valid structures found - too many cis bonds were formed in the search process. \n\
Please try rerunning the program and hope for better luck. \n\
Additionally, increasing the numberInitStructs variable at the top of the main.py \n\
file will give you more chances to find non-cis structures\n\n"
        shutil.rmtree("noCis")
        shutil.rmtree("diffs")
        sys.exit()
    else:
        seqFile = open("finSeqs.txt", "a+")
        seqFile.write(newSeq + "\n")
        seqFile.close()

        struct1 = "s" + maxpair.split("s")[1] 
        struct2 = "s" + maxpair.split("s")[2] 

        print "\n%s and %s are going to be used as your pair of initial structures.\n \
They are written to s1_int_oneh.pdb and s2_int_oneh.pdb. \n \
They have an aligned-backbone rmsd of %f (nm)\n" % (struct1, struct2, maxrmsd)
        os.rename("noCis/%s_int_noh.pdb" % struct1, "s1_int_oneh.pdb")
        os.rename("noCis/%s_int_noh.pdb" % struct2, "s2_int_oneh.pdb")

    os.system("rm *.out")

#Move the final desired sequences and relevant log files to a new directory with the 
#name of the sequence you are creating.
#
#Params: newSeq - string: The name of the sequence you are creating
def moveBest(newSeq):
    if os.path.exists(newSeq):
        keep = -1
        while keep != "yes" and keep != "no":
            keep = raw_input("Your created sequence already has a directory made for it.\n\
Would you like to replace this directory ('yes' or 'no')\n\
If 'yes', the directory will be replaced. \n\
If 'no', the program will quit. \n")
        if keep == "yes": 
            shutil.rmtree(newSeq)
            os.mkdir(newSeq)
        if keep == "no":
            os.system("rm *.log") 
            os.system("rm *int_oneh.pdb")
            os.system("rm *.dihed") 
            os.system("rm *.xvg")
            sys.exit()
    else:
        os.mkdir(newSeq)

    print "Moving your structure s1_int_oneh.pdb and log files to %s/" % newSeq
    print "Moving your structure s2_int_oneh.pdb and log files to %s/" % newSeq
            
    os.system("mv *.log %s/" % newSeq)
    os.system("mv *.pdb %s/" % newSeq)
    os.system("mv *.dihed %s/" % newSeq)
    os.system("mv *.xvg %s/" % newSeq)

#Function to create the noCis and diffs directories, Removing noCis and iffs directories 
#if they already exist
def mkNoCisDiffs():
    if os.path.exists("noCis"):
        shutil.rmtree("noCis")
    if os.path.exists("diffs"):
        shutil.rmtree("diffs")
    os.mkdir("noCis")
    os.mkdir("diffs")
