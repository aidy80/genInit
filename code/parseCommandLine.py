#Aidan Fike
#July 21, 2019

#File with function to parse command line arguments and extract relevant information. Also
#informs the user if they make a mistake

import sys
import os
import dihedDir

usage = "\nUSAGE: To use this program run:\n\n\
python main.py isRandom name/proteinLength threshold\n\n\
Where isRandom is whether or not you want the sequence you are creating to be \n\
randomly generated. If 'True', you want a random sequence, then the next argument \n\
should be the number of amino acids in your new sequence. Only (G, A, V, F, S, N, R, D) \n\
will be used for random sequences. Instead, if isRandom is given as 'False', the next \n\
argument should be the name of the protein you wish to create. This protein name can \n\
be make of any of the 21 common amino acids. \n    \
'threshold' is the minimium backbone-aligned rmsd difference allowed between your \n\
two generated structures. This should be given in units of [nm]. \n\
Lastly, you can run this program with gromacs4 rather than the gromacs5/gromacs2018 default. \n\
by adding a -gro4 flag after your other command line options.\n\n\
A couple examples are:\n    \
For a random 8mer with a minimum rmsd of .19nm created using gromacs4: \n        \
python main.py True 8 .19 -gro4\n   \
For a desired sequence, GNSRVGGGGG, with a minimum rmsd of .3nm using gromacs5/2018 you may want to use: \n        \
python main.py False GNSRVGGGGG .3 \n"

print "\nNOTE: To successfully run this program, you will need to \n\
load modules for chimera64/1.6.2 and a version of gromacs 5 or 2018. \n\
Alternatively, you can also use gromacs4 (only tested with 4.6.7), \n\
but to do so, you must include the -gro4 flag when you run the program. \n\
If you already are following these requirements, ignore this note.\n"

#Function to parse through all of the command line arguments. 
#
#Return: newSeq: string - The sequence whose structures will be created in this program
#        thresh: float - The backbone rmsd threshold
def getArgvInfo():
    newSeq = ""
    gro4 = False
    if len(sys.argv) < 4 and not (len(sys.argv) == 2 and sys.argv[1] == "-h"):
        print usage,"\n"
        print "ERROR: You didn't input enough command line arguments\n"
        sys.exit()

    if sys.argv[1] == "-h":
        print usage, "\n"
        sys.exit()

    if sys.argv[1] == "True":
        numAminos = sys.argv[2]
        if not str.isdigit(numAminos):
            print usage, "\n"
            print "ERROR: You didn't correctly specify the number of amino acids you wanted in your new random sequence.\n"
            sys.exit()

        numAminos = int(numAminos)

        newSeq = dihedDir.createNewSeq(numAminos) 

    elif sys.argv[1] == "False":
        newSeq = sys.argv[2]
        if not (sys.argv[2].isupper() and sys.argv[2].isalpha() and dihedDir.isAllowedSeq(sys.argv[2])):
            print usage, "\n"
            print "ERROR: The sequence you wanted to create is not possible with this program. \n\
It contains uncommon amino acids. Next time, please input a sequence with only uppercase letters representing\n\
allowed amino acids (A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y)\n"
            sys.exit()
        elif len(newSeq) == 1:   
            print "ERROR: This program does not work for sequences of only 1 amino acid, sorry\n"
            sys.exit()


    else:
        print usage, "\n"
        print "ERROR: You didn't choose isRandom correctly. Choose either 'True' or 'False' for\n\
        whether the desired sequence should be generated randomly (True) or given manually in the next argument (False)'\n"
        sys.exit()

    thresh = sys.argv[3]
    if not isFloat(thresh):
        print usage, "\n"
        print "\nERROR: %s is not a valid decimal-formatted number. \n\
Next time you run, this argument should be used to input the threshold rmsd value you want. \n\
This will become the minimium backbone-aligned backbone rmsd between your two structures\n" % thresh
        sys.exit()

    if "-gro4" in sys.argv:
        gro4 = True
    elif len(sys.argv) > 4:
        print usage, "\n"
        print "ERROR: The command %s was no recognized\n" % sys.argv[4]
        sys.exit()

    if len(sys.argv) > 5:
        for arg in sys.argv[5:]:
            print usage, "\n"
            print "ERROR: The command %s was no recognized\n" % arg
            sys.exit()
         
    return newSeq, float(thresh), gro4

#Check if a given user input is a float
def isFloat(x):
    try:
        float(x)
        return True
    except:
        return False
