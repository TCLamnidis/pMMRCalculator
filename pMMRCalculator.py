#!/usr/bin/env python3
import sys, argparse
from operator import add

VERSION = "1.0.2"

## A function to return the number of lines of a file
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

## A function to return the number of lines, as well as the last line with a SNP on chr 22. (assumes .snp file is sorted by ascending chr)
def file_len_snp(fname):
    switch=True
    with open(fname) as f:
        for i, l in enumerate(f):
            if switch and int(l.strip().split()[1]) > 22:
                stop_line=i
                switch=False
                # print(stop_line)
            pass
    return (i + 1, stop_line)

## A function to return the number of genotypes per line in a .geno file. 
def file_width(fname):
    with open(fname) as f:
        for i in f:
            return(len(i.strip()))
            break

## A function to calculate the pairwise mismatch rate of all pairs of genotypes from a string of genotypes.
def pMMR(genos, pMMRTable, ntot):
    for ind1 in range(0,len(genos)-1):
        if genos[ind1]=="9":
            continue
        geno1=int(genos[ind1])/2
        for ind2 in range(ind1+1,len(genos)):
            # print (ind1, ind2)
            newidx=ind2-(ind1+1)
            if genos[ind2] == "9":
                continue
            else:
                geno2 = int(genos[ind2]) / 2
                ntot[ind1][newidx] += 1
                pMMRTable[ind1][newidx] += geno1 * (1 - geno2) + (1 - geno1) * geno2
    return (pMMRTable, ntot)

## MAIN ##
## Parse arguments
parser = argparse.ArgumentParser(description="Calculate the pairwise mismatch rate of genotyped between all individuals in the input eigenstrat dataset.")
parser.add_argument("-i", "--Input", type=str, metavar="<INPUT FILES PREFIX>", required=False, help="The desired input file prefix. Input files are assumed to be <INPUT PREFIX>.geno, <INPUT PREFIX>.snp and <INPUT PREFIX>.ind .")
parser.add_argument("-o", "--Output", type=argparse.FileType('w'), metavar="<OUTPUT FILE>", required=False, help="The desired output file name. Omit to print to stdout.")
parser.add_argument("-s", "--Suffix", type=str, metavar="<INPUT FILES SUFFIX>", default='', required=False, help="The desired input file suffix. Input files are assumed to be <INPUT PREFIX>.geno<INPUT SUFFIX>, <INPUT PREFIX>.snp<INPUT SUFFIX> and <INPUT PREFIX>.ind<INPUT SUFFIX> .")
parser.add_argument("-v", "--version", action="store_true", help="Print the version of the script and exit.")
args = parser.parse_args()

## Print version and exit
if args.version:
    print(VERSION, file=sys.stderr)
    sys.exit(0)

## Check that an input file prefix was provided.
if args.Input == None:
    raise IOError("No input file prefix provided.")

## Open input files
IndFile = open(args.Input+".ind"+args.Suffix, "r")
GenoFile = open(args.Input+".geno"+args.Suffix, "r")

## Set output to stdout if the option is omitted.
if args.Output == None:
    args.Output = sys.stdout

## Calculate statistics for each input file.
linesGeno=[file_len(args.Input+".geno"+args.Suffix), file_width(args.Input+".geno"+args.Suffix)]
linesSnp=file_len_snp(args.Input+".snp"+args.Suffix)
linesInd=file_len(args.Input+".ind"+args.Suffix)

## Set last line with a SNP before chrX.
stop_line=linesSnp[1]

##Check geno and snp compatibility
if linesGeno[0] != linesSnp[0]:
    raise IOError("Input .snp and .geno files do not match.")

##Check geno and ind compatibility
if linesGeno[1] != linesInd:
    raise IOError("Input .ind and .geno files do not match.")

## Read individual name sfrom the .ind file.
Inds=[]
for line in IndFile:
    Ind=line.strip().split()[0]
    Inds.append(Ind)

## Initialise empty list with Ind-1 elements since last individual is already tested with everything.
ntot=[0 for x in Inds[:-1]]
pMMRTable=[0 for x in Inds[:-1]]
## Then replace each element in that list with a list of length equal to the pairwise comparisons per individual.
# print (Inds)
for i in range(len(Inds)-1):
    ntot[i]=[0 for x  in Inds[i+1:]]
    pMMRTable[i]=[0 for x  in Inds[i+1:]]

## Read through .geno file and calculate pMMR. every 25k lines, print debug info. Stop calculations after genotypes in autosomes are exhausted.
lineCount=0
for line in GenoFile:
    lineCount+=1
    (pMMRTable, ntot) = pMMR(line.strip(), pMMRTable, ntot)
    if lineCount % 25000 == 0:
        print("Processed",lineCount,"SNPs", sep=" ", file=sys.stderr)
    if lineCount >= stop_line:
        break
# print (ntot)
# print (pMMRTable)

## Print results
print ("Ind1","Ind2","nSNPs","nMismatch","pMismatch", sep="\t",file=args.Output)
for ind1 in range(0,len(Inds)-1):
    for ind2 in range(ind1+1,len(Inds)):
        newidx=ind2-(ind1+1)
        print (Inds[ind1],Inds[ind2], ntot[ind1][newidx], pMMRTable[ind1][newidx], "{:.5f}".format(pMMRTable[ind1][newidx]/ntot[ind1][newidx]), sep="\t", file=args.Output)
# print (Inds[1])
# print("INDS: ",linesInd, file=sys.stdout)
# print("SNPs: ",linesSnp, file=sys.stdout)
# print("GENOs: ",linesGeno, file=sys.stdout)
