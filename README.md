# pMMRCalculator
A small python tool to calculate pairwise mismatch rate between all individuals in an EigenStrat dataset.

_Mathematical equations added to README using [this tool](https://www.codecogs.com/latex/eqneditor.php)._

## Available options:
```
usage: pMMRCalculator.py [-h] [-i <INPUT FILES PREFIX>] [-o <OUTPUT FILE>]
                         [-s <INPUT FILES SUFFIX>] [-v] [-j]

Calculate the pairwise mismatch rate of genotyped between all individuals in
the input eigenstrat dataset.

Available options:
  -h, --help            show this help message and exit
  -i <INPUT FILES PREFIX>, --Input <INPUT FILES PREFIX>
                        The desired input file prefix. Input files are assumed
                        to be <INPUT PREFIX>.geno, <INPUT PREFIX>.snp and
                        <INPUT PREFIX>.ind .
  -o <OUTPUT FILE>, --Output <OUTPUT FILE>
                        The desired output file name. Omit to print to stdout.
  -s <INPUT FILES SUFFIX>, --Suffix <INPUT FILES SUFFIX>
                        The desired input file suffix. Input files are assumed
                        to be <INPUT PREFIX>.geno<INPUT SUFFIX>, <INPUT
                        PREFIX>.snp<INPUT SUFFIX> and <INPUT PREFIX>.ind<INPUT
                        SUFFIX> .
  -v, --version         Print the version of the script and exit.
  -j, --json            Create additional json formatted output file named
                        <OUTPUT FILE>.json . [Default:
                        'pmmrcalculator_output.json']
```

## Usage example:
```bash
pMMRCalculator.py -i test.input -o test.output.txt
```

## Example output:
```
Ind1	Ind2	nSNPs	nMismatch	pMismatch
Indiviudal1	Individual2	585385	137736.0	0.23529
Indiviudal1	Individual3	585199	140524.5	0.24013
Indiviudal1	Individual4	584162	140474.0	0.24047
Indiviudal2	Individual3	583507	140741.5	0.24120
Indiviudal2	Individual4	586613	141040.0	0.24043
Indiviudal3	Individual4	586005	140373.0	0.23954
```
The output will include a header and five coloumns. These columns contain the ID 
of the two individuals in each pairwise comparison, the number of overlapping 
SNPs between the two individuals, the cumulative mismatch proportions across all 
SNPs, and the pairwise mismatch rate between the individuals.

## Notes
The number of overlapping SNPs (`nSNPs`) shown here corresponds to the 
intersection of non-missing genotypes of the two individuals. 

The cumulative mismatch proportions are calculated with the following formula:
<br>
<br>
<a href="https://www.codecogs.com/eqnedit.php?latex=\bg_white&space;\fn_phv&space;\large&space;\mathit{n_{mismatch}=\sum_{S}^{&space;}\frac{G_{iS}}{2}&space;\times&space;(1&space;-&space;\frac{G_{jS}}{2})&space;&plus;&space;(1&space;-&space;\frac{G_{iS}}{2})&space;\times&space;\frac{G_{jS}}{2}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\bg_white&space;\fn_phv&space;\large&space;\mathit{n_{mismatch}=\sum_{S}^{&space;}\frac{G_{iS}}{2}&space;\times&space;(1&space;-&space;\frac{G_{jS}}{2})&space;&plus;&space;(1&space;-&space;\frac{G_{iS}}{2})&space;\times&space;\frac{G_{jS}}{2}}" title="\large \mathit{n_{mismatch}=\sum_{S}^{ }\frac{G_{iS}}{2} * (1 - \frac{G_{jS}}{2}) + (1 - \frac{G_{iS}}{2}) * \frac{G_{jS}}{2}}" /></a>

Where _G<sub>iS</sub>_ is the genotype of individual _i_ at SNP _S_, and _G<sub>jS</sub>_ is the genotype of individual _j_ at SNP _S_.
