# networkExpansion
Program for chemical network expansion. The software ilustrates the idea described in the paper `Computational synthesis design for controlled-degradation and revalorization.` by Anna Żądło-Dobrowolska et al.
Link to the paper and DOI will be added later.
## Usage 
### Options:

 -r   file with reactions
 
 -s   file with substrates
 
 -g   number of synthetic generations to run

 -o   output file, product of reactions together with number of generations with be stored in this file

 -v  verbose, print additional info to standard output
  
### Example usage: 
`python ./networkExpansion.py -r reaction.db -s substrates.smi  -g 1 -o outfile.txt`

both files: reaction.db (example of database) and substrates.smi (list of substrates) are provided and allow to perform simple calculation which produces one product after one generation.

## Note about reaction database

Reaction database should be provided as text file where each reaction is in new line. Reaction is stored as json string with following keys:
* idx - reaction unique name/identificator.
* rxSmarts - reaction SMARTS
* incompatGroups - list of functional groups, in SMARTS format, which are not allowed in substrates
* bannedProducts - list of functional groups, in SMARTS format, which are not allowed in product

See reaction.db for example.
