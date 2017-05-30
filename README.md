# MLST
R Scipts used to determine MLST loci for Dichelobacter nodosus scheme
Database avalible at https://pubmlst.org/dnodosus/

Used to detemine number of sequence types from a random permutation of loci

Serogroup determination using IPCRESS
Slater, G.S.C. & Birney, E., 2005. Automated generation of heuristics for biological sequence comparison. BMC bioinformatics, 6, p.31

ipcress --input ipcress_serogroup.txt -p -P --sequence isolate.fasta > isolate_serogroup.txt
