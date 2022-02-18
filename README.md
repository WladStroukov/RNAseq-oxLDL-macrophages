# RNAseq analysis of macrophages exposed to oxLDL  

## Background

  1. Description of biological background
  Uptake of excessive cholestorol in the form of oxidized LDL through binding to the scavenger receptor on macrophages 
may lead to their transformation to foam cells and can contribute to the disease pathogenesis of atherosclerosis. 

  2. Experimental procedure

  3. Sequencing by Novogene]  
Used kits for RNA extraction  
Sequencer  


## Scripts
some of the scipts were written on windows and may contain DOS linebreaks. If this raises an error they can be 
converted to UNIX line breaks using the dos2unix command: `$ dos2unix script.run`  
./Scripts should be set as work directory and execute scripts from here.


  1. Quality control (by Novogene)
  2. Build STAR index GrCh38 version 104
  3. Map reads using STAR
  4. feature counts: use featurecounts script to remove first line and select column 1,7:end and save in new file. Make sure the parseCountMatrix script is exectuable using `$ chmod +x ./Scripts/parseCountMatrix.sh`
  5. DGE analysis with DESeq2 + functional analysis

[Explanation of scripts and processing steps]

## Environment 

[Required software and libraries]
[KCL HPC rosalind - add acknowledgements]


[set up environment]

