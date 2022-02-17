#RNAseq analysis of macrophages exposed to oxLDL  

## Background

[Description of biological background and experimental procedures]  
[Sequencing by Novogene]


## Scripts
some of the scipts were written on windows and may contain DOS linebreaks. If this raises an error they can be 
converted to UNIX line breaks using the dos2unix command: `$ dos2unix script.run`
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

