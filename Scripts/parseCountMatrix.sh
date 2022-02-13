tail -n +2 ../Data/Monocytes-oxLDL_featurecounts.txt | cut -d $'\t' -f1,7- > ../Data/countmatrix.txt
