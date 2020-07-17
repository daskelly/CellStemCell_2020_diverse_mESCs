#!/bin/bash

# Before running this script, pull the singularity image:
# singularity pull library://daskelly/remote-builds/tidyqtl2-r4.0.0:1.0.0
# in this directory

PREFIX=mapping_code
FILENAME=$PREFIX.r
OUTFILE=$PREFIX.RData
echo "source('map_qtl.r', echo=TRUE)" > $FILENAME
echo >> $FILENAME
echo "source('mediation.r', echo=TRUE)" >> $FILENAME
echo "save(eQTL, caQTL, mediation_results, file='$OUTFILE')" >> $FILENAME
singularity run --vm-ram 8192 tidyqtl2-r4.0.0_1.0.0.sif $FILENAME

rm -f $FILENAME
echo "Your output is saved in $OUTFILE"
