#!/bin/bash

################################################################################
#  My apologies to anyone having trouble running this script.
#
#  By the time I realized how complex this pipeline was going to become, I had
#  already finished writing most of it.  I made one abortive attempt to rewrite
#  it as a Ruffus pipeline in Python, in the hopes of improving it's
#  portability, robustness, and performance, while simultaneously making the
#  code more readable, but that turned out to be more work than I decided was
#  worth investing in a script that I couldn't really imagine anyone else using.
#
#  I suppose that as in nature, sometimes these things just happen:
#  http://en.wikipedia.org/wiki/Recurrent_laryngeal_nerve#Evidence_of_evolution
#
#  If you are having problems running this script, please do not hesistate to
#  contact me for help at ted.codes@gmail.com.
################################################################################

################################################################################
#  Check that fragmentation scheme has been provided
################################################################################
die () { echo "$@" 1>&2 ; exit 1 ; }

usage="3 arguments required:"
usage=$usage" [Input FASTA file]"
usage=$usage" [Output directory prefix]"
usage=$usage" [Fragmentation scheme]"
[ "$#" -eq 3 ] || die "${usage}, $# provided"

fragcheck="Third argument must be numeric"
echo $3 | grep -E -q '^[0-9]+$' || die "${fragcheck}, '$3' provided"

dirpref=$2
filepref=$(basename $dirpref)
scheme=$3

################################################################################
#  Create and move into top-level working directory
################################################################################
dir0=$(pwd)
echo "Current working directory is:"
echo $dir0
echo "All required scripts are assumed to be present."
echo

dir1=${dir0}/${dirpref}_${scheme}
echo "Creating and moving into subdirectory:"
echo $dir1
mkdir -p $dir1
cd $dir1
echo

################################################################################
#  Initialize log file
################################################################################
log=${dir1}/${filepref}_${scheme}.log
echo "Further status messages will be written to:"
echo $log
echo

date > $log
echo "Provided sequence partitioning scheme is:" >> $log
echo $scheme >> $log
echo >> $log

date >> $log
echo "Project directory is:" >> $log
echo $dir1 >> $log
echo >> $log

################################################################################
#  Check that the provided FASTA file exists and is not empty
################################################################################
eck=${dir0}/$1
if [ -s "${eck}" ]
then
    echo "Expanded CEGMA KOGs (ECK) database file found:" >> $log
    ls -lh $eck >> $log
    echo >> $log
else
    die "Expanded CEGMA KOGs (ECK) database ($eck) does not exist or is empty."
fi

################################################################################
#  Create test data sets by fragmenting ECK sequences
################################################################################
date >> $log
echo "Splitting ECK seqs..." >> $log
nice ${dir0}/eckTestData.py $eck $scheme ${filepref}
echo >> $log

date >> $log
echo "Created files:" >> $log
ls -lh ${filepref}_*.fasta >> $log
echo >> $log

################################################################################
#  For each data set...
################################################################################
for distribution in {"shf","ord"}{"_rnd","_evn"}
do
    ############################################################################
    #  Determine FASTA file name from complete fragmentation scheme
    #  Create subdirectory and change into it to analyze this data set
    ############################################################################
    dir2=${dir1}/${distribution}
    date >> $log
    echo "Creating and moving into subdirectory" >> $log
    echo $dir2 >> $log
    mkdir -p $dir2
    fasta=$(ls ${filepref}_*_${distribution}.fasta)
    mv $fasta $dir2/
    cd $dir2
    echo >> $log

    ############################################################################
    #  ...format as a BLAST database
    ############################################################################
    date >> $log
    echo "Formatting BLASTp databases" >> $log
    nice makeblastdb \
         -in $fasta \
         -dbtype prot \
         2>&1 >> $log
    echo >> $log

    ############################################################################
    #  For each E-value...
    ############################################################################
    #for e in $(seq 5 -2 3)
    for e in 5
    do

        ########################################################################
        #  Determine BLAST output name from E-value
        #  Create subdirectory and change into it to analyze this BLAST run
        ########################################################################
        cutoff="1e-${e}"
        dir3=${dir2}/${cutoff}
        date >> $log
        echo "Creating and moving into subdirectory" >> $log
        echo $dir3 >> $log
        mkdir -p $dir3
        blastp=${fasta%.fasta}_${cutoff}.blastp
        cd $dir3
        echo >> $log

        ########################################################################
        #  ...BLAST data set against itself
        ########################################################################
        date >> $log
        echo "BLAST'ing ${fasta} with max evalue cutoff 1e-${e}" >> $log
        nice blastp \
             -query ${dir2}/${fasta} \
             -db ${dir2}/${fasta} \
             -out $blastp \
             -outfmt '7 std qlen slen' \
             -evalue $cutoff \
             -soft_masking true \
             2>&1 >> $log
        echo >> $log

        ########################################################################
        #  ...generate a set of MCL-formatted abc graphs
        ########################################################################
        date >> $log
        echo "Generating abc graphs from ${blastp} for MCL" >> $log
        nice ${dir0}/blast2graph.py \
             $blastp \
             ${blastp%.blastp} \
             --fasta ${dir2}/$fasta
        echo >> $log

        ########################################################################
        #  Move connected component FASTA files into new directory
        ########################################################################
        date >> $log
        echo "Creating and moving connected component FASTA files into" >> $log
        echo "${dir3}/comp_fastas" >> $log
        mkdir -p ${dir3}/comp_fastas
        mv ${blastp%.blastp}_comp*.fasta ${dir3}/comp_fastas
        echo "Moving into subdirectory" >> $log
        echo ${dir3}/comp_fastas
        cd comp_fastas
        echo >> $log

        ########################################################################
        #  Generate MAFFT alignments and RAxML trees for each connected
        #  component in the BLAST graph
        ########################################################################
        for file in *.fasta
        do
            date >> $log
            echo "Generating MAFFT alignment for" >> $log
            echo $file >> $log
            mfa=${file%.fasta}.mfa
            nice mafft \
                  --auto \
                  --reorder \
                  --treeout \
                  $file \
                  > $mfa \
                  2>> $log
            echo >> $log

            date >> $log
            echo "Converting MAFFT FASTA alignment to Phylip format" >> $log
            echo $mfa "-->" >> $log
            echo $phy >> $log
            phy=${file%.fasta}.phy
            nice ${dir0}/fasta2phylip.py \
                 --cleanup \
                 $mfa \
                 $phy
            echo >> $log

# Don't have time at the moment
#            date >> $log
#            echo "Generating RAxML trees from MAFFT alignment" >> $log
#            nice raxmlHPC-PTHREADS-SSE3 \
#                 -s $phy \
#                 -n ${phy%.phy} \
#                 -m PROTGAMMAWAG \
#                 -p 42 \
#                 -T 8
#            echo >> $log
        done

        ########################################################################
        #  Finished analyzing connected component FASTA files
        ########################################################################
        date >> $log
        echo "Finished analyzing connected component FASTA files" >> $log
        echo "Moving back up into" >> $log
        echo $dir3 >> $log
        cd $dir3
        echo >> $log

        ########################################################################
        #  Analyze raw and normalized data sets seperately
        ########################################################################
        # Normalization (DiMensioNeD or DiMensionLesS)
        for norm in {"norm_wo","no_norm"}{"_dmnd","_dmls"}
        do
            ####################################################################
            #  Create subdirectory and change into it
            ####################################################################
            dir4=${dir3}/${norm}
            abc_pref=${blastp%.blastp}_${norm}
            date >> $log
            echo "Creating and moving into subdirectory:" >> $log
            echo $dir4 >> $log
            mkdir -p $dir4
            mv ${abc_pref}_???.abc ${dir4}/
            cd $dir4
            echo >> $log

            ####################################################################
            #  For each graph...
            ####################################################################
            for metric in bit bpl bsr pev
            do
                abc=${abc_pref}_${metric}.abc
#                mci=${abc_pref}_${metric}.mci
                date >> $log
                echo "Creating clusters from:" >> $log
                echo $abc >> $log
#                echo "Converting abc to mci file"
#                nice mcxload \
#                    --stream-mirror \
#                    --write-binary \
#                    -abc $abc \
#                    -o $mci
#                echo >> $log

                ################################################################
                #  For each inflation parameter...
                ################################################################
                for I in $(seq -w 11 1 60)
                do

                    ############################################################
                    #  ...cluster graph using MCL
                    ############################################################
                    i1=${I:0:1}
                    i2=${I:1}
                    date >> $log
                    echo "Running MCL with inflation parameter ${i1}.${i2}" \
                         >> $log
                    nice mcl \
                         $abc --abc \
                         -I ${i1}.${i2} \
                         -o ${abc%.abc}_I${I}.mcl \
                         2>&1 >> $log
                    echo >> $log

                done
            done

            ####################################################################
            #  Announce completion of this batch of MCL jobs
            ####################################################################
            date >> $log
            echo "Finished all MCL jobs for ${blastp}!" >> $log
            echo >> $log

            ####################################################################
            #  Generate barcharts using ggplot2 in R
            ####################################################################
            date >> $log
            echo "Compiling Rtab file for ${abc_pref} clusters" >> $log
            nice ${dir0}/mcl2rtab.py \
                 ${abc_pref} \
                 ${abc_pref}_???_I??.mcl
            echo >> $log

            date >> $log
            echo "Generating stacked barchart PDFs for:" >> $log
            echo ${abc_pref} >> $log
            nice ${dir0}/ClusteringBarcharts.R \
                 ${abc_pref}_*.Rtab
            echo >> $log

            ####################################################################
            #  Generate files for interactively viewing the results using either
            #  Cytoscape or Gephi
            ####################################################################
            date >> $log
            echo "Generating GML files for Cytoscape" >> $log
            nice ${dir0}/graphs2gml.py \
                 --gexf \
                 --graphml \
                 --compress bz2 \
                 --blast ${dir3}/${blastp} \
                 --graphs ${abc_pref}_???.abc \
                 --clusterings ${abc_pref}_???_I??.mcl \
                 --out_pref $abc_pref
            echo >> $log

            ####################################################################
            #  Announce end of analyses for this BLAST file
            ####################################################################
            date >> $log
            echo "Moving back into directory:" >> $log
            echo $dir3 >> $log
            cd $dir3
            echo >> $log

        done

        ########################################################################
        #  Move into new directory for next BLASTp E-value cutoff
        ########################################################################
        date >> $log
        echo "Done processing ${fasta} with E-value cutoff ${cutoff}" >> $log
        echo >> $log

        date >> $log
        echo "Moving back into directory:" >> $log
        echo $dir2 >> $log
        cd $dir2
        echo >> $log

    done

    ############################################################################
    #  Print closing messages for data set
    ############################################################################
    date >> $log
    echo "Finished all BLASTp jobs for ${fasta}!" >> $log
    echo  >> $log

    date >> $log
    echo "Moving back into directory:" >> $log
    echo $dir1 >> $log
    cd $dir1
    echo >> $log

done

################################################################################
#  Consolidate stacked barchart PDFs into top-level working directory
################################################################################
mkdir -p ${dir1}/barcharts
date >> $log
echo "Moving all PDFs into:" >> $log
echo ${dir1}/barcharts >> $log
mv */*/*/*.pdf ./barcharts/
echo >> $log

################################################################################
#  Consolidate Cytoscape files into top-level working directory
################################################################################
mkdir -p ${dir1}/cytoscape
date >> $log
echo "Moving all Cytoscape & Gephi files into:" >> $log
echo ${dir1}/cytoscape >> $log
mv */*/*/*.g*.bz2 ./cytoscape/
echo >> $log

################################################################################
#  Move back into original working directory
################################################################################
date >> $log
echo "Moving back into original working directory:" >> $log
echo $dir0 >> $log
cd $dir0
echo >> $log

date >> $log
echo "Finished!" >> $log
echo >> $log

exit 0

