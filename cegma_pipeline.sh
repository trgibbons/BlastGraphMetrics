#!/bin/bash

################################################################################
#  Check that fragmentation scheme has been provided
################################################################################
if [ "${#1}" == 6 ]
then
    scheme=$1
else
    echo "Must specify a fragmentation scheme using a six-digit integer"
    exit 0
fi

################################################################################
#  Check environment
#  Initialize variables and top-level working directory
################################################################################
dir0=$(pwd)  # presumably ~/Dropbox/github/BlastGraphMetrics/
echo "Current working directory is:"
echo $dir0
echo

cegma=${dir0}/cegma/cegma.fasta
if [ -s "${cegma}" ]
then
    echo "CEGMA file found!"
    ls -lh $cegma
    echo
fi

dir1=${dir0}/cegma_${scheme}
echo "Creating and moving into subdirectory:"
echo $dir1
mkdir $dir1
cd $dir1
echo

################################################################################
#  Create test data sets by fragmenting CEGMA sequences
################################################################################
echo "Splitting CEGMA seqs..."
time nice ${dir0}/cegmaTestData.py $cegma $scheme
echo

echo "Created files:"
ls -lh cegma_${scheme}_???_???.fasta
echo

################################################################################
#  For each data set...
################################################################################
for distribution in {"ord","shf"}{"_evn","_rnd"}
do
    ############################################################################
    #  Determine FASTA file name from complete fragmentation scheme
    #  Create subdirectory and change into it to analyze this data set
    ############################################################################
    dir2=${dir1}/${distribution}
    echo "Creating and moving into subdirectory"
    echo $dir2
    mkdir $dir2
    fasta=cegma_${scheme}_${distribution}.fasta
    mv $fasta $dir2/
    cd $dir2
    echo

    ############################################################################
    #  ...format as a BLAST database
    ############################################################################
    echo "Formatting BLASTp databases"
    time nice makeblastdb -in $fasta -dbtype prot
    echo

    ############################################################################
    #  For each E-value...
    ############################################################################
    for e in $(seq 5 -2 3)
    do

        ########################################################################
        #  Determine BLAST output name from E-value
        #  Create subdirectory and change into it to analyze this BLAST run
        ########################################################################
        cutoff="1e-${e}"
        dir3=${dir2}/${cutoff}
        echo "Creating and moving into subdirectory"
        echo $dir3
        mkdir $dir3
        blastp=${fasta%.fasta}_${cutoff}.blastp
        cd $dir3
        echo

        ########################################################################
        #  ...BLAST data set against itself
        ########################################################################
        echo "BLAST'ing ${fasta} with max evalue cutoff 1e-${e}"
        time nice blastp \
            -query ${dir2}/${fasta} \
            -db ${dir2}/${fasta} \
            -out $blastp \
            -outfmt '7 std qlen slen' \
            -evalue $cutoff \
            -soft_masking true
        echo

        ########################################################################
        #  ...generate a set of MCL-formatted abc graphs
        ########################################################################
        echo "Generating abc graphs from ${blastp} for MCL"
        time nice ${dir0}/blast2graph.py $blastp ${blastp%.blastp}
        echo

        ########################################################################
        #  Analyze raw and normalized data sets seperately
        ########################################################################
        for normalization in raw nrm
        do

            ####################################################################
            #  Create subdirectory and chang into it
            ####################################################################
            dir4=${dir3}/${normalization}
            abc_pref=${blastp%.blastp}_${normalization}
            echo "Creating and moving into subdirectory:"
            mkdir $dir4
            mv ${abc_pref}_???.abc ${dir4}/
            cd $dir4
            echo

            ####################################################################
            #  For each graph...
            ####################################################################
            for metric in bit bpr bsr evl
            do
                abc=${abc_pref}_${metric}.abc
                echo "Creating clusters from ${abc}..."

                ################################################################
                #  For each inflation parameter...
                ################################################################
                #for I in $(seq -w 11 1 20; seq -w 25 5 50)
                for I in $(seq -w 11 1 50)
                do

                    ############################################################
                    #  ...cluster graph using MCL
                    ############################################################
                    i1=${I:0:1}
                    i2=${I:1}
                    echo "\tInflation: ${i1}.${i2}..."
                    time nice mcl \
                        $abc --abc \
                        -I ${i1}.${i2} \
                        -o ${abc%.abc}_I${I}.mcl
                    echo

                done
            done

            ####################################################################
            #  Print closing messages and generate Rtab file for BLAST file
            ####################################################################
            echo "Finished all MCL jobs for ${blastp}!"
            echo

            echo "Compiling Rtab file for ${abc_pref} clusters"
            time ${dir0}/mcl2rtab.py ${abc_pref} ${abc_pref}_???_I??.mcl
            echo

            echo "Generating heatmap PDFs for ${abc_pref}"
            time ${dir0}/ClusteringHeatmaps.R ${abc_pref}_*.Rtab

            echo "Moving back into directory:"
            echo $dir3
            cd $dir3
            echo

        done

        ########################################################################
        #  Move into new directory for next BLASTp E-value cutoff
        ########################################################################
        echo "Done processing ${fasta} with E-value cutoff ${cutoff}"
        echo
        echo "Moving back into directory:"
        echo $dir2
        cd $dir2
        echo

    done

    ############################################################################
    #  Print closing messages for data set
    ############################################################################
    echo "Finished all BLASTp jobs for ${fasta}!"
    echo

    echo "Moving back into directory:"
    echo $dir1
    cd $dir1
    echo

done

################################################################################
#  Consolidate PDFs into top-level working directory
################################################################################
echo "Moving all PDFs into:"
echo ${dir1}/heatmaps
mkdir ${dir1}/heatmaps
mv */*/*/*.pdf ./heatmaps/
echo

################################################################################
#  Move back into original working directory
################################################################################
echo "Moving back into original working directory:"
echo $dir0
cd $dir0
echo

echo "Finished!"
echo

exit 0

