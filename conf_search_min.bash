#!/bin/bash -l
#$ -cwd -V
#$ -pe smp 16
#$ -N Z-hexathiocane
#$ -q common.q

export OMP_NUM_THREADS=$NSLOTS
export MKL_NUM_THREADS=$NSLOTS

module add xtb/6.3.2
module add anaconda3/5.0.1
source activate sci37


SMILESFILE=smiles-list.txt
NAMESFILE=names-list.txt


if test -f "$NAMESFILE"; then
    readarray names < $NAMESFILE
    nameexist=true
else
    echo "No $NAMEFILE file found. Naming shall be numeric."
fi
if test -f "$SMILESFILE"; then
    readarray lines < $SMILESFILE
    for i in "${lines[@]}"
    do
        echo $i
        for j in "${!lines[@]}"
        do
            if [ "${lines[$j]}" = "${i}" ] ; then
                if [ "$nameexist" = true ] ; then
                    dirname="${names[$j]}"
                else
                    dirname=$j
                fi
                mkdir $dirname
                echo $dirname
                cd $dirname
                python /home/uccabsa/scripts/conf_search_min.py -s $i -n 15000
                #mkdir crest
                #cp 1_*.xyz crest/
                #cd crest/
                #/home/uccabsa/Programs/crest/crest 1_*.xyz
                cd ../
            fi
        done
    done
else
    python /home/uccabsa/scripts/conf_search_min.py -s $1 -n 15000
    # mkdir crest
    # cp 1.xyz crest/
    # cd crest
    # /home/uccabsa/Programs/crest/crest 1.xyz
    cd ../
fi
