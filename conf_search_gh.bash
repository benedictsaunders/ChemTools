#!/bin/bash -l

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
                python conf_search.py -s $i -n 15000
                cd ../
            fi
        done
    done
else
    while getopts s:n: flag
    do
        case "${flag}" in
            s) smiles=${OPTARG};;
            n) numconfs=${OPTARG};;
        esac
    done
    python conf_search_min.py -s $smiles -n $numconfs
    cd ../
fi
