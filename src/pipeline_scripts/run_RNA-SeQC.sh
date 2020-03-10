#!/bin/bash

# module purge
# module load RNA-SeQC/v1.1.8-Java-1.7.0_80
# module load parallel/20150322-foss-2015a

FLAGS=(false false false false)
SINGLEEND=0
while getopts :r:t:i:o:e opt; do
    case $opt in
        r)
            REFERENCE=$OPTARG
            FLAGS[0]=true
            export REFERENCE
            ;;
        t)
            ANNOTATION=$OPTARG
            FLAGS[1]=true
            export ANNOTATION
            ;;
        i)
            INPUT+=($OPTARG)
            FLAGS[2]=true
            ;;
        o)
            OUT=$OPTARG
            FLAGS[3]=true
            ;;
	e)
	    SINGLEEND=1
	    ;;
        :)
            echo "Option -$OPTARG requires an argument" >&2
            exit 1
            ;;
        \?)
            echo "Unknown option -$OPTARG " >&2
            exit 1
            ;;
        *)
            echo "Error occured for -$OPTARG" >&2
            exit 1
            ;;
    esac
done

export SINGLEEND

for FLAG in ${FLAGS[@]}; do
    if [ $FLAG == false ]; then
        echo -e "Usage:
\t-i bam file to analyse, option can be used multiple times
\t-t annotation file
\t-r reference genome
\t-o output folder" >&2
        exit 1
    fi
done

TMPLIST=$(mktemp)
TMPOUT=$(mktemp)

# create input list from input parameters. one file each
for FILE in ${INPUT[@]}; do
    ID=$(basename $FILE)
    echo "$ID|$FILE|none" >> $TMPLIST
    echo "$OUT/RNA-SeQC" >> $TMPOUT
done

# create output directory (RNA-SeQC will fail if it doesn't exist)
mkdir -p $OUT

run_RNA_SeQC () {
    LIST=$1
    OUT=$2
    echo $LIST $OUT
    if [ $SINGLEEND -eq 1 ]
    then
	rna-seqc \
            -o $OUT \
            -s $LIST \
            -t $ANNOTATION \
            -r $REFERENCE \
	    -singleEnd
    else
	rna-seqc \
            -o $OUT \
            -s $LIST \
            -t $ANNOTATION \
            -r $REFERENCE
    fi
}

export -f run_RNA_SeQC

parallel --xapply 'run_RNA_SeQC {}' :::: $TMPLIST :::: $TMPOUT

rm $TMPLIST $TMPOUT

