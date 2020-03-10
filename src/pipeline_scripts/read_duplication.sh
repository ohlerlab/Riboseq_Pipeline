#!/bin/bash

# module purge
# module load R/3.2.0-foss-2015a-bare
# module load RSeQC/2.6.2-foss-2015a-Python-2.7.9
# module load parallel/20150322-foss-2015a

FLAGS=(false false)

while getopts :i:o: opt; do
    case $opt in
        i)
            INPUT+=($OPTARG)
            FLAGS[0]=true
            ;;
        o)
            OUT=$OPTARG
            FLAGS[1]=true
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

for FLAG in ${FLAGS[@]}; do
    if [ $FLAG == false ]; then
        echo -e "Usage:
\t-i bam file (can be used multiple times)
\t-o output folder" >&2
        exit 1
    fi
done

TMPINPUT=$(mktemp)
TMPOUT=$(mktemp)

for FILE in ${INPUT[@]}; do
    ID=$(basename $FILE)
    echo "$FILE" >> $TMPINPUT
    echo "$OUT/$ID" >> $TMPOUT
done

mkdir -p $OUT

run_read_duplication () {
    IN=$1
    OUT=$2
    read_duplication.py \
        -i $IN \
        -o ${OUT} \
        2>/dev/null
}

export -f run_read_duplication

parallel --xapply run_read_duplication :::: $TMPINPUT :::: $TMPOUT

rm $TMPINPUT $TMPOUT
