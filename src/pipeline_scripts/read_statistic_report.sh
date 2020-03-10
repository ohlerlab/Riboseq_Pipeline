#!/bin/bash

# Script for parsing STAR output logfile and fastqc reports,
# independent from any folder structure.

set -xv
set -e

DEBUG=0
FLAGS=(false false false)

while getopts :l:g:o: opt; do
    case $opt in
        l)
            STAR_LOG=$OPTARG
            FLAGS[0]=true
            ;;
        g)
            FASTQC_ORIGINAL=$OPTARG
            FLAGS[1]=true
            ;;
        o)
            OUT=$OPTARG
            FLAGS[2]=true
            ;;
        :)
            echo "Option -$OPTARG requires an argument" >&2
            exit 1
            ;;
        \?)
            echo "Unknown option -$OPTARG" >&2
            exit 1
            ;;
        *)
            echo "Unimplemented option -$OPTARG" >&2
            exit 1
            ;;
    esac
done

for FLAG in ${FLAGS[@]}; do
    if [ $FLAG == false ]; then
        echo -e "Usage:
\t-l STAR log file
\t-g Path to folder where zipped FASTQC reports are stroed (original)
\t-o Filename of output file" >&2
        exit 1
    fi
done

NAMES=$(mktemp)
DATA=$(mktemp)

_get_total_reads () {
    unzip -p $1 */fastqc_data.txt \
        | awk -F $'\t' '{if ($1=="Total Sequences") print $2}'
}

_get_read_length () {
    unzip -p $1 */fastqc_data.txt \
        | awk -F $'\t' '{if ($1=="Sequence length") print $2}'
}

get_total_reads () {
    for FILE in $1/*.zip; do
        _get_total_reads $FILE
    done | awk '{sum+=$1}END{print sum}'
}

get_read_length () {
    MIN=1000
    MAX=0
    for FILE in $1/*.zip; do
        RL=$(_get_read_length $FILE)
        CUR_MIN=$(echo $RL | cut -d- -f 1)
        CUR_MAX=$(echo $RL | cut -d- -f 2)
        if [ $MIN -gt $CUR_MIN ]; then
            MIN=$CUR_MIN
        fi
        if [ $MAX -lt $CUR_MAX ]; then
            MAX=$CUR_MAX
        fi
    done
    echo "$MIN-$MAX"
}

get_unmapped () {
    # unmapped reads are
    # [input reads] - ([uniquely mapped] + [multiple mapped] + [too many mapping positions])
    # there was a function without [too many mapping positions] but it didn't make any sense to me anymore
    echo $(sed -n '1p' $1) \
        - $(sed -n '3p' $1) \
        - $(sed -n '17p' $1) \
        - $(sed -n '19p' $1) \
        | bc -ql
}

export -f _get_total_reads
export -f _get_read_length
export -f get_total_reads
export -f get_read_length
export -f get_unmapped

awk 'BEGIN{FS="\t"}{print $1}' $STAR_LOG \
    | sed '1,5d;8d;23d;28d;s/^ *//;s/ |//' \
    > $NAMES

awk 'BEGIN{FS="\t"}{print $2}' $STAR_LOG \
    | sed '1,5d;8d;23d;28d;s/%//' \
    > $DATA

UNMAPPED=$(get_unmapped $DATA)
UNMAPPED_REL=$(echo "scale=2;print 0;100*$UNMAPPED/$(sed -n '1p' $DATA)" | bc -ql)

if [ $DEBUG -eq 1 ]; then
    echo This is a check number whether we calculated number of unmapped reads correctly
    echo \
        $UNMAPPED_REL =~ \
        $( \
            echo \
                $(sed -n '21p' $DATA) \
                + $(sed -n '22p' $DATA) \
                + $(sed -n '23p' $DATA) \
            | bc -lq \
        )
fi

TOTAL_READS_ORIGINAL=$(get_total_reads $FASTQC_ORIGINAL)

sed -i '20a\Number of unmapped reads' $NAMES
sed -i "20a\\$UNMAPPED" $DATA

sed -i '21a\Unmapped reads in %' $NAMES
sed -i "21a\\$UNMAPPED_REL" $DATA

sed -i '1i\Number of original reads' $NAMES
sed -i "1i\\$TOTAL_READS_ORIGINAL" $DATA

paste $NAMES $DATA > $OUT

rm -f $NAMES $DATA

