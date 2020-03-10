set -uev

files=''

pattern=${1:-''}

if [ "$pattern" ]; then 
	echo "filtering for files matching $pattern"
	files=$(snakemake -n --rerun-incomplete | grep -e 'rule all' -A1 | tail -n1 | sed 's/[,:]/\n/g' | grep -v 'input'  |  grep -Pe "$pattern" | xargs echo ) 
	
	if [ -z "$files" ]; then
		echo "no files matching"
		exit
	fi
fi

set -x 
snakemake -j 50  -k -p --restart-times 1 --max-jobs-per-second 5 -s Snakefile --cluster-config ../src/config_pipeline.json  --rerun-incomplete --use-conda --drmaa=" -cwd -V -l h_vmem={cluster.h_vmem} -l h_rt={cluster.h_rt} -pe {cluster.pe} -j yes -P medium -o sge_log" $files
