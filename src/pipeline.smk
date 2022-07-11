
import glob
import pandas as pd
from pathlib import Path
from functools import partial
# from ipdb import set_trace

singularity: "docker://dermotharnett/riboseq_pipeline"

#print('warning - sampling the first 100k reads');HEADIFTEST = '| head -n 400000'
HEADIFTEST = ''


def is_nonempty(file):
  assert Path(file).stat().st_size
def is_over_size(file,n):
  assert Path(file).stat().st_size > n
def newfolder(file,newroot):
  file = Path(file)
  assert snakedir in file.parents , str(snakedir) + " doesn't seem to be in parents of " + str(file)
  return str(Path(*(newroot,)+file.relative_to(snakedir).parts[1:]))

#shell.executable("/bin/bash")
shell.prefix("set -e pipefail;")
# user set parameter


################################################################################
########load and check  pipeline configuration
################################################################################
 

configfile: "../config/config.yaml"

# PROJECTFOLDER=config['PROJECTFOLDER'] not used
# TMPDIR = Path('../tmp') not used

seqfilesdf = pd.read_csv(config['sample_config'],dtype=str).set_index("sample_id", drop=False)
### Creates Boolean
seqfilesdf.isriboseq=seqfilesdf.isriboseq.map({'True':True,'False':False})
#sampledf = pd.read_csv(config['sample_parameter']).set_index("sample_id", drop=False) old code from seperated config files, delete if script runs ok


###??? what is the use if this whole block till line 63?
assert seqfilesdf.sample_id.is_unique
lacks_mate_info = (not 'mate' in seqfilesdf.columns) or (seqfilesdf.mate.isna().all())
if lacks_mate_info: seqfilesdf['mate'] = '1'
lacks_pair_id = (not 'pair_id' in seqfilesdf.columns) or (seqfilesdf.pair_id.isna().all())
if lacks_pair_id:
    assert seqfilesdf['mate'].isin(['1']).all()
    seqfilesdf['pair_id']=seqfilesdf['sample_id']
lacks_file_id = not 'file_id' in seqfilesdf.columns or (seqfilesdf.file_id.isna().all())
if lacks_file_id: seqfilesdf['file_id']=seqfilesdf.pair_id+'_R'+seqfilesdf.mate+'.fastq.gz'
assert (~seqfilesdf.file_id.isna()).all()

assert isinstance(seqfilesdf.iloc[0,1],str), "file column should be a string in read_files.csv"
assert 'file_id' in seqfilesdf.columns
assert not (pd.Series([Path(f).name for f in seqfilesdf.file_id]).duplicated().any()),"files need unique filenames"

seqfilesdf.mate = seqfilesdf.mate.fillna('1')
assert set(seqfilesdf.mate).issubset(set(['1','2']))

seqfilesdf.pair_id = seqfilesdf.pair_id.fillna(seqfilesdf.sample_id+'.fastq.gz')
seqfilesdf.file_id = seqfilesdf.file_id.fillna(seqfilesdf.sample_id+'_'+seqfilesdf.mate+'.fastq.gz')

assert seqfilesdf.file_id.is_unique,"pairid + mate combo must be unique"

print('found '+str(len(seqfilesdf.sample_id))+' sample ids')
print('for '+str(len(seqfilesdf))+' files')
n_paired = str(sum(seqfilesdf.groupby('pair_id').size() > 1))
print('of these, '+n_paired+' were paired')

for i in seqfilesdf.file: assert Path(i).stat().st_size > 100, "This file isn't of size 100 or more " +i

#For creating a pair ID column
# pairids = [Path(f).name.replace('_R2_','_R1_').replace('.fastq.gz','') for f in seqfilesdf.file]
# pairids = [re.sub('_[A-Za-z0-9]{8,20}$','',p) for p in pairids]
# #insert tinot he read file
# seqfilesdf.insert(seqfilesdf.shape[1],'pair_id', pairids, False)
#seqfilesdf.to_csv(config['sample_files'])


###??? this line we can take out when we combine sample_files and sample_parameters
#make sure our ids and files look ok
assert set(seqfilesdf.sample_id) == set(seqfilesdf.sample_id), "Sample IDs need to be setequal in "+config['sample_files']+" and "+config['sample_parameter'] +": \n"+"seqfile ids " + seqfilesdf.sample_id[0:3] + "... \n" +"seqfile ids " + seqfilesdf.sample_id[0:3] + "... "

for sample in seqfilesdf.sample_id:
  for f in seqfilesdf.loc[[sample],'file']:
    fp = Path(f)
    assert fp.exists, f
    assert 'fastq.gz' in fp.name or 'fq.gz' in fp.name or 'fastq' in fp.name , f

#define samples
samples = list(seqfilesdf['sample_id'].unique())
fastqs = list(seqfilesdf['file'].unique())

assert seqfilesdf.isriboseq.isin([True,False]).all()
ribosamples = seqfilesdf.sample_id[seqfilesdf.isriboseq]
rnasamples = seqfilesdf.sample_id[~seqfilesdf.isriboseq]

#for trimming CDS for riboseq
REF_orig=config['REF_orig']
GTF_orig=config['GTF_orig']
ORFEXT_FASTA=Path(GTF_orig.replace('.gtf','.orfext.fa')).name

assert(Path(GTF_orig).exists()), GTF_orig + ", the GTF file, doesn't exist"

#local copies of the annotation
REF = Path(Path(re.sub(string=REF_orig,pattern='.(b)?gz$',repl=r'')).name)
GTF = Path(Path(re.sub(string=GTF_orig,pattern='.(b)?gz$',repl=r'')).name)

# RNAFASTA = GTF.with_suffix('.fa') not used
# CODINGFASTA=GTF.with_suffix('.coding.fa') not used
# PROTEINFASTA=GTF.with_suffix('.protein.fa') not used
# CDSFASTA=GTF.with_suffix('.cds.fa') not used
# BED=GTF.with_suffix('.bed') not used
# PCFASTA=config['PCFASTA'] not used
PCFASTA_tname = config['PCFASTA'].replace('.fa','.shortheader.fa')

ALIGNERS = ['star']
ALIGNER_TO_USE = 'star'

#this configures whether we want to trim ids from e.g. ENSG0000001.1 to ENSG000001
if config['TRIM_IDS']: 
  mod_id_sed_cmd = r''' sed -r 's/((transcript_id|gene_id|protein_id|ID|Parent|exon_id|havana_gene|havana_transcript)\W+\w+)\.[0-9]+/\1/g' '''
else:
  mod_id_sed_cmd = ' cat '

assert seqfilesdf.libtype.str.match('[IOM]?[SU][FR]?').all()


rule all:
  input:
    seqfilesdf.file,
    expand("{aligner}/data/{sample}/{sample}.bam", aligner = ALIGNERS, sample = samples),
    # ("multiqc/multiqc_report.html"),
    expand("ORFquant/{sample}/.done", sample = ribosamples),
    expand('riboseqc/data/{sample}/.done', sample=ribosamples),
    expand('salmon/data/{sample}/quant.sf',sample=rnasamples),
    expand('ribostan/{sample}/{sample}.ribostan.tsv', sample=ribosamples)


## copy_ref: create a plain text copy of the reference genome file in
## case it is in a compressed format and use SAMtools to generate its
## index. These will be used later for mapping with STAR.

rule copy_ref:
  input: REF_orig
  output: REF,str(REF)+'.fai'
  shell: """
      zless {REF_orig} > {output[0]}
      samtools faidx {output[0]}
      """


## link_in_anno: check the annotation file for consistency with the
## reference genome file in naming chromosomes.

rule link_in_anno:
  input: REF_orig=REF_orig,GTF=GTF,REFAI=str(REF)+'.fai'
  output: touch('annocheck.done')
  run:
    from pathlib import Path
    REFchrs = pd.read_csv(input.REFAI,sep='\t',header=None,names=['Chr','length','cumlength','V4','V5'],usecols=[0,1])
    REFchrs.in_ref = True
    #awk '{if(! /#/ ){f[$1]=$5}}END{for(i in f){ print i,f[i]}}'
    import subprocess
    from subprocess import Popen, PIPE
    from io import StringIO


    cmd=r"""awk '{if(! /#/ ){f[$1]=$5}}END{for(i in f){ print i,f[i]}}' """+str(GTF_orig)
    a=subprocess.Popen(cmd,stdout = subprocess.PIPE,shell=True)
    GTFchrs = pd.read_csv(StringIO(a.communicate()[0].decode('ascii')),sep=' ',header=None,names=['Chr','longest_in_GTF'])

    allchrs = pd.merge(REFchrs,GTFchrs,how='outer')
    allchrs['outofbounds'] = allchrs.longest_in_GTF.gt(allchrs.length) 
    allchrs['missing_in_REF'] = allchrs.length.isna()
    problems = allchrs.outofbounds | allchrs.missing_in_REF
    refonly = allchrs.longest_in_GTF.isna()
    assert not problems.any(),print("\n\n\n GTF, Chromosomes Don't quite match here: \n\n",allchrs[problems|refonly],"\n\n")

    cmd=r"""awk '{if(! /#/ ){f[$1]=$5}}END{for(i in f){ print i,f[i]}}' """+str(GTF_orig)
    a=subprocess.Popen(cmd,stdout = subprocess.PIPE,shell=True)
    GTFchrs = pd.read_csv(StringIO(a.communicate()[0].decode('ascii')),sep=' ',header=None,names=['Chr','longest_in_GTF'])

    allchrs = pd.merge(REFchrs,GTFchrs,how='outer')
    allchrs['outofbounds'] = allchrs.longest_in_GTF.gt(allchrs.length) 
    allchrs['missing_in_REF'] = allchrs.length.isna()
    problems = allchrs.outofbounds | allchrs.missing_in_REF
    assert not problems.any(),print("\n\n\n GTF, Chromosomes Don't quite match here: \n\n",allchrs[problems],"\n\n")


## link_in_files: make a preprocessed_reads/ directory for each sample
## and place a symbolic link to the unprocessed raw data (i.e.
## *.fastq.gz file) in each directory to have a nice project structure.

def name_preprocessed_reads(wc):
  filedf = seqfilesdf[seqfilesdf.sample_id==wc['sample']]  
  assert (filedf.file_id  == wc['fastq']).sum()==1, "Fastq file isn't amongst the file ids"
  filedf = filedf[filedf.file_id==wc['fastq']]
  return filedf.file

rule link_in_files:
  input: name_preprocessed_reads
  output: 'preprocessed_reads/{sample}/{fastq}'
  run:  
    sample = wildcards['sample']
    fastq = wildcards['fastq']
    shell(r"""
      mkdir -p $(dirname {output})
      ln -sf $(readlink -f {input} {output} )
    """)


## cutadapt_reads: use cutadapt to trim off the adapter sequence
## specified as ADAPTERSEQ in config.yaml from all reads. Resulting
## trimmed reads that are shorter than MINREADLENGTH or longer than
## MAXREADLENGTH are discarded. Furthermore, bases on the 3' end with
## quality scores lower than QUALLIM are trimmed off.

ADAPTERSEQ=config['ADAPTERSEQ']
MINREADLENGTH=config['MINREADLENGTH']
MAXREADLENGTH=config['MAXREADLENGTH']
QUALLIM=config['QUALLIM']

rule cutadapt_reads:
  input: 'preprocessed_reads/{sample}/{fastq}'
  output: 'cutadapt_reads/{sample}/{fastq}'
  log: 'cutadapt_reads/{sample}/{fastq}.cutadaptstats.txt'
  shell: """    
       set -ex
       
       mkdir -p cutadapt_reads/{wildcards.sample}/
        zless {input} \
           {HEADIFTEST} \
           | cutadapt \
             -a {ADAPTERSEQ} \
            --minimum-length {MINREADLENGTH} \
            --maximum-length {MAXREADLENGTH} \
            -q {QUALLIM} - \
        2> {log}  \
        | gzip  > {output}.tmp



      #delete the zipped file if it's empty
      gzip -l {output}.tmp  | awk 'NR==2 {{exit( $2 != 0) }}' && rm {output}.tmp

      mv {output}.tmp {output}
"""


## collapse_reads: collapse duplicate reads based on their UMI sequence
## (if you are using UMIs).

rule collapse_reads:
    input: 'cutadapt_reads/{sample}/{fastq}'
    output: 'collapse_reads/{sample}/{fastq}'
    params: 
    shell: r"""
       set -evx
       mkdir -p collapse_reads/{wildcards.sample}/
       zcat {input}  \
         | ../src/pipeline_scripts/collapse_reads.pl {wildcards.sample} \
         2> collapse_reads/{wildcards.sample}/{wildcards.fastq}.collreadstats.txt \
         | gzip > {output}
     """


## trim_reads: trim off UMI sequences from reads after collapsing
## duplicates (if you are using UMIs).

rule trim_reads:
    input: 'collapse_reads/{sample}/{fastq}'
    output: 'trim_reads/{sample}/{fastq}'
    params:
      outdir = lambda wc,output: Path(output[0]).parent
    shell: r"""
        set -evx
        OUTDIR=$(dirname {output})
        mkdir -p  $OUTDIR
        zcat {input} > {input}.nozip
        ../src/pipeline_scripts/remove8N.pl {input}.nozip {output}
        gzip -f {output}
        mv {output}.gz {output}
     """


##########################################################
######### trna/rRNA filtering with STAR
##########################################################


## number_contaminants: read in fasta file of tRNA and rRNA contaminant
## sequences specified as Contam in config.yaml. Assign a number to
## each unique sequence. Create a new fasta file of these numbered
## sequences within the pipeline structure.

rule number_contaminants:
  input: contaminants=config['Contam']
  output: 'tRNA_rRNA_index/tRNA_rRNA_index.fa',
  threads: 8
  params:
    fafile = lambda wc,output: output[0].replace('.done','')+'.fa'
  shell: r"""
      #number the contaminant sequences
    R -e '
      library(Biostrings); 
      seq = readDNAStringSet("{input.contaminants}");
      seq = unique(seq);
      seq = setNames(seq,paste0(seq_along(seq),"_",names(seq)));
      writeXStringSet(seq,"{output}");
      write.table(col.names=F,row.names=F,data.frame(names(seq)),"tRNA_rRNA_index.names.txt");
    ' 
    """


## make_trna_rrna_indices: use STAR to generate index of tRNA and rRNA
## contaminant sequences. We are handling the contaminant fasta file as
## a reference genome file.

rule make_trna_rrna_indices:
  input: 'tRNA_rRNA_index/tRNA_rRNA_index.fa'
  output: touch('tRNA_rRNA_index/tRNA_rRNA_index.done'),
  # conda: '../envs/star'
  threads: 8
  params:
    outprefix = lambda wc,output: output[0].replace('/tRNA_rRNA_index.done',''), 
  shell: r"""
    STAR \
    --genomeSAindexNbases 9 \
    --runThreadN {threads} \
    --runMode genomeGenerate \
    --genomeDir {params.outprefix} \
    --genomeFastaFiles {input}
    """


## filter_tRNA_rRNA: use STAR to map reads to contaminant sequences.
## Use SAMtools to discard reads that successfully mapped to
## contaminant sequences.

rule filter_tRNA_rRNA: 
  input: 'trim_reads/{sample}/{fastq}','tRNA_rRNA_index/tRNA_rRNA_index.done'  
  output: 'filter_reads/{sample}/{fastq}'
  threads: 8
  params:
    # indexname = lambda wc,input: input[1].replace('.done',''), not used
    genomedir = lambda wc,input: input[1].replace('/tRNA_rRNA_index.done',''),
    outdir = lambda wc,output: os.path.dirname(output[0]),
    nozip = lambda wc,output: output[0].replace('.gz','')
  shell: r"""
    #set -evx
    set -x
    [ -f {params.outdir} ] && rm -rf {params.outdir}
   
    mkdir -p {params.outdir}

    STAR \
      --genomeDir {params.genomedir} \
      --runThreadN {threads} \
      --readFilesCommand zcat \
      --outSAMunmapped Within \
      --outMultimapperOrder Random \
      --outFilterMultimapNmax 1 \
      --alignSJoverhangMin 8 \
      --outTmpDir {params.outdir}/_tmpSTAR \
      --alignSJDBoverhangMin 1 \
      --genomeLoad NoSharedMemory \
      --outSAMattributes NH HI AS NM MD \
      --outSAMtype SAM \
      --outFileNamePrefix {output[0]}.filtered_reads.sam \
      --outReadsUnmapped Fastx \
      --readFilesIn {input[0]}

    mv {output[0]}.filtered_reads.samUnmapped.out.mate1 {params.nozip}
    gzip {params.nozip}
    mv {output[0]}.filtered_reads.samAligned.out.sam {output[0]}.filtered_reads.sam

    samtools view -bh {output[0]}.filtered_reads.sam \
    | samtools sort -@ {threads} > {output[0]}.filtered_reads.bam
    samtools index {output[0]}.filtered_reads.bam
    
      #those which mismatch twice should not be included
    samtools view -hb {output[0]}.filtered_reads.bam \
    | bamtools filter -tag XM:2-10 -in - -out /dev/stdout \
    | samtools view -H > {output[0]}.mm.sam
    #>> {params.outdir}/unmapped.sam
   
      #group the idx columns stuff is from 
    samtools idxstats {output[0]}.filtered_reads.bam \
    | perl -lanpe 's/^(\S+)_[^_\s]+\t/$1\t/' > {output[0]}.idxtmp

    samtools stats {output[0]}.filtered_reads.bam > {output[0]}.filtered_reads.bam.stats
    #samtools stats {output[0]}.filtered_reads.bam
    """
    #rename unmapped to the output


## link_processed_reads: create a single directory for each sample with
## all relevant fastq files of processed reads inside. This rule is the
## signal splitter where we go from sample to individual fastq files.

def choose_processed_reads(wc,config=config):
  #correct zcat strings based on read pairs
  filedf = (seqfilesdf.loc[[wc['sample']]])
  isrna = ~seqfilesdf.loc[wc['sample'],'isriboseq']
  filedf = filedf[filedf.file_id==wc.fileid]
  if isrna:
    return [config['FILT_RNA_FOLDER']+'/'+wc['sample']+'/'+f for f in filedf.file_id]
  else:
    return [config['FILT_RIBO_FOLDER']+'/'+wc['sample']+'/'+f for f in filedf.file_id]

rule link_processed_reads:
  input: choose_processed_reads 
  output: 'processed_reads/{sample}/{fileid}'
  run: 
    shell(r"""
        echo choose_processed_reads
        mkdir -p $(dirname {output} )
        for i in $(readlink -f {input} ); do  ln -rifs $i  {output} ; done
    """)
    assert Path(output[0]).stat().st_size > 100

################################################################################
########Annotation
################################################################################


## makeGTF: create a plain text GTF file copy of the annotation file
## even if it is in a compressed format.

rule makeGTF:
  input: GTF=GTF_orig
  output: GTF
  shell: r""" 
      # set -x
      #with filtering output all sequences
      zless {input.GTF} \
      | {mod_id_sed_cmd} \
      | gffread -F -T -o {GTF}

    """


## make_utrs: extract all features marked as 5'UTRs or 3'UTRs from the
## annotation files and produce separate GTF files for those.
 
rule make_utrs:
  input: GTF=GTF_orig
  output: fputrs='fputrs.gtf',tputrs='tputrs.gtf'
  # script: 'make_utrfiles.R'
  shell: r"""
      set -ex
      #with filtering output all sequences
      cat {input.GTF}  \
      | awk -v OFS="\t"  '{{if($3=="five_prime_UTR"){{         ;print $0}}}}' \
      | sed -r  's/((transcript_id|gene_id|protein_id|ID=\w+|Parent)\W+\w+)\.[0-9]+/\1/g' \
      > {output.fputrs} 

      cat {input.GTF} \
      | awk -v OFS="\t"  '{{if($3=="three_prime_UTR"){{         ;print $0}}}}' \
      | sed -r  's/((transcript_id|gene_id|protein_id|ID=\w+|Parent)\W+\w+)\.[0-9]+/\1/g' \
      > {output.tputrs}
      """


################################################################################
########Fastqc
################################################################################


## fastqc: run processed reads (reads that have been adapater-trimmed,
## duplicate-collapsed and UMI-trimmed, and contaminant-filtered)
## through FastQC.

def get_sample_fastqs(wc,mate='1',folder='processed_reads',seqfilesdf=seqfilesdf):
   #correct zcat strings based on read pairs
  filedf = (seqfilesdf.loc[[wc['sample']]])
  filedf = filedf.loc[filedf.mate==mate,]
  isrna = ~seqfilesdf.loc[wc['sample'],'isriboseq']
  assert isrna.all() | (~isrna).all()
  isrna = isrna.all()
  folder =  config['FILT_RNA_FOLDER']if isrna else config['FILT_RIBO_FOLDER']
  matefiles = [folder+'/'+wc['sample']+'/'+f for f in filedf.file_id]
  return matefiles

get_sample_fastqs2 = partial(get_sample_fastqs,mate='2')

rule fastqc:
     input:
        lfastqs=get_sample_fastqs,
        rfastqs=get_sample_fastqs2,
     output: touch('fastqc/data/{sample}/.done')
     threads: 4
     log:'fastqc/reports/{sample}/fastqc.log'
     params:
        outdir = lambda wc: 'fastqc/data/'+wc.sample+'/'
     shell: '''
          OUTDIR=$(dirname {output[0]})
          mkdir -p {params.outdir}
          wait $(for i in {input.lfastqs} {input.rfastqs}; do $( fastqc -o {params.outdir} $i ) & done) 
        '''


## collect_fastqc: summarize FastQC results of all samples in a single
## tab-separated file.
###??? why don't we do this anymore :(
rule collect_fastqc:
     input:
          all_results = expand("fastqc/data/{sample}/.done", sample=samples)
     output:
          result='fastqc/summary/fastqc_summary.tsv',
          log='fastqc/summary/fastqc_summary.log'
     shell:
          r"""
          set -e
          mkdir -p $(dirname {output.result}) 
          {SCRIPTDIR}/collect_fastqc_results.sh -i fastqc/ \
          > {output.result} \
          2> {output.log} 
          """


################################################################################
########STAR
################################################################################


## star_index: use STAR to index the reference genome and annotation, a
## necessary step before mapping with STAR.

rule star_index:
 input: REF=ancient(REF),GTF=ancient(GTF)
 output: touch('starindex/starindex.done')
 threads: 8
 shell: r"""
   STAR \
     --runThreadN {threads} \
     --runMode genomeGenerate \
     --genomeDir $(dirname {output}) \
     --sjdbGTFfile {input.GTF} \
     --genomeFastaFiles {input.REF}
   """


## star: use STAR to map processed reads (reads that have been
## adapater-trimmed, duplicate-collapsed and UMI-trimmed, and
## contaminant-filtered) to the indexed reference genome. Use SAMtools
## to sort and index the resulting bam files.

#mating paired end reads
def get_file_string(wc,seqfilesdf=seqfilesdf):
# #correct zcat strings based on read pairs
 filedf = (seqfilesdf.loc[[wc['sample']]])

 mate1files = filedf.loc[filedf.mate==1,'file'].str.replace('input','processed_reads').values
 mate2files = filedf.loc[filedf.mate==2,'file'].str.replace('input','processed_reads').values
 lfilestring = '<(zcat '+' '.join(mate1files)+')'
 rfilestring = '<(zcat '+' '.join(mate2files)+')' if mate2files   else '' 
 filestring = lfilestring+' '+rfilestring
 return(filestring)

rule star:
     input:
          lfastqs=get_sample_fastqs,
          rfastqs=get_sample_fastqs2,
          STARINDEX='starindex/starindex.done',
     output:
          done = touch('star/data/{sample,[^/]+}/.done'),bam='star/data/{sample}/{sample}.bam',bai='star/data/{sample}/{sample}.bam.bai'
     threads: 8
     #please check params, there is a lot of stuff i did.
     params:
        sample= lambda wc: wc['sample'],
        #filestring = lambda wc: get_file_string(wc,seqfilesdf),
        GEN_DIR=lambda wc,input: input.STARINDEX.replace('starindex.done',''),
        #only used for remap (now remap='')
        # markdup = lambda wc: '' if seqfilesdf.isriboseq[wc['sample']] else '-m', not used?
        platform = 'NotSpecified',
        outputdir = lambda wc,output: os.path.dirname(output[0]),
        repdir = lambda wc,output: os.path.dirname(output[0]).replace('data','reports'),
        #tophatindex =lambda wc,input: input['bowtie_index'].replace('.done',''),
        halfthreads = lambda wc,threads: threads/2,
        sortmem = lambda wc,threads: str(int(5000/(threads/2)))+'M',
        #remap = '1' if seqfilesdf.isriboseq[wc['sample']] else ''
        # remap = '',
        lfilestring = lambda wc,input: '<(zcat '+' '.join(input.lfastqs)+')',
        rfilestring = lambda wc,input: '<(zcat '+' '.join(input.rfastqs)+')' if input.rfastqs   else '' ,
     shell: r"""
        set -x
        MY_TMP_DIR=$(mktemp -d)
        
        trap "set -x; rm -rf ${{MY_TMP_DIR}}" EXIT KILL TERM INT HUP

        mkdir -p $MY_TMP_DIR
        mkdir -p $MY_TMP_DIR/star
        mkdir -p $MY_TMP_DIR/tophat2



        #--outSAMmultNmax 20 --winAnchorMultimapNmax 50 --outFilterMultimapNmax 20 --genomeDir {input.STARINDEX} \

        STAR \
              --genomeDir {params.GEN_DIR} \
              --runThreadN {threads} \
              --outSAMunmapped Within \
              --outFilterType BySJout \
              --outMultimapperOrder Random \
              --alignSJoverhangMin 8 \
              --alignSJDBoverhangMin 1 \
              --outFilterMismatchNmax 999 \
              --outFilterMismatchNoverLmax 0.04 \
              --alignIntronMin 20 \
              --alignIntronMax 1000000 \
              --alignMatesGapMax 1000000 \
              --genomeLoad NoSharedMemory \
              --quantMode GeneCounts \
              --outSAMattributes NH HI AS NM MD \
              --outSAMtype BAM Unsorted\
              --outSAMattrRGline \"ID:{sample}\" \"SM:{sample}\" \"PL:{params.platform}\" \
              --outFileNamePrefix ${{MY_TMP_DIR}}/star/ \
              --outReadsUnmapped Fastx \
              --readFilesIn {params.lfilestring} {params.rfilestring}
          
         mv ${{MY_TMP_DIR}}/star/Aligned.out.bam ${{MY_TMP_DIR}}/star/all.bam

         samtools sort \
          -m {params.sortmem} \
          -T ${{MY_TMP_DIR}} \
          -o {params.outputdir}/{wildcards.sample}.bam \
          -@ {params.halfthreads}\
          ${{MY_TMP_DIR}}/star/all.bam
      
        samtools index {params.outputdir}/{wildcards.sample}.bam 

        mkdir -p {params.repdir}
        samtools stats {params.outputdir}/{wildcards.sample}.bam > {params.repdir}/{wildcards.sample}.bamstats.txt
        samtools flagstat {params.outputdir}/{wildcards.sample}.bam > {params.repdir}/{wildcards.sample}.flagstat.log
        samtools idxstats {params.outputdir}/{wildcards.sample}.bam > {params.repdir}/{wildcards.sample}.idxstats.log
        
        cp ${{MY_TMP_DIR}}/star/ReadsPerGene.out.tab {params.outputdir}/ReadsPerGene.out.tab
        cp ${{MY_TMP_DIR}}/star/SJ.out.tab {params.outputdir}/
        cp ${{MY_TMP_DIR}}/star/{{Log.final.out,Log.out}} {params.repdir}/
        if [ $iftophat ] ;then cp ${{MY_TMP_DIR}}/tophat2/align_summary.txt {params.repdir} ;fi

          """


################################################################################


## make_picard_files: generate a text file containing sequence
## alignment map headers and coordinates of rRNA sequences to be used
## by Picard. This is done by using SAMtools to extract headers from a
## bam file produced by STAR and running grep to find rRNA genes on the
## reference annotation GTF file. Additionally, use gtfToGenePred to
## convert the GTF into a genePred table.

rrna_intervals = 'qc/picard_rrna_intervals.txt'
refflat = Path('qc') / Path(GTF).with_suffix('.refflat').name
#refflat = snakedir/ 'qc' / Path(config['GFF_orig']).with_suffix('.refflat').name

rule make_picard_files:
  input: GTF=GTF,bam=ALIGNER_TO_USE+'/data/'+samples[0]+'/'+samples[0]+'.bam'
  output: intervals=rrna_intervals,refflat=refflat
  shell:r"""
        set -x 
          
        samtools view -H {input.bam} > {output.intervals}
        
        grep -Pe 'gene_type..rRNA.' {input[0]} \
        | awk '$3 =="transcript"' \
        | cut -f 1,4,5,7,9 \
        | perl -lane ' /transcript_id "([^"]+)"/ or die "notranscript_id on $."; print join "\t", (@F[0,1,2,3], $1) ' \
        | sort -k1V -k2n -k3n - >> {output.intervals}
        
        gtfToGenePred -geneNameAsName2 {input.GTF} {GTF}.genepred
        cat {input.GTF}.genepred | awk -vOFS="\t" '{{print $1,$0}}' > {output.refflat}

  """


## qc: use Picard and custom scripts to generate quality control
## reports and plots of mapped reads. Results are placed in the qc/
## directory.

rule qc:
     input:
          fastqc='fastqc/data/{sample}/.done',
          bam=ALIGNER_TO_USE+'/data/{sample}/{sample}.bam',
          refflat = refflat,
          rrna_intervals = rrna_intervals,
          readstatscript= config['rnaseqpipescriptdir']+"read_statistic_report.sh",
          read_duplication= config['rnaseqpipescriptdir']+"read_duplication.sh",
     output:
          'qc/data/{sample}/read_alignment_report.tsv',done=touch('qc/data/{sample}/.done')
     resources:
     params:
        singleendflag = lambda wc: ' -e ' if seqfilesdf.libtype.str.match('^[SU]')[wc.sample]  else '',
        scriptdir = lambda wc: config['rnaseqpipescriptdir']
     shell: """
          set -exv
          
        OUTDIR=$(dirname {output.done})
        mkdir -p qc/reports/{wildcards.sample}/


        {params.scriptdir}/run_RNA-SeQC.sh -i {input.bam} {params.singleendflag} -t {GTF} -r {REF} -o ${{OUTDIR}}/

        {input.readstatscript} \
         -l star/reports/{wildcards.sample}/Log.final.out \
         -g $(dirname {input.fastqc}) \
         -o ${{OUTDIR}}/read_alignment_report.tsv \
         &> qc/reports/{wildcards.sample}/{wildcards.sample}_qc.log 

          if [ ! -s qc/data/{wildcards.sample}/read_alignment_report.tsv ]; then 
            echo 'alignments report is empty'
            rm qc/data/{wildcards.sample}/read_alignment_report.tsv
            exit
          fi

         picard CollectRnaSeqMetrics -Xms4G \
          I={input.bam} \
          O=${{OUTDIR}}/{wildcards.sample}_picard_qc.txt \
          REF_FLAT={refflat} \
          STRAND=FIRST_READ_TRANSCRIPTION_STRAND \
          RIBOSOMAL_INTERVALS={rrna_intervals}
        
        picard CollectAlignmentSummaryMetrics -Xms4G \
          INPUT={input.bam} \
          OUTPUT=${{OUTDIR}}/{wildcards.sample}.picard.alignmentmetrics.txt \
          R={REF}

        {input.read_duplication} \
         -i {input.bam} \
         -o ${{OUTDIR}}/duplication/

          """


## multiqc: use MultiQC to collate reports produced by the following
## rules: fastqc, star, qc, and salmon into one nice html report.

def get_multiqc_dirs(wildcards,input):
      reportsdirs = list(input)
      reportsdirs=[s.replace(ALIGNER_TO_USE+'/data',ALIGNER_TO_USE+'/reports') for s in reportsdirs]
      reportsdirs=[s.replace('tophat2/data','tophat2/reports') for s in reportsdirs]
      reportsdirs=[os.path.dirname(s) for s in list(reportsdirs)]
      assert len(reportsdirs) > 0
      return(reportsdirs)

rule multiqc:
  input:
      expand("fastqc/data/{sample}/.done", sample = samples),
      expand(ALIGNER_TO_USE+"/data/{sample}/{sample}.bam", sample = samples),
      expand("qc/data/{sample}/.done", sample = samples),
      # expand("tophat2/data/{sample}/.done", sample = samples),
      # [f.replace('input','filter_reads') for f in  seqfilesdf.file[ribosamples]],
      expand("salmon/data/{sample}/.done", sample = samples),
      # 'sample_file.txt'
  #conda: '../envs/multiqc'
  params: 
    multiqcscript = config['multiqcscript'],
    sample_reads_file=seqfilesdf.file,
    reportsdirs= get_multiqc_dirs,
    sampnames = '--sample-names '+config.get('samplenamesfile') if config.get('samplenamesfile') else ''
  output:
    'multiqc/multiqc_report.html'
  shell:r"""
      cat {params.sample_reads_file} | sed 's/.fastq.gz//g' | sed 's/\t.*\//\t/g' \
      | awk -vOFS='\t' 'BEGIN{{print "fastqname","samplename"}}{{sumsamp[$1] = sumsamp[$1]+1;print $2,$1"_fq"sumsamp[$1]}}' \
      > multiqc/samplenames.txt

      {params.multiqcscript} {params.reportsdirs} -fo $(dirname {output[0]}) {params.sampnames}
      """

RIBOSEQCPACKAGE = config['RIBOSEQCPACKAGE']
rule make_riboseqc_anno:
  input : 
    GTF,
    REF,str(REF)+'.fai'
  output: GTF.with_suffix('.matchchrs.gtf_Rannot')
  params:
    annobase = Path(GTF).name.replace('.gtf','').replace('.gz',''),
    gtfmatchchrs = lambda wc,input: input[0].replace('.gtf','.matchchrs.gtf')
  shell:r"""
    set -x
    awk -vOFS="\t" '{{print $1,0,$2}}' {REF}.fai | bedtools intersect -b - -a {GTF} > {params.gtfmatchchrs} 
    mkdir -p $(dirname {output[0]})
    R -e 'devtools::load_all("{RIBOSEQCPACKAGE}");args(prepare_annotation_files) ;prepare_annotation_files(annotation_directory=".",gtf_file="{params.gtfmatchchrs}",annotation_name="{params.annobase}",forge_BS=FALSE, genome_seq=FaFile("{REF}"))'
 """

rule run_riboseqc:
   input: GTF.with_suffix('.matchchrs.gtf_Rannot'),bam=ALIGNER_TO_USE+'/data/{sample}/{sample}.bam'
   output: touch('riboseqc/data/{sample}/.done'),'riboseqc/data/{sample}/_for_ORFquant','riboseqc/reports/{sample}/riboseqcreport.html',
   threads: 4
   params:
     annofile = lambda wc,input: input[0].replace('annot.done',Path(GTF_orig).name.replace('.gtf','.matchchrs.gtf_Rannot')),
     outname = lambda wc,output: output[0].replace('.done',''),
     report_file = lambda wc: 'riboseqc/reports/'+wc['sample']+'/'+'riboseqcreport.html',
   shell:r"""
         set -x
         mkdir -p {params.outname}
         mkdir -p riboseqc/reports/{wildcards.sample}
         R -e 'devtools::load_all("{RIBOSEQCPACKAGE}");RiboseQC::RiboseQC_analysis("{params.annofile}", bam="{input.bam}",rescue_all_rls=TRUE,dest_names="{params.outname}", genome_seq = "{REF}", report_file="{params.report_file}")'
     """

#####Added pulling ORFquant from the config
ORFquantPACKAGE = config['ORFquantPACKAGE']
rule run_ORFquant:
  input : 'riboseqc/data/{sample}/.done',annofile=GTF.with_suffix('.matchchrs.gtf_Rannot')
  output: touch('ORFquant/{sample}/.done')
  params: 
    for_ORFquantfile = lambda wc,input: 'c("'+('","'.join(['riboseqc/data/'+s+'/''_for_ORFquant' for s in [wc['sample']] ]))+'")',

    outputdir = lambda wc,output: output[0].replace('.done','')
  threads: 10
  shell:r"""
    set -ex
      mkdir -p {params.outputdir}

      R -e 'devtools::load_all("{RIBOSEQCPACKAGE}");devtools::load_all("{ORFquantPACKAGE}");run_ORFquant(for_ORFquant_file = {params.for_ORFquantfile},annotation_file = "{input.annofile}", n_cores = {threads},prefix="{params.outputdir}") '
      """

##########################################################################
########Transcript Alignmnents
################################################################################

rule make_orf_fasta:
  input: gtf=GTF_orig,fasta=REF
  output: fasta=ORFEXT_FASTA
  shell:r"""set -ex
  mkdir -p $(dirname {output.fasta})
  R -e 'devtools::load_all("{RIBOSTANPACKAGE}");make_ext_fasta(gtf="{input.gtf}",fa="{input.fasta}",out="{output.fasta}")'
  """

rule star_transcript_index:
  input: fasta = lambda wc:  ORFEXT_FASTA if wc['sequences'] == 'ORFext' else PCFASTA_tname
  output: directory('StarIndex/{sequences}')
  # conda: "../envs/star"
  threads: 15
  shell:r"""
    mkdir -p {output}
    STAR  \
    --runThreadN {threads} \
    --runMode genomeGenerate \
    --genomeDir {output} \
    --genomeFastaFiles {input.fasta} \
    --genomeSAindexNbases 11 \
    --genomeChrBinNbits 12
  """ 

rule star_transcript:
  input:
    fastq = get_sample_fastqs,
    transcriptindexfold = 'StarIndex/{sequences}',
  output: 
    bam='star/{sequences}/data/{sample}/{sample}.bam',
    bai='star/{sequences}/data/{sample}/{sample}.bam.bai',
  # conda: "../envs/star"
  params:
    bamnosort=lambda wc,output: output.bam.replace('.bam','nosort.bam')
  shell:r"""

  mkdir -p star_transcript/data/{wildcards.sample}/{wildcards.sample}.fastq.gz

  STAR --runThreadN {threads} --genomeDir {input.transcriptindexfold} \
    --readFilesIn <(zcat {input.fastq}) \
    --outFileNamePrefix star_transcript/{wildcards.sample}.transcript_ \
    --outSAMtype BAM Unsorted \
    --outSAMmode NoQS \
    --outSAMattributes NH NM \
    --seedSearchLmax 10 \
    --outFilterMultimapNmax 255 \
    --outFilterMismatchNmax 1 \
    --outFilterIntronMotifs RemoveNoncanonical

  mv star_transcript/{wildcards.sample}.transcript_Aligned.out.bam {params.bamnosort}
  samtools sort {params.bamnosort} -o {output.bam}
  rm {params.bamnosort}
  samtools index {output.bam}

  """


rule ntrim_pc_fasta:
  input: config['PCFASTA']
  output: PCFASTA_tname
  shell:r"""perl -lanpe 's/\|.*$//' {input}   > {output}"""

################################################################################
########Quantification
################################################################################
 

rule make_salmon_index:
  threads: 8
  input: config['PCFASTA']
  output: salmonindex = touch('salmonindex/.done')
  # conda: '../envs/salmon'
  shell:r""" salmon index -p {threads}  -k 21 -t {input} -i salmonindex"""

rule salmon:
  input:fastqs=get_sample_fastqs,salmonindex='salmonindex/.done'
  params:
    # salmonindex = lambda wc,input: 'salmonindexribo' if wc['sample'] in ribosamples else input.salmonindex.replace('.done',''),
    salmonindex = 'salmonindex',
    lib = lambda wc: seqfilesdf.loc[wc['sample'],'libtype']
  output:
      quant = touch('salmon/data/{sample}/quant.sf')
  threads: 4
  # conda:'../envs/salmon'
  shell:r"""
      set -ex
      mkdir -p salmon/reports/{wildcards.sample}
      mkdir -p salmon/data/{wildcards.sample}
      salmon quant \
      -p {threads} \
      -l {params.lib} \
      --seqBias \
      -i {params.salmonindex} \
      -r <(zcat {input.fastqs} ) \
      --output salmon/data/{wildcards.sample} \
      --validateMappings
"""


RIBOSTANPACKAGE = config['RibostanPACKAGE']

rule ribostan:
  input:
    ribobam = 'star/ORFext/data/{sample}/{sample}.bam',
    ribofasta = ORFEXT_FASTA
  threads:4
  # conda: '../envs/ribostan'
  output: efile = 'ribostan/{sample}/{sample}.ribostan.tsv'
  shell: r"""
    mkdir -p $( dirname {output.efile} ) 
    R -e 'devtools::load_all("{RIBOSTANPACKAGE}");get_exprfile("{input.ribobam}", "{input.ribofasta}", "{output.efile}")'
  """


################################################################################
########Downstream analysis
################################################################################
rule read_countdata:
  input:
    salmon = expand('salmon/data/{sample}/quant.sf', sample=ribosamples),
    ribostan = expand('ribostan/{sample}/{sample}.ribostan.tsv', sample=rnasamples),
    GTF=GTF,
    REF=REF
  threads:4
  output: 
    'r_data/tx_countdata.rds',
    'r_data/iso_tx_countdata.rds',
    'r_data/gtf_gr.rds'
  shell: r"""
    Rscript ../src/read_countdata.R --gtf={input.GTF} --fafile={input.REF}
  """

rule run_deseq:
  threads:4
  input:  
    'r_data/tx_countdata.rds',
    'r_data/iso_tx_countdata.rds',
    'r_data/gtf_gr.rds'
  output:
    'r_data/dds.rds'
  shell: r"""
  Rscript ../src/run_deseq.R 
  """

rule deseq_report:
  input:
   'r_data/tx_countdata.rds','r_data/dds.rds'
  threads:4
  output: report = 'deseq_report/ribo_de_report.html'
  shell: r"""
  R -e '
  rmarkdown::render(
          "ribo_de_report.Rmd",
          output_format = "html",
          output_file="{output.report}",
          output_dir = dirname("{output.report}"),
          intermediates_dir = dirname("{output.report}"),
          knit_root_dir = dirname("{output.report}")
  )
  '
  """


  ### Make UMI and no UMI version
  ### 'Rsamtools','rtracklayer','GenomicAlignments','BSgenome','GenomicFiles','reshape2','DT','ggpubr','viridis','GenomicFeatures' for riboseqc
