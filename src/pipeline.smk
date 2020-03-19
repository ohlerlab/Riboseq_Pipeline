
import glob
import pandas as pd
from pathlib import Path
from functools import partial
from ipdb import set_trace

#### conda install -c bioconda subread needed
#### conda install -c bioconda ucsc-gtftogenepred same
#### conda install -c bioconda fastqc

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
 

configfile: "../src/config.yaml"
config['root'] = Path(config['root'])

snakedir = Path(config['root']).resolve() / 'pipeline'
assert snakedir.exists
TMPDIR = Path('../tmp')

seqfilesdf = pd.read_csv(config['sample_files'],dtype=str).set_index("sample_id", drop=False)
sampledf = pd.read_csv(config['sample_parameter']).set_index("sample_id", drop=False)

assert sampledf.sample_id.is_unique
if not 'mate' in seqfilesdf.columns: seqfilesdf['mate'] = '1'
if not 'pair_id' in seqfilesdf.columns:
    assert seqfilesdf['mate'].isin(['1']).all()
    seqfilesdf['pair_id']=seqfilesdf['sample_id']

if not 'file_id' in seqfilesdf.columns: seqfilesdf.insert(seqfilesdf.shape[1],'file_id',seqfilesdf.pair_id+'_R'+seqfilesdf.mate+'.fastq.gz',False)


assert sampledf.sample_id.is_unique
assert isinstance(seqfilesdf.iloc[0,1],str), "file column should be a string in read_files.csv"
assert not (pd.Series([Path(f).name for f in seqfilesdf.file]).duplicated().any()),"files need unique filenames"

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

#make sure our ids and files look ok
assert set(seqfilesdf.sample_id) == set(sampledf.sample_id), "Sample IDs need to be setequal in "+config['sample_files']+" and "+config['sample_parameter'] +": \n"+"seqfile ids " + seqfilesdf.sample_id[0:3] + "... \n" +"seqfile ids " + sampledf.sample_id[0:3] + "... "

for sample in sampledf.sample_id:
  for f in seqfilesdf.loc[[sample],'file']:

    fp = Path(f)
    assert fp.exists, f
    assert 'fastq.gz' in fp.name or 'fq.gz' in fp.name , f


#define samples
samples = list(sampledf['sample_id'].unique())
fastqs = list(seqfilesdf['file'].unique())
ribosamples = sampledf.sample_id[sampledf.assay.str.lower()=='ribo']
rnasamples = sampledf.sample_id[sampledf.assay.str.lower()=='total']

#for trimming CDS for riboseq
STARTCODTRIM=config['STARTCODTRIM']
STOPCODTRIM=config['STOPCODTRIM']
REF_orig=config['REF_orig']
GTF_orig=config['GTF_orig']

assert(Path(GTF_orig).exists()), GTF_orig + ", the GTF file, doesn't exist"

#local copies of the annotation
REF = Path(Path(re.sub(string=REF_orig,pattern='.(b)?gz$',repl=r'')).name)
GTF = Path(Path(re.sub(string=GTF_orig,pattern='.(b)?gz$',repl=r'')).name)

RNAFASTA = GTF.with_suffix('.fa')
CODINGFASTA=GTF.with_suffix('.coding.fa')
PROTEINFASTA=GTF.with_suffix('.protein.fa')
CDSFASTA=GTF.with_suffix('.cds.fa')
BED=GTF.with_suffix('.bed')


#For now we'll just use star
#ALIGNERS = ['hisat2','star']
ALIGNERS = ['star']
ALIGNER_TO_USE = 'star'

#this configures whether we want to trim ids from e.g. ENSG0000001.1 to ENSG000001
if config['TRIM_IDS']: 
  mod_id_sed_cmd = r''' sed -r 's/((transcript_id|gene_id|protein_id|ID|Parent|exon_id|havana_gene|havana_transcript)\W+\w+)\.[0-9]+/\1/g' '''
else:
  mod_id_sed_cmd = ' cat '


rule all:
  input:
    seqfilesdf.file,
    expand("{aligner}/data/{sample}/{sample}.bam", aligner = ALIGNERS, sample = samples),
    ("multiqc/multiqc_report.html"),
    expand('feature_counts_readrange/data/{sample}/{gcol_generegion}/{readrange}/feature_counts', sample=ribosamples, gcol_generegion='gene_id__CDS', readrange=config['RNALENRANGE']),
    expand('feature_counts_readrange/data/{sample}/{gcol_generegion}/{readrange}/feature_counts', sample=rnasamples, gcol_generegion='gene_id__CDS', readrange=config['RIBOLENRANGE']),
    expand("ORFquant/{sample}/.done", sample = ribosamples),
    expand('riboseqc/data/{sample}/.done', sample=ribosamples)

MINREADLENGTH=config['MINREADLENGTH']
MAXREADLENGTH=config['MAXREADLENGTH']
QUALLIM=config['QUALLIM']
REMOVE8NBIN=config['REMOVE8NBIN']

rule copy_ref:
  input: REF_orig
  output: REF,str(REF)+'.fai'
  shell: """
      zless {REF_orig} >  {output[0]}
      samtools faidx {output[0]}
      """

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


ADAPTERSEQ=config['ADAPTERSEQ']

rule cutadapt_reads:
  input: 'preprocessed_reads/{sample}/{fastq}'
  output: 'cutadapt_reads/{sample}/{fastq}'
  #conda: '../envs/cutadapt'
  log: 'cutadapt_reads/{sample}/{fastq}.cutadaptstats.txt'
  shell: """    
       set -ex
       
       mkdir -p cutadapt_reads/{wildcards.sample}/
        zcat {input} \
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

rule collapse_reads:
    input: 'cutadapt_reads/{sample}/{fastq}'
    output: 'collapse_reads/{sample}/{fastq}'
    run:
        sample = wildcards['sample']
        collapse_reads_script = config['collapse_reads_script'],

        shell(r"""
       set -evx

       mkdir -p collapse_reads/{sample}/
     
       zcat {input}  \
         | {collapse_reads_script} {wildcards.sample} \
         2> collapse_reads/{wildcards.sample}/{wildcards.fastq}.collreadstats.txt \
         | gzip > {output}
     """)
        is_over_size(output[0],100)

#
rule trim_reads:
    input: 'collapse_reads/{sample}/{fastq}'
    output: 'trim_reads/{sample}/{fastq}'
    params:
      outdir = lambda wc,output: Path(output[0]).parent
    run:
        sample = wildcards['sample']
        shell(r"""
          set -evx
          OUTDIR=$(dirname {output})
          mkdir -p  $OUTDIR
          zcat {input} > {input}.nozip
          ../src/pipeline_scripts/remove8N.pl {input}.nozip {output}
          gzip -f {output}
          mv {output}.gz {output}
     """)
#


##########################################################
######### trna/rRNA filtering with bowtie
##########################################################
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

rule make_trna_rrna_indices:
  input: 'tRNA_rRNA_index/tRNA_rRNA_index.fa'
  output: touch('tRNA_rRNA_index/tRNA_rRNA_index.done'),
  conda: '../envs/star'
  threads: 8
  params:
    outprefix = lambda wc,output: output[0].replace('/tRNA_rRNA_index.done',''), 
  shell: r"""
    STAR \
    --runThreadN {threads} \
    --runMode genomeGenerate \
    --genomeDir {params.outprefix} \
    --genomeFastaFiles {input}
    """

rule filter_tRNA_rRNA: 
  input: 'trim_reads/{sample}/{fastq}','tRNA_rRNA_index/tRNA_rRNA_index.done'  
  output: 'filter_reads/{sample}/{fastq}'
  conda: '../envs/star'
  threads: 8
  params:
    indexname = lambda wc,input: input[1].replace('.done',''),
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

folder2filter = 'cutadapt_reads' if config.get('no_UMIs',False) else 'trim_reads'

#this rule is the 'signal spliter where we go from sample to indiv fastqs
def choose_processed_reads(wc,config=config):
  #correct zcat strings based on read pairs
  filedf = (seqfilesdf.loc[[wc['sample']]])
  isrna = not 'ribo' in sampledf.loc[wc['sample'],'assay'] #this should be made capital insensitive
  filedf = filedf[filedf.file_id==wc.fileid]
  #import ipdb; ipdb.set_trace()
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
  
rule makeGTF:
  input: GTF=GTF_orig
  output: GTF
  conda: '../envs/gffread'
  #conda: '~/miniconda3/envs/seq/bin/gffread'
  shell: r""" 
      # set -x
      #with filtering output all sequences
      zless {input.GTF} \
      | {mod_id_sed_cmd} \
      | gffread -F -T -o {GTF}

    """
 
rule make_utrs:
  input: GTF=GTF_orig
  output: fputrs='fputrs.gtf',tputrs='tputrs.gtf'
  # script: 'make_utrfiles.R'
  run:
    shell(r"""
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

     
      """) 


################################################################################
########Quality checking
################################################################################
  
def get_sample_fastqs(wc,mate='1',folder='processed_reads',seqfilesdf=seqfilesdf):
   #correct zcat strings based on read pairs
  filedf = (seqfilesdf.loc[[wc['sample']]])
  filedf = filedf.loc[filedf.mate==mate,]
  matefiles = [folder+'/'+wc['sample']+'/'+f for f in filedf.file_id]
  return(matefiles)

get_sample_fastqs2 = partial(get_sample_fastqs,mate='2')

#TODO
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

#note that it soft clips by default
#k -number of multimaps reported  

################################################################################
########GSNAP
################################################################################

#make a bedgraph with 1 for regions that have NO coverage
rule mask_bg:
  input: ALIGNER_TO_USE+'/data/{sample}/{sample}.bam'
  output: ALIGNER_TO_USE+'/data/{sample}/{sample}.{strand}.nocov.bg'
  shell: r"""
     bamCoverage -b {input} -o  {output}.tmp -bs 200  -of  bedgraph --normalizeUsing None -v
     cat {output}.tmp | grep -w "0$" | sed s'/0$/1/g' > {output} 
     rm {output}.tmp
   """

rule masked_fasta:
  input:  bgs=expand(ALIGNER_TO_USE+'/data/{sample}/{sample}.{strand}.nocov.bg',sample=rnasamples,strand=['forward','reverse']),REF=REF
  # input:  REF=REF,RNABAMs=ALIGNER_TO_USE+'/data/{sample}/{sample}.bam'
  # output: 'masked_fasta/{cell_line}/masked.fa'
  output: 'masked_fasta/masked.fa','masked_fasta/mask.bedgraph'
  params: bgs = lambda wc,input: '-i ' + ' -i '.join(input.bgs)
  shell: r"""
    bedtools unionbedg {params.bgs} | grep -vwe "0" > masked_fasta/mask.bedgraph
    bedtools maskfasta -fi {input.REF} -bed masked_fasta/mask.bedgraph -fo {output}
  """

# rule chrname_vcfs:
#   input: lambda wc: VCF_DICT[wc.celltype]
#   output: 'vcf/{celltype}/{celltype}.vcf'
#   shell: r""" 

#       awk '{{ 
#         if($0 !~ /^#/) 
#             print "chr"$0;
#         else if(match($0,/(##contig=<ID=)(.*)/,m))
#             print m[1]"chr"m[2];
#         else print $0 
#       }}' {input} > {output}

#       """

################################################################################
########STAR
################################################################################

rule star_index:
 input: REF=ancient(REF),GTF=ancient(GTF)
 output: touch('starindex/starindex.done')
 threads: 8
 conda: "../envs/star"
 shell: r"""
   STAR \
     --runThreadN {threads} \
     --runMode genomeGenerate \
     --genomeDir $(dirname {output}) \
     --sjdbGTFfile {input.GTF} \
     --genomeFastaFiles {input.REF}
   """

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

#changing from working with run/shell to only shell and putting everything in params. now have to

rule star:
     input:
          lfastqs=get_sample_fastqs,
          rfastqs=get_sample_fastqs2,
          STARINDEX='starindex/starindex.done',
     output:
          done = touch('star/data/{sample,[^/]+}/.done'),bam='star/data/{sample}/{sample}.bam',bai='star/data/{sample}/{sample}.bam.bai'
     threads: 8
     conda: "../envs/star"
     #please check params, there is a lot of stuff i did.
     params:
        sample= lambda wc: wc['sample'],
        #filestring = lambda wc: get_file_string(wc,seqfilesdf),
        GEN_DIR=lambda wc,input: input.STARINDEX.replace('starindex.done',''),
        #only used for remap (now remap='')
        markdup = lambda wc: '' if sampledf.assay[wc['sample']] == 'ribo' else '-m',
        platform = 'NotSpecified',
        outputdir = lambda wc,output: os.path.dirname(output[0]),
        repdir = lambda wc,output: os.path.dirname(output[0]).replace('data','reports'),
        #tophatindex =lambda wc,input: input['bowtie_index'].replace('.done',''),
        halfthreads = lambda wc,threads: threads/2,
        sortmem = lambda wc,threads: str(int(5000/(threads/2)))+'M',
        #remap = '1' if sampledf.assay[wildcards['sample']] == 'ribo' else ''
        remap = '',
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
########Hisat
################################################################################


#prepare the junctions for hisat to align with. These should be in files labelled per cell line
rule hisat_splice_files:
  # input: get_junction_files
  input: GTF
  # output: 'hisat_splice_files/{cell_line}/hisat_junctions.tsv'
  output: 'hisat_splice_files/hisat_junctions.tsv'
  conda: '../envs/hisat2'
  shell: r""" 
  mkdir -p $(dirname {output} )
  rm -f {output}
  hisat2_extract_splice_sites.py {GTF} >> {output} 

  """




def gethisatfasta(wc):
  return 'masked_fasta/masked.fa' if (wc['ismasked'] == 'masked') else str(REF)

#for a given individual, fetch the vcf file, and the appropriate splice sites, 
#and then construct a specific index
rule hisat_index:
  input: 
    # juncfile='hisat_splice_files/{cell_line}/hisat_junctions.tsv',
    juncfile= ['hisat_splice_files/hisat_junctions.tsv'],
    # REF=lambda wc: 'masked_fasta/{cell_line}/masked.fa' if sampledf[wildcards['cell_line']],
    # vcf='vcf/{celltype}/{celltype}.vcf',
    # vcfs = expand('vcf/{celltype}/{celltype}.vcf',celltype=VCF_DICT.keys()),
    REF  = ancient(gethisatfasta),
  # output: dir('hisat_index/{cell_line}/{rnamasked}/')
  output: 'hisat_index/{ismasked}'
  params: 
    # snps = lambda wc, input:  ','.join(input.vcfs),
    juncs = lambda wc,input: '--ss ' +(' --snp '.join(input.juncfile))
  #so 20x10 = 200GB of memory total
  threads: 20
  conda: '../envs/hisat2'
  shell:r"""
       # set -xe
      mkdir -p  {output}
      hisat2-build --wrapper basic-0 -p 20  --ss {input.juncfile} {input.REF} {output}
       """

       # if [-z "{snps}"] hisat2_extract_snps_haplotypes_VCF.py {input.REF} {params.snps} {output[0]}_snpfile
       # hisat2-build --wrapper basic-0 -p 20 --snp {output[0]}_snpfile.snp --hap {output[0]}_snpfile.haplotype --ss {input.juncfile} {input.REF} {output}



#hsat index needs to be masked for the riboseq, but not for the rnaseq
def get_hsat_index(wc):
  masked = sampledf.loc[wc['sample'],'assay'] == 'ribo'
  index = 'hisat_index/masked/' if masked else 'hisat_index/notmasked/' 
  return(index)

def get_hisat_strandopt(wc):
  optstring = ' --rna-strandness '

  paired = sampledf.loc[wc['sample'],'library_layout']=='PAIRED'
  strandedness = sampledf.loc[wc['sample'],'protocol']
  revstrand = strandedness=='reverse'
  
  if strandedness=='no':
    strandstring= ' '
  elif (not paired) & (not revstrand):
    strandstring=optstring+'F'
  elif (not paired) & (revstrand):
    strandstring=optstring+'R'
  elif (paired) & (not revstrand):
    strandstring=optstring+'FR'
  elif (paired) & (revstrand):
    strandstring=optstring+'RF'

  return(strandstring)

rule hisat:
     input:
          lfastqs=get_sample_fastqs,
          rfastqs=get_sample_fastqs2,
          # hisatindex=get_hsat_index,
          hisatindex='hisat_index/notmasked'
     output:
          'hisat2/data/{sample}/{sample}.bam'
     threads: 8
     conda: "../envs/hisat2"
     #please check params, there is a lot of stuff i did.
     params:
        sample= lambda wc: wc['sample'],
        lfastqlist = lambda wc,input: '-U ' + (','.join(input.lfastqs)) if len(input.rfastqs)==0 else '-1 ' + (','.join(input.lfastqs)),
        rfastqlist = lambda wc,input: '-2 ' + (','.join(input.rfastqs)) if len(input.rfastqs)!=0 else '',
        #only used for remap (now remap='')
        #markdup = lambda wc: '' if sampledf.assay[wildcards['sample']] == 'ribo' else '-m'
        strandstring = get_hisat_strandopt,
        platform = 'NotSpecified',
        outputdir = lambda wc,output: os.path.dirname(output[0]),
        repdir = lambda wc,output: os.path.dirname(output[0]).replace('data','reports'),
        #tophatindex =lambda wc,input: input['bowtie_index'].replace('.done',''),
        halfthreads = lambda wc,threads: threads/2,
        sortmem = lambda wc,threads: str(int(5000/(threads/2)))+'M',
        noalign = lambda wc,output: output[0].replace('.bam','noalign.fastq.gz')
        
     shell: r"""

        mkdir -p $(dirname {output})
        mkdir -p  {params.repdir}
        hisat2 \
          --no-softclip \
          -k 20 \
          {params.strandstring} \
          --un-gz {params.noalign} \
          --summary-file {params.repdir}/hisat_summary.txt \
          --threads 8  \
          -x {input.hisatindex}/ \
          {params.lfastqlist} {params.rfastqlist} \
          -S {output[0]}.tmp

        mv {output[0]}.tmp {output[0]}

       MY_TMP_DIR=$(mktemp -d)

        samtools sort \
        -m {params.sortmem} \
        -@ {params.halfthreads}\
        -T ${{MY_TMP_DIR}} \
        -o {output[0]}.sort \
        {output[0]}


        mv {output[0]}.sort {output[0]}

        samtools index {output[0]}

        mkdir -p {params.repdir}/{wildcards.sample}
        samtools stats {output[0]} > {params.repdir}/{wildcards.sample}/{wildcards.sample}.bamstats.txt
        samtools flagstat {output[0]} > {params.repdir}/{wildcards.sample}/{wildcards.sample}.flagstat.log
        samtools idxstats {output[0]} > {params.repdir}/{wildcards.sample}/{wildcards.sample}.idxstats.log

          """


################################################################################

  
rrna_intervals = 'qc/picard_rrna_intervals.txt'
refflat = snakedir/ 'qc' / Path(GTF).with_suffix('.refflat').name
#refflat = snakedir/ 'qc' / Path(config['GFF_orig']).with_suffix('.refflat').name

rule make_picard_files:
  input: GTF=GTF,bam=ALIGNER_TO_USE+'/data/'+samples[0]+'/'+samples[0]+'.bam'
  output: intervals=rrna_intervals,refflat=refflat
  #conda: '../envs/picard'
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
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################


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
     conda: '../envs/picard'
     resources:
     params:
        singleendflag = lambda wc: ' -e ' if sampledf.loc[wc['sample'],'library_layout'] != 'PAIRED' else '',
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

##### -o ${{OUTDIR}}/read_alignment_report.tsv not created
 
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################


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
      expand("feature_counts/data/{sample}/feature_counts", sample = samples),


      # 'sample_file.txt'
  #conda: '../envs/multiqc'
  params: 
    multiqcscript = config['multiqcscript'],
    sample_reads_file=config['sample_files'],
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




#this is going to count reads in each library over 5'UTRS, CDS, and 3' UTRs
rule readlenfilt:
  input: ALIGNER_TO_USE+'/data/{sample}/{sample}.bam'
  output:  'readlenfilt/data/{sample}/{readrange}/readlenfilt.bam'
  threads:4
  run:
    minreadlen,maxreadlen = wildcards['readrange'].split('_')
    readrangebam = output[0]
    shell(r"""
          set -ex
          samtools view -h $(dirname {input})/{wildcards.sample}.bam \
          | awk '((length($10) >= {minreadlen})&&(length($10) <= {maxreadlen})) || $1 ~ /^@/' \
          | samtools view -F 4 -S -b - > {readrangebam}.tmp
          #

         MY_TMP_DIR=$(mktemp -d)

         samtools sort \
          -@ {threads}\
          -m 2G \
          -T ${{MY_TMP_DIR}} \
          -o {readrangebam} \
          {readrangebam}.tmp

          samtools index {readrangebam}
      """)
       

rule feature_counts_readrange:
     input:
          'tputrs.gtf','fputrs.gtf',bam = 'readlenfilt/data/{sample}/{readrange}/readlenfilt.bam',GTF=GTF
     output:
          'feature_counts_readrange/data/{sample,[^/]+}/{gcol_generegion}/{readrange}/feature_counts'
     threads: 2
     log: r"""feature_counts_readrange/reports/{sample}/{gcol_generegion}/{readrange}/feature_counts.log"""
     run:
          groupcol,generegions = wildcards['gcol_generegion'].split('__')
          minreadlen,maxreadlen = wildcards['readrange'].split('_')
          if (generegions in ['CDS']): GTF,featuretype = input.GTF,'CDS'            
          elif (generegions in ['tputrs']): GTF,featuretype = tputrs,'exon'
          elif (generegions in ['fputrs']): GTF,featuretype = fputrs,'exon'
          else: GTF = GTF

          
          protocol = sampledf.loc[wildcards['sample'],'protocol']
          if (protocol == 'no'):
               protocol = 0
          elif (protocol == 'yes'):
               protocol = 1
          elif (protocol == 'reverse'):
               protocol = 2
          else:
               sys.exit('Protocol not known!')

          library = sampledf.loc[wildcards['sample'],'library_layout']

          if (library == 'PAIRED'):
               library = '-p'
          else:
               library = ''

          countmultimappers = ' ' 
          
          if (generegions=='tRNAs'):
            featuretype = 'tRNA'
            countmultimappers = '-M --fraction'
          

          sample = wildcards['sample']
          
          shell(r"""
          set -ex
          mkdir -p feature_counts_readrange/data/{sample}/{generegions}/{wildcards.readrange}/
          mkdir -p feature_counts/reports/{wildcards.sample}/
          
          featureCounts \
            -O \
            -T {threads} \
            -t {featuretype} -g {groupcol} \
            -a {GTF} \
            -s {protocol} {library} {countmultimappers} \
            -o {output} \
            {input.bam}             

            # &> feature_counts/reports/{wildcards.sample}/{wildcards.sample}.feature_counts.log

          """)

def get_readrange(wc):
    if wc['sample'] in ribosamples:
          selregion='gene_id__CDS'
          selreadrange='25_31'
    else:
          selregion='gene_id__CDS'
          selreadrange='20_1000'
    fcountfile = 'feature_counts_readrange/data/'+wc['sample']+'/'+selregion+'/'+selreadrange+'/feature_counts'
    return   fcountfile

rule feature_counts:
     input:
          get_readrange,
          ALIGNER_TO_USE+'/data/{sample}/{sample}.bam',GTF,'tputrs.gtf','fputrs.gtf',
     output:
          'feature_counts/data/{sample,[^/]+}/feature_counts'
     threads: 2
     run:       
        bamfile = ALIGNER_TO_USE+'/data/'+wildcards['sample']+'/'+wildcards['sample']+'.bam' 

        shell(r"""

          mkdir -p feature_counts/reports/{wildcards.sample}
          mkdir -p feature_counts/data/{wildcards.sample}

          head -n2 {input[0]} \
          | tail -n1 \
          | awk -v FS="\t" -v OFS="\t"  '{{$7= "{bamfile}";print $0}}' \
          >    feature_counts/data/{wildcards.sample}/feature_counts

          tail -n+3 {input[0]} \
          >>   feature_counts/data/{wildcards.sample}/feature_counts

          cp {input[0]}.summary {output[0]}.summary
          """)

# #this is going to count reads in each library over 5'UTRS, CDS, and 3' UTRs
rule aggregate_feature_counts:
  input : expand("feature_counts/data/{sample}/feature_counts", sample = samples),
  output: 'feature_counts/all_feature_counts'
  run:
    fcountfiles = list(input)
    shell(r""" 

       #( (sed '2q;d' {fcountfiles[0]} | cut -f1 && tail -n+3 {fcountfiles[0]}| sort -k1 | cut -f1) > {output})
       tail -n+3 {fcountfiles[0]}| sort -k1 | cut -f1 > {output}
       
       #now for eahc fcount table, join it to the ids
       for fcountfile in $(echo {fcountfiles}); do

          tail -n+3 $fcountfile| sort -k1 | cut -f1,7 | join {output} - | sort -k1 > {output}tmp
          mv {output}tmp {output}
       
       done

      echo "feature_id {samples}" | cat - {output} > {output}tmp
      mv {output}tmp {output}
    
      """)
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################


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
    R -e 'library("RiboseQC"); prepare_annotation_files(annotation_directory=".",gtf_file="{params.gtfmatchchrs}",annotation_name="{params.annobase}",forge_BS=FALSE, genome_seq=FaFile("{REF}"))'
 """

#R -e 'if (! "RiboseQC" %in% installed.packages()) devtools::install("{RIBOSEQCPACKAGE}",upgrade="never")'
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################


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
         R -e 'library("RiboseQC");RiboseQC_analysis("{params.annofile}", bam="{input.bam}",rescue_all_rls=TRUE,dest_names="{params.outname}", genome_seq = "{REF}", report_file="{params.report_file}")'
     """

rule segment_periodicity:
  input: '../ext_data/segments.gtf',GTF,'riboseqc/data/{sample}/_for_ORFquant',ALIGNER_TO_USE+'/data/{sample}/{sample}.bam'
  output: 'segment_scores/{sample}/segment_scores.tsv'
  shell: r''' Rscript ../src/segment_periodicity.R {input} {output} '''

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
      R -e 'library("ORFquant");run_ORFquant(for_ORFquant_file = {params.for_ORFquantfile},annotation_file = "{input.annofile}", n_cores = {threads},prefix="{params.outputdir}") '
        
      """
######## Here loading ORFquant can be done directly within R after installing it.
###       R -e 'devtools::load_all("{ORFquantPACKAGE}");run_ORFquant(for_ORFquant_file = {params.for_ORFquantfile},annotation_file = "{params.annofile}",genome_seq = "{REF}", n_cores = {threads},prefix="{params.outputdir}") '
###      R -e 'if (is.element("ORFquant",installed.packages()[,1]) == 0) {devtools::install("{ORFquantPACKAGE}")}'
################################################################################
########Mappaability
################################################################################


rule mappability_reads:
  input: REF,RNAFASTA
  output: 'mappability_reads/mappability_{kmer}/mappability_{kmer}.fastq.gz'
  threads: 8
  shell: r"""

  samtools faidx {RNAFASTA}
  cut -f1,2 {RNAFASTA}.fai  > {output}.trsizes

  #do this for entire genome
  set -x
  #these one liners, in order, cat the chromosome sizesinto awk to create bed files spanning the chromosomes, 
  #create fastas from a bed file, stick together every second line, generate the tiles of kmer size, modify the names of hte
  #tiles, and then finally turn them into fastq format

  cat {output}.trsizes \
    | awk '{{print $1"\t"1"\t"$2}}' \
    | bedtools getfasta -s -fi {RNAFASTA} -bed -  \
    | sed '$!N;s/\n/ /' \
    | perl -lane '$i=0;while($i<=(length($F[1])-{wildcards.kmer})){{print $F[0] , "$i\n", substr $F[1],$i,{wildcards.kmer} ; $i = $i +1}}' \
    | perl -lan -F'[\:\-|\)|\(]' -e 'if( /^>/){{print $F[0],"_",$F[1]+$F[4],"_",$F[1]+$F[4]+{wildcards.kmer}}}else{{print @F}}' \
    | perl -lanpe 's/>/@/ ; s/^([^@]+)/\1\n+/; if($1){{$a="I" x length($1); s/\+/+\n$a/}}' \
    | gzip > {output}

    """

################################################################################
########Run Rseq
################################################################################
rule rseq:
  input: '../src/rseq_design.yaml','../src/build_project.R',expand('feature_counts/data/{sample}/feature_counts{summary}',sample=samples,summary=['','.summary'])
  output: directory('run_rseq/subanalyses')
  shell: r"""
  set -ex
  mkdir -p {output[0]}

  cd ..

  Rscript src/build_project.R

  """
