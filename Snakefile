import glob
import os
import errno
from itertools import chain
from os.path import join
import glob
import re
from collections import defaultdict
import numpy as np
import pandas as pd

include: config['config_path']
import pandas as pd
df = pd.read_table(TO_MERGE_TSV)
df = df.sort_values(by=['replicate', 'condition', 'assay'])
#print(df.head())
df_grouped = df.groupby(['assay', 'condition', 'replicate'])
group_keys = df_grouped.groups.keys()

ASSAYS = []
CONDITIONS = []
REPLICATES = []
SAMPLES = []
for x in group_keys:
    ASSAYS.append(x[0])
    CONDITIONS.append(x[1])
    REPLICATES.append(int(x[2]))
    SAMPLES.append('{}-{}-{}'.format(x[1], x[0], x[2]))

print(SAMPLES)
RIBO_SAMPLES = [sample for sample in SAMPLES if 'ribo' in sample]
RNA_SAMPLES = [sample for sample in SAMPLES if 'rna' in sample]

workdir: OUT_DIR

def merge_fastq_input(wildcards):
    assay = wildcards.assay
    condition = wildcards.condition
    replicate = int(wildcards.replicate)
    sample_list = df_grouped.get_group((assay, condition, replicate)).sampleID.tolist()
    return [os.path.join(RAWDATA_DIR, '') + '{}.fastq.gz'.format(sample) for sample in sample_list]

#RIBOTRICER_ANNOTATION = '/home/cmb-panasas2/skchoudh/github_projects/C_albicans_project/notebooks/June_06_GTF_analysis/C_albicans_SC5314_version_A22-s07-m01-r27_features.UTR5_CDS_UTR3_cleaned_ribotricer_longest_candidate_orfs.tsv'
RIBOTRICER_ANNOTATION = '/home/cmb-panasas2/skchoudh/github_projects/C_albicans_project/notebooks/July_07_GTF_analysis/C_albicans_SC5314_version_A22-s07-m01-r27_features.UTR5_CDS_UTR3_cleaned_with_noveltx_ribotricer_longest_candidate_orfs.tsv'
RIBOTRICER_ANNOTATION_SEQ = '/home/cmb-panasas2/skchoudh/github_projects/C_albicans_project/notebooks/July_07_GTF_analysis/C_albicans_SC5314_version_A22-s07-m01-r27_features.UTR5_CDS_UTR3_cleaned_with_noveltx_ribotricer_longest_candidate_orfs_seq.tsv'
REGIONS = ['super_uORF', 'annotated', 'super_dORF']
GENE_BED = '/home/cmb-panasas2/skchoudh/genomes/C_albicans_SC5314/Assembly22' + '/annotation/' + '/C_albicans_SC5314_version_A22-s07-m01-r27_features.bed12'
UTR5_BED = '/home/cmb-panasas2/skchoudh/github_projects/C_albicans_project/Bruno_et_al_2014_data/Bruno_et_al_UTR5.A21toA22_hapAB_merged.collapsed.named_uniquehits.bed'
UTR3_BED = '/home/cmb-panasas2/skchoudh/github_projects/C_albicans_project/Bruno_et_al_2014_data/Bruno_et_al_UTR3.A21toA22_hapAB_merged.collapsed.named_uniquehits.bed'
CDS_BED = '/home/cmb-panasas2/skchoudh/genomes/C_albicans_SC5314/Assembly22' + '/annotation/' + '/C_albicans_SC5314_version_A22-s07-m01-r27_features.gtfTogenePred.genePredToBed6'
SRC_DIR = '/home/cmb-panasas2/skchoudh/github_projects/C_albicans_project/scripts'
GTF = '/home/cmb-panasas2/skchoudh/github_projects/C_albicans_project/notebooks/June_06_GTF_analysis/C_albicans_SC5314_version_A22-s07-m01-r27_features.UTR5_CDS_UTR3_cleaned.gtf'
RRNA_FASTA = '/home/cmb-panasas2/skchoudh/genomes/C_albicans_SC5314/Assembly22/fasta/C_albicans_SC5314_version_A22-s07-m01-r27_other_features_no_introns.fasta'
#STAR_INDEX = '/home/cmb-panasas2/skchoudh/genomes/C_albicans_SC5314/Assembly22/star_mrna_index_UTR5_CDS_UTR3_cleaned_gtf_ribopod'
STAR_INDEX = '/home/cmb-panasas2/skchoudh/genomes/C_albicans_SC5314/Assembly22/star_mrna_index_r27_features_gtf'
RRNA_INDEX = '/home/cmb-panasas2/skchoudh/genomes/C_albicans_SC5314/Assembly22/star_rrna_index_ribopod'


# Suffix of all fastq
COMMON_SUFFIX = '_R1_001'
REGIONS = ['super_uORF,uORF', 'super_uORF,uORF,overlap_uORF', 'super_uORF,overlap_uORF', 'annotated', 'super_dORF', 'uORF']
THREADS = 16 
# All this data is forward stranded

def mkdir_p(path):
  """Python version mkdir -p

  Parameters
  ----------

  path : str
  """
  if path:
    try:
      os.makedirs(path) 
    except OSError as exc:  # Python >2.5
      if exc.errno == errno.EEXIST and os.path.isdir(path):
        pass
      else:
        raise

mkdir_p(os.path.join(OUT_DIR, 'slurm-logs'))

TOTAL_GENOME_SIZE = 28605418
SA_INDEX_Nbases = int(np.floor(min(14, np.log2(TOTAL_GENOME_SIZE) / 2.0 - 1)))


rule all:
    input:
      #STAR_INDEX + '/' + 'chrName.txt',
      #RRNA_INDEX + '/' + 'SA'
      expand('preprocessed/merged_fastq/{condition}-{assay}-{replicate}.fastq.gz', zip, condition=CONDITIONS, assay=ASSAYS, replicate=REPLICATES),
      expand('bams_unique/{sample}.bam', sample=SAMPLES),
      expand('ribotricer_rna_results/{sample}_translating_ORFs.tsv', sample=RNA_SAMPLES),
      expand('ribotricer_rna_results/region_counts/{region}/{sample}_counts_cnt.txt', region=REGIONS, sample=RNA_SAMPLES),
      expand('ribotricer_ribo_results/{sample}_translating_ORFs.tsv', sample=RIBO_SAMPLES),
      expand('ribotricer_ribo_results/region_counts/{region}/{sample}_counts_cnt.txt', region=REGIONS, sample=RIBO_SAMPLES),
      expand('pickled_data/counts/mRNA/{sample}.pickle', sample=SAMPLES),
      expand('pickled_data/counts/mRNA_collapsed/{sample}.pickle', sample=SAMPLES),
      expand('pickled_data/raw_fastq_reads/{sample}.txt', sample=SAMPLES),
      expand('pickled_data/trimmed_fastq_reads/{sample}.txt', sample=SAMPLES),
      expand('pickled_data/counts/rRNA/{sample}.pickle', sample=SAMPLES),
      expand('pickled_data/counts/tRNA/{sample}.pickle', sample=SAMPLES),
      expand('pickled_data/counts/snRNA/{sample}.pickle', sample=SAMPLES),
#expand('pickled_data/counts/snoRNA/{sample}.pickle', sample=SAMPLES),

rule create_index:
    input:
        GENOME_FASTA,
        GTF
    output: STAR_INDEX + '/' + 'chrName.txt'
    threads: THREADS
    shell:
        r'''mkdir -p {STAR_INDEX} && STAR --runThreadN {threads}\
            --genomeSAindexNbases {SA_INDEX_Nbases}\
            --runMode genomeGenerate\
            --genomeDir {STAR_INDEX}\
            --genomeFastaFiles {input[0]}\
            --sjdbGTFfeatureExon exon\
            --sjdbGTFfile {input[1]}'''

rule create_rrna_index:
    input:
        RRNA_FASTA
    output: RRNA_INDEX + '/' + 'SA'
    threads: THREADS
    shell:
        r'''mkdir -p {RRNA_INDEX} && STAR --runThreadN {threads}\
            --runMode genomeGenerate \
            --genomeDir {RRNA_INDEX} \
            --genomeFastaFiles {input}
        '''

rule merge_fastqs:
    input: merge_fastq_input
    output: 'preprocessed/merged_fastq/{condition}-{assay}-{replicate}.fastq.gz'
    shell: 
        r'''cat {input} > {output}'''


rule perfom_trimming:
    input:
        R1 = 'preprocessed/merged_fastq/{sample}.fastq.gz'
    params:
        out_dir='preprocessed/trimmed',
        phred_cutoff=5
    output:
        'preprocessed/trimmed/{sample}_trimmed.fq.gz'
    shell:
        '''
            trim_galore -o {params.out_dir} -q {params.phred_cutoff} {input.R1}
        '''

rule map_star_rRNA:
    input:
        R1='preprocessed/trimmed/{sample}_trimmed.fq.gz',
        index=RRNA_INDEX
    output:
        bam='mapped_rRNA/bams/{sample}.bam',
        fastq='unmapped_rRNA/fastq/{sample}.unmapped_rRNA.fastq.gz'
    params:
        prefix = 'mapped_rRNA/{sample}.map_rRNA',
        fastq_prefix='unmapped_rRNA/fastq/{sample}.unmapped_rRNA.fastq',
        unmapped = 'unmapped_rRNA/fastq',
        starlogs = 'mapped_rRNA/starlogs'
    threads: THREADS
    shell:
        r'''
        STAR --runThreadN {threads}\
             --genomeDir {input.index}\
             --outFilterMismatchNmax 2\
             --outFileNamePrefix {params.prefix} --readFilesIn {input.R1}\
             --outReadsUnmapped Fastx\
             --readFilesCommand zcat\
             --outSAMtype BAM Unsorted\
             --outTmpDir /tmp/{wildcards.sample}rRNA\
             && samtools sort -o {params.prefix}Aligned.sortedByCoord.out.bam -T /tmp/{wildcards.sample} {params.prefix}Aligned.out.bam && rm {params.prefix}Aligned.out.bam\
             && mv {params.prefix}Aligned.sortedByCoord.out.bam {output.bam} && mkdir -p {params.starlogs} && mv {params.prefix}Log.final.out {params.prefix}Log.out {params.prefix}Log.progress.out {params.starlogs}\
             && mkdir -p {params.unmapped} && mv {params.prefix}Unmapped.out.mate1 {params.fastq_prefix} && gzip {params.fastq_prefix} && samtools index {output.bam}
        '''


rule map_star_mRNA:
    input:
        R1='unmapped_rRNA/fastq/{sample}.unmapped_rRNA.fastq.gz',
        index=STAR_INDEX
    output:
        bam = 'mapped_mRNA/bams_star/{sample}.bam',
        fastq = 'unmapped_mRNA/fastq/{sample}.unmapped_mRNA.fastq.gz'
    params:
        prefix = 'mapped_mRNA/{sample}.map_mRNA',
        fastq_prefix = 'unmapped_mRNA/fastq/{sample}.unmapped_mRNA.fastq',
        unmapped = 'unmapped_mRNA/fastq',
        starlogs = 'mapped_mRNA/starlogs'
    threads: THREADS
    shell:
        r'''
        STAR --runThreadN {threads}\
        --genomeDir {input.index}\
        --outFilterMismatchNmax 2\
        --readFilesCommand zcat\
        --outFileNamePrefix {params.prefix} --readFilesIn {input.R1}\
        --alignEndsType EndToEnd\
        --outSAMtype BAM SortedByCoordinate\
        --limitBAMsortRAM 1352618371\
        --outReadsUnmapped Fastx &&\
        mv {params.prefix}Aligned.sortedByCoord.out.bam {output.bam} &&\
        mkdir -p {params.starlogs} &&\
        mv {params.prefix}Log.final.out {params.prefix}Log.out {params.prefix}Log.progress.out {params.prefix}SJ.out.tab {params.starlogs} &&\
        mkdir -p {params.unmapped} &&\
        mv {params.prefix}Unmapped.out.mate1 {params.fastq_prefix} &&\
        gzip {params.fastq_prefix} &&\
        samtools index {output.bam}
        '''


rule infer_experiment:
    input: 'mapped_mRNA/bams_star/{sample}.bam'
    output: 'inferred_experiment/{sample}.txt'
    shell:
        r'''
        riboraptor infer-protocol --bam {input} --refseq {GENE_BED} > {output}
        '''

rule filter_uniq_mapping:
    ## Based on https://groups.google.com/d/msg/rna-star/FhVk4f61Vcg/osx6DLTpAgAJ
    input: 'mapped_mRNA/bams_star/{sample}.bam'
    output: 'mapped_mRNA/bams_NH1/{sample}.NH1.bam'
    threads: THREADS
    shell:
        r'''
            bamtools filter -tag NH:1 -mapQuality ">=255" -mapQuality "<256" -in {input} -out {output}.unsorted && samtools sort -T /tmp/{wildcards.sample}_uniq_mapping -o {output} {output}.unsorted && samtools index {output}
        '''

        
rule filter_twice_mapping:
    ## Based on https://groups.google.com/d/msg/rna-star/FhVk4f61Vcg/osx6DLTpAgAJ
    input: 'mapped_mRNA/bams_star/{sample}.bam'
    output: 'mapped_mRNA/bams_NH2/{sample}.NH2.bam'
    threads: THREADS
    shell:
        r'''
            bamtools filter -tag NH:2 -mapQuality "<255" -in {input} -out {output}.unsorted && samtools sort -T /tmp/{wildcards.sample}_twice_mapping -o {output} {output}.unsorted && samtools index {output}
        '''

rule collapse_to_uniq_mapping:
    input: 'mapped_mRNA/bams_NH2/{sample}.NH2.bam'
    output: 'mapped_mRNA/bams_NH2ToNH1/{sample}.NH2ToNH1.bam'
    threads: THREADS
    shell:
        r'''samtools index {input} && python {SRC_DIR}/count_twice_map_uniq_reads.py {input} {output}.unsorted && samtools sort -T /tmp/{wildcards.sample}_collapse -o {output} {output}.unsorted && samtools index {output}
        '''

rule merge_NH2NH1_NH1:
    input: 
        NH2NH1='mapped_mRNA/bams_NH2ToNH1/{sample}.NH2ToNH1.bam',
        NH1='mapped_mRNA/bams_NH1/{sample}.NH1.bam'
    output: 'mapped_mRNA/bams_collapsed/{sample}.bam'
    threads: THREADS
    shell:
        r'''bamtools merge -in {input.NH2NH1} -in {input.NH1} -out {output}.unsorted && samtools sort -T /tmp/{wildcards.sample} -o {output} {output}.unsorted && samtools index {output}'''

rule extract_uniq_mapping:
    input: 'mapped_mRNA/bams_collapsed/{sample}.bam'
    output: 'bams_unique/{sample}.bam'
    params:
        tmp_dir = '/tmp'
    threads: 16
    shell:
        r'''
        samtools view -b -q 255 \
        {input} -o {output}.temp \
        && samtools sort {output}.temp -o {output} \
        -T {params.tmp_dir}/{wildcards.sample}_sort \
        && rm -rf {output}.temp \
        && samtools index {output}
        '''


rule sort_unique_bams:
    input: 'bams_unique/{sample}.bam'
    output: 'bams_sortedByCoord/{sample}.sortedByCoord.bam'
    shell:
        r'''samtools sort {input} -T /tmp/ -O bam -o {output}'''


rule add_xs:
    input: 'bams_sortedByCoord/{sample}.sortedByCoord.bam'
    output: 'bams_sortedByCoord/{sample}.sortedByCoord.XS.bam'
    shell:
        r'''samtools view -h {input} | awk -v strType=1 -f {SRC_DIR}/tagXSstrandedData.awk | samtools view -bS - > {output}'''

rule predict_orfs_rna:
  input:
    bam = 'bams_unique/{sample}.bam',
  output: 'ribotricer_rna_results/{sample}_translating_ORFs.tsv'
  params:
    prefix = 'ribotricer_rna_results/{sample}',
  shell:
    r'''ribotricer_rna detect-orfs --report_all --bam {input.bam} --prefix {params.prefix} --ribotricer_index {RIBOTRICER_ANNOTATION}'''

rule count_orfs_rna:
  input: 'ribotricer_rna_results/{sample}_translating_ORFs.tsv'
  output: 'ribotricer_rna_results/region_counts/{region}/{sample}_counts_cnt.txt'
  params:
    prefix = 'ribotricer_rna_results/region_counts/{region}/{sample}_counts',
    region = '{region}'
  shell:
    r'''ribotricer_rna count-orfs --report_all --ribotricer_index {RIBOTRICER_ANNOTATION} --features {params.region} --detected_orfs {input} --prefix {params.prefix}'''


rule predict_orfs:
    input: 'bams_unique/{sample}.bam'
    output: 'ribotricer_ribo_results/{sample}_translating_ORFs.tsv'
    params:
      prefix = 'ribotricer_ribo_results/{sample}',
    shell:
      r'''ribotricer detect-orfs --report_all --bam {input} --prefix {params.prefix} --ribotricer_index {RIBOTRICER_ANNOTATION}'''

rule count_orfs:
    input: 'ribotricer_ribo_results/{sample}_translating_ORFs.tsv'
    output: 'ribotricer_ribo_results/region_counts/{region}/{sample}_counts_cnt.txt'
    params:
      prefix = 'ribotricer_ribo_results/region_counts/{region}/{sample}_counts',
      region = '{region}'
    shell:
      r'''ribotricer count-orfs --report_all --ribotricer_index {RIBOTRICER_ANNOTATION} --features {params.region} --detected_orfs {input} --out {output}'''

rule count_rrna:
    input: 'mapped_rRNA/bams/{sample}.bam'
    output: 'pickled_data/counts/rRNA/{sample}.pickle'
    shell:
        r'''python /home/cmb-panasas2/skchoudh/github_projects/C_albicans_project/scripts/count_reads.py {input} rRNA {output}'''

rule count_trna:
    input: 'mapped_rRNA/bams/{sample}.bam'
    output: 'pickled_data/counts/tRNA/{sample}.pickle'
    shell:
        r'''python /home/cmb-panasas2/skchoudh/github_projects/C_albicans_project/scripts/count_reads.py {input} tRNA {output}'''

rule count_snrna:
    input: 'mapped_rRNA/bams/{sample}.bam'
    output: 'pickled_data/counts/snRNA/{sample}.pickle'
    shell:
        r'''python /home/cmb-panasas2/skchoudh/github_projects/C_albicans_project/scripts/count_reads.py {input} snRNA {output}'''

rule count_snorna:
    input: 'mapped_rRNA/bams/{sample}.bam'
    output: 'pickled_data/counts/snoRNA/{sample}.pickle'
    shell:
        r'''python /home/cmb-panasas2/skchoudh/github_projects/C_albicans_project/scripts/count_reads.py {input} snoRNA {output}'''

rule count_mrna:
    input: 'mapped_mRNA/bams_star/{sample}.bam'
    output: 'pickled_data/counts/mRNA/{sample}.pickle'
    shell:
        r'''python /home/cmb-panasas2/skchoudh/github_projects/C_albicans_project/scripts/count_reads.py {input} mRNA {output}'''

rule count_mrna_collaped:
    input: 'bams_unique/{sample}.bam'
    output: 'pickled_data/counts/mRNA_collapsed/{sample}.pickle'
    shell:
        r'''python /home/cmb-panasas2/skchoudh/github_projects/C_albicans_project/scripts/count_reads.py {input} mRNA {output}'''

rule count_fastq_reads:
    input: 'preprocessed/merged_fastq/{sample}.fastq.gz'
    output: 'pickled_data/raw_fastq_reads/{sample}.txt'
    shell:
        r'''zcat {input} | echo $((`wc -l`/4)) > {output}'''

rule count_trimmed_fastq_reads:
    input: 'preprocessed/trimmed/{sample}_trimmed.fq.gz',
    output: 'pickled_data/trimmed_fastq_reads/{sample}.txt'
    shell:
        r'''zcat {input} | echo $((`wc -l`/4)) > {output}'''
