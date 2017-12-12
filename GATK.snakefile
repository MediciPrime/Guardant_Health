configfile:  "config.yaml"

targets = []
for sample in config['samples']:
    targets.extend(
        expand(
            'data/variants/{reference}/{sample}/{sample}.g.vcf',
            sample=sample,
            reference=config['samples'][sample]['reference']
        )
    )

rule all:
    input:
        targets

rule cutadapt:
    input:
        r1 = 'data/raw_reads/{sample}_R1.fastq.gz',
        r2 = 'data/raw_reads/{sample}_R2.fastq.gz'
    output:
        r1 = 'data/trimmed_reads/{sample}_R1_trimmed.fq',
        r2 = 'data/trimmed_reads/{sample}_R2_trimmed.fq'
    params:
        five_prime = 'GTTCAGAGTTCTACAGTCCGACGATC',
        three_prime = 'TGGAATTCTCGGGTGCCAAGG'
    shell:
        'cutadapt -a {params.five_prime} -A {params.three_prime} '
        '-q 20 -m 25 -o {output.r1} -p {output.r2} '
        '{input.r1} {input.r2}'

rule bwa_index:
    input:
        ref_fasta = 'data/reference/{reference}/{reference}.fasta'
    output:
        'data/reference/{reference}/{reference}.fasta.bwt'
    shell:
        'bwa index {input.ref_fasta}'

rule bwa_align:
    input:
        ref_fasta = 'data/reference/{reference}/{reference}.fasta',
        r1 = 'data/trimmed_reads/{sample}_R1_trimmed.fq',
        r2 = 'data/trimmed_reads/{sample}_R2_trimmed.fq'
    output:
        'data/sam_files/{reference}/{sample}/{sample}.sam'
    params:
        sample = config['references'][reference]['sample']
    threads: 10
    shell:
        'bwa mem {input.ref_fasta} {input.r1} {input.r2} > {output}'

def get_bar(wc):
    barcode = config['samples'][wc.sample]['barcode']
    return barcode
        
rule bwa_picard_rg:
    input:
        sam_file = 'data/sam_files/{reference}/{sample}/{sample}.sam'
    output:
        'data/bam_files/{reference}/{sample}/{sample}_rg_sorted.bam'
    params:
        barcode=get_bar
    shell:
        'picard AddOrReplaceReadGroups I={input.sam_file} O={output} '
        'SORT_ORDER=coordinate RGLB={sample} RGPL=illumina RGPU={params.barcode} '
        'RGSM={wildcards.sample}'
        
rule bwa_picard_md:
    input:
        bam_sorted = 'data/bam_files/{reference}/{sample}/{sample}_rg_sorted.bam'
    output:
        bam_dedup = 'data/bam_files/{reference}/{sample}/{sample}_rg_dedupped.bam',
        output_met = 'data/bam_files/{reference}/{sample}/{sample}_output.metrics'
    shell:
        'picard MarkDuplicates I={input.bam_sorted} O={output.bam_dedup} '
        'CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M={output.output_met}'

rule fasta_index:
    input:
        ref = 'data/reference/{reference}/{reference}.fasta'
    output:
        ref_index = 'data/reference/{reference}/{reference}.fasta.fai'
    shell:
        'samtools faidx {input.ref}'

rule fasta_dictionary:
    input:
        ref = 'data/reference/{reference}/{reference}.fasta'
    output:
        ref_dict = 'data/reference/{reference}/{reference}.dict'
    shell:
        'picard CreateSequenceDictionary R={input.ref}'      
        
rule gatk_haplotype:
    input:
        ref = 'data/reference/{reference}/{reference}.fasta',
        bam_split = 'data/bam_files/{reference}/{sample}/{sample}_rg_dedupped.bam',
        ref_index = 'data/reference/{reference}/{reference}.fasta.fai',
        ref_dict = 'data/reference/{reference}/{reference}.dict'
    output:
        'data/variants/{reference}/{sample}/{sample}.g.vcf'
    shell:
        'java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R {input.ref} '
        '-I {input.bam_split} --genotyping_mode DISCOVERY -stand_call_conf 30 '
        '-stand_emit_conf 10 -o {output}'
