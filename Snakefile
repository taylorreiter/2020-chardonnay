
SAMPLES = ['BM19_64_6_UMI_S75', 'BM19_65_4_UMI_S76', 'BM19_65_6_UMI_S77', 
           'BM19_67_4_UMI_S78', 'BM19_67_5_UMI_S79', 'BM19_67_6_UMI_S80', 
           'BM19_68_4_UMI_S81', 'BM19_68_5_UMI_S82', 'BM19_68_6_UMI_S83', 
           'BM19_69_1_UMI_S113', 'BM19_69_4_UMI_S84', 'BM19_69_5_UMI_S85',
           'BM19_69_6_UMI_S86', 'BM19_70_4_UMI_S87', 'BM19_70_5_UMI_S88', 
           'BM19_70_6_UMI_S89', 'BM19_71_4_UMI_S114', 'BM19_71_5_UMI_S90',
           'BM19_71_6_UMI_S91', 'BM19_73_1_UMI_S92', 'BM19_73_2_UMI_S93', 
           'BM19_73_3_UMI_S94', 'BM19_73_4_UMI_S95', 'BM19_73_5_UMI_S96', 
           'BM19_73_6_UMI_S97', 'BM19_74_1_UMI_S98', 'BM19_74_2_UMI_S99', 
           'BM19_74_3_UMI_S100', 'BM19_74_4_UMI_S101', 'BMCHAR_IN_UMI_S112']

rule all:
   input:
       expand('outputs/dedup/{sample}.dedup.bam', sample = SAMPLES)

rule umi_extract:
    input: "inputs/raw/{sample}_L003_R1_001.fastq.gz"
    output: "outputs/umi_extract/{sample}.fq.gz"
    conda: "envs/umitools.yml"
    shell:'''
    umi_tools extract --stdin={input} --bc-pattern=NNNNNN --log={output}.log --stdout {output}
    '''

rule fastqc:
    input: "outputs/umi_extract/{sample}.fq.gz"
    output: "outputs/fastqc_umi/{sample}.html"
    params: outdir = "outputs/fastqc_umi"
    conda: "envs/fastqc.yml"
    shell:'''
    fastqc -o {params.outdir} -t 8 --nogroup {input}
    '''

rule bbduk_trim:
    input: 
        reads = "outputs/umi_extract/{sample}.fq.gz",
        polya = "inputs/polya.fa",
        adapters = "inputs/truseq_rna.fa.gz"
    output: "outputs/bbduk/{sample}.fq"
    conda: "envs/bbmap.yml"
    shell:'''
    bbduk.sh in={input.reads} out={output} ref={input.polya},{input.adapters} \
    k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20
    '''

rule fastqc_trim:
    input: "outputs/bbduk/{sample}.fq"
    output: "outputs/fastqc_trim/{sample}.html"
    params: outdir = "outputs/fastqc_trim"
    conda: "envs/fastqc.yml"
    shell:'''
    fastqc -o {params.outdir} -t 8 --nogroup {input}
    '''

#rule multiqc_trim:


rule download_genome:
    output: 'inputs/genome/GCF_000146045.2_R64_genomic.fna.gz'
    shell:'''
    wget -O {output} https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Saccharomyces_cerevisiae/latest_assembly_versions/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz
    '''

rule decompress_genome:
    input: 'inputs/genome/GCF_000146045.2_R64_genomic.fna.gz'
    output: 'inputs/genome/GCF_000146045.2_R64_genomic.fna'
    shell: "gunzip {input}"

rule download_gtf:
    output: 'inputs/genome/GCF_000146045.2_R64_genomic.gtf.gz'
    shell:'''
    wget -O {output} https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Saccharomyces_cerevisiae/latest_assembly_versions/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gtf.gz
    '''

rule decompress_gtf:
    input: 'inputs/genome/GCF_000146045.2_R64_genomic.gtf.gz'
    output: 'inputs/genome/GCF_000146045.2_R64_genomic.gtf'
    shell: "gunzip {input}"

rule star_index_genome:
    input:
        genome = 'inputs/genome/GCF_000146045.2_R64_genomic.fna',
        gtf = 'inputs/genome/GCF_000146045.2_R64_genomic.gtf'
    params: input_dir = 'inputs/genome' 
    output: 'inputs/genome/SAindex'
    shell:'''
    STAR --runThreadN 1 --runMode genomeGenerate --genomeDir {params.input_dir} \
         --genomeFastaFiles {input.genome} --sjdbGTFfile {input.gtf} --sjdbOverhang  99
    '''

rule star_align:
    #v2.5.2a
    input:
        reads = 'outputs/bbduk/{sample}.fq',
        genome_index = 'inputs/genome/SAindex'
    output: 'outputs/star/{sample}Aligned.sortedByCoord.out.bam' 
    params: 
        out_prefix = lambda wildcards: 'outputs/star/' + wildcards.sample,
        genome_dir = 'inputs/genome'
    conda: 'envs/star.yml'
    shell:'''
    STAR --runThreadN 2 --genomeDir {params.genome_dir}      \
        --readFilesIn {input.reads} --outFilterType BySJout  \
        --outFilterMultimapNmax 20 --alignSJoverhangMin 8    \
        --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.6 --alignIntronMin 20 \
        --alignIntronMax 1000000 --alignMatesGapMax 1000000  \
        --outSAMattributes NH HI NM MD --outSAMtype BAM      \
        SortedByCoordinate --outFileNamePrefix {params.out_prefix}
    '''

rule index_bam:
    input: 'outputs/star/{sample}Aligned.sortedByCoord.out.bam'
    output: 'outputs/star/{sample}Aligned.sortedByCoord.out.bam.bai'
    conda: 'envs/samtools.yml'
    shell:'''
    samtools index {input}
    '''
    
rule dedup_UMI:
    input: 
        bam = 'outputs/star/{sample}Aligned.sortedByCoord.out.bam',
        bai = 'outputs/star/{sample}Aligned.sortedByCoord.out.bam.bai'
    output: 'outputs/dedup/{sample}.dedup.bam'
    conda: 'envs/umitools.yml'
    shell:'''
    umi_tools dedup -I {input.bam} --output-stats=deduplicated -S {output}
    '''

#rule htseq_count:?
#    # v0.6.0
#    shell:'''
#    htseq-count -m intersection-nonempty -s yes -f bam -r pos $bam {input.gtf} > $bam_dir/read_counts.txt
#    '''


#rule deseq2
