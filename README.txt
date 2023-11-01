TCS: Tiled-ClickSeq batch scripts
Last Modified: Oct-23


Tiled Batch
        Written by Andrew Routh 2020-2023

USAGE: ./TCS_batch [OPTIONS] FILE PRIMER

e.g. ./TCS_batch /path/to/data/mydata_R1.fastq Virus_primers.BED

Required Arguments:
    File        Enter full path of R1 file
    Primer      (BED file of primers)

Optional Arguments:
    -h show this help text

    -p Perform custom stages; select combination of P, M, R, C, D. No whitespace, e.g. 'PM'.
            (default = PM: Preprocess and map)
        P Only perform data preprocessing
        M Only default data mapping
        Q Map all Split files to find primer targets
        S Perform Strict data processing, mapping and reconstruction
        V Do ViReMa mapping
        T Do Minority variant calling using PARCL script
        C Clean/Remove temporary/previous files and directories (SplitFASTQ folder and bowtie2 mappings to pre-concensus genome)

    -t set threads (default: 1)

    -g provide base genome

    -m provide metadata file with root vs genome names




Requirements:

bowtie
bowtie2
samtools
fastp
pilon
fastx_toolkit
cutadapt
pypy
bedtools	

Optional:
umi_tools
hisat2



