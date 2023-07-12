#!/bin/bash

for d in v1.6 v1.6/GRCh37 v1.6/GRCh38
do
    if [ ! -d $d ]; then
        mkdir $d
    fi
done


# for GRCh37 / hg19
cd v1.6/GRCh37
wget -c https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh37/annotationsGRCh37_v1.6.tar.gz
wget -c https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh37/whole_genome_SNVs.tsv.gz
wget -c https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh37/whole_genome_SNVs.tsv.gz.tbi
wget -c https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh37/InDels.tsv.gz
wget -c https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh37/InDels.tsv.gz.tbi
tar -kzxvf annotationsGRCh37_v1.6.tar.gz

# for GRCh38 / hg38
cd ../../v1.6/GRCh38
wget -c https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/annotationsGRCh38_v1.6.tar.gz
wget -c https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz
wget -c https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz.tbi
wget -c https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.indel.tsv.gz
wget -c https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.indel.tsv.gz.tbi
tar -kzxvf annotationsGRCh38_v1.6.tar.gz


