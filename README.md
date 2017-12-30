
# Preparation Steps

## Download genome from 23andMe to ./
unzip genome*
chmod 400 genome*
mkdir Raw
chmod 400 Raw/
sudo mv genome_YOURS.txt Raw/

## Download MacArthur ClinVar
mkdir clinvar
git clone https://github.com/macarthur-lab/clinvar.git
gzip -d output/b37/single/clinvar_allelese.single.b37.tsv.gz
ln -s output/b37/single/clinvar_allelese.single.b37.tsv clinvar/clinvar_allelse.tsv

# Analysis

## Primary Cross-referencing
perl crossReferenceClinVar.pl Raw/genome_YOURS.txt clinvar/clinvar_alleles.tsv m YOUR.result.tsv 1>YOUR.result.het.tsv 2>YOUR.result.hom.tsv;

## Move Results
sudo mv YOUR.result.*.tsv Output/;
chmod 400 Output/*

## Cleaner View
cut -f 1-11,16,18,20-23,27,29-32 Output/YOUR.result.h* | sort -u | less -S

## Summary Zygosity
cut -f 1,2 Output/YOURS.result.tsv | sort | uniq -c | awk '{print $1"\t"$2"\t"$3;}' | sort -k1n

## Genotyping

perl germlineFault.pl Raw/genome_YOURS_v4_Full_20170204085740.txt Raw/genome_YOUR_MOM_v4_Full_20170211153523.txt Raw/genome_YOUR_DAD_v4_Full_20170211150722.txt m YOURS.faults.tsv

grepbc 1:1 YOURS.faults.tsv 4 2 > YOURS.faults.homMother.tsv
grepbc -:- YOURS.faults.homMother.tsv 5 2 > YOURS.faults.homMother.noFather.tsv

grepbc 1:1 YOURS.faults.tsv 5 2 > YOURS.faults.homFather.tsv
grepbc -:- YOURS.faults.homFather.tsv 4 2 > YOURS.faults.homFather.noMother.tsv

(NOTE: grepbc can be found in my personal_helpful repo: https://github.com/AdamDS/personal_helpful)
