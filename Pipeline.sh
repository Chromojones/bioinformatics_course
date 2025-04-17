#!/bin/bash #

######## Create a repository in the home directory ########

cd ~/
mkdir dnaseq
cd dnaseq/
mkdir data meta results logs scripts
mkdir data/trimmed data/untrimmed
touch README.md
echo "Created repository containing:"
ls -Rh
cd ~/dnaseq

######## Download and unzip the FASTQ files ########
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
mv NGS0001.R1.fastq.qz NGS0001.R1.fastq.gz
mv NGS0001.R2.fastq.qz NGS0001.R2.fastq.gz
mv hg19.fa.gz ~/dnaseq/data/
echo "Downloaded the following files:"
ls NGS* anno*
mv NGS* ~/dnaseq/data/untrimmed
mv anno* ~/dnaseq/data

######## Install packages for DNA-Seq using conda ########
cd ~/

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x ./Miniconda3-latest-Linux-x86_64.sh
. ./Miniconda3-latest-Linux-x86_64.sh

source ~/.bashrc
conda config --add channels defaults --add channels bioconda --add channels conda-forge
conda config --show channels
conda --version

conda install samtools
conda install bwa
conda install freebayes
conda install picard
conda install mosdepth
sudo apt install libvcflib-tools
conda install trimmomatic
conda install fastqc
conda install bcftools

######## Perform basic quality assessment ########

mkdir ~/dnaseq/results/fastqc_untrim
cd ~/dnaseq/data/untrimmed
fastqc -t 4 *.fastq
mv *fastqc* ~/dnaseq/results/fastqc_untrim/
echo "FastQC completed, reports deposited in Results directory for safekeeping"
ls -lRh ~/dnaseq/results/
cd ~/

######## Preview our FastQC Report ########

for zip in ~/dnaseq/results/fastqc_untrim/*.zip
do
unzip -l $zip >> ~/dnaseq/logs/fastqcsummary.txt
done
head ~/dnaseq/logs/fastqcsummary.txt
tail ~/dnaseq/logs/fastqcsummary.txt
sudo apt install w3m w3m-img 
w3m ~/dnaseq/results/fastqc_untrim/NGS0001.R1_fastqc.html | head -n 100

######## Trim reads and repeat quality assessment ########

trimmomatic PE  \
-threads 4 \
-phred33 \
~/dnaseq/data/untrimmed/NGS0001.R1.fastq.gz ~/dnaseq/data/untrimmed/NGS0001.R2.fastq.gz \
-baseout ~/dnaseq/data/trimmed/NGS0001_trimmed_R.fastq \
ILLUMINACLIP:~/anaconda3/pkgs/trimmomatic-0.39-*/share/trimmomatic-*/adapters/NexteraPE-PE.fa:2:30:10 \
TRAILING:3 MINLEN:50

ls ~/dnaseq/data/trimmed/
mkdir ~/dnaseq/results/fastqc_trim
cd ~/dnaseq/data/trimmed
fastqc -t 4 *.fastq
mv *fastqc* ~/dnaseq/results/fastqc_trim/
echo "FastQC completed, reports deposited in Results directory for safekeeping"
ls -lRh ~/dnaseq/results/fastqc_trim
cd ~/
for zip in ~/dnaseq/results/fastqc_trim/*.zip
do
unzip -l $zip >> ~/dnaseq/logs/fastqctrimsummary.txt
done
head ~/dnaseq/logs/fastqctrimsummary.txt
tail ~/dnaseq/logs/fastqctrimsummary.txt
w3m ~/dnaseq/results/fastqc_trim/NGS0001_trimmed_R_1P_fastqc.html | head -n 100

######## Alignment: Create an index and run BWAMEM  ########

mkdir ~/dnaseq/data/reference
mv ~/dnaseq/data/hg19.fa.gz  ~/dnaseq/data/reference
bwa index ~/dnaseq/data/reference/hg19.fa.gz
mkdir ~/dnaseq/data/alignment/
bwa mem -t 4 -v 1 -R '@RG\tID:11V6WR1.111.D1375ACXX.1.NGS0001\tSM:NGS0001\tPL:ILLUMINA\tLB:nextera-NGS0001-blood\tDT:2020-09-28\tPU:11V6WR1' -I 250,50  ~/dnaseq/data/reference/hg19.fa.gz ~/dnaseq/data/trimmed/NGS0001_trimmed_R_1P.fastq ~/dnaseq/data/trimmed/NGS0001_trimmed_R_2P.fastq > ~/dnaseq/data/alignment/NGS0001.sam
ls ~/dnaseq/data/alignment/
cd ~/dnaseq/data/alignment/

######## Convert SAM to BAM ########

samtools view -h -b NGS0001.sam > NGS0001.bam
rm NGS0001.sam
echo "SAM file removed."
samtools sort NGS0001.bam > NGS0001_sorted.bam
samtools index NGS0001_sorted.bam
ls

######## Mark Duplicates and Filtering ########

picard MarkDuplicates I=NGS0001_sorted.bam O=NGS0001_sorted_marked.bam M=marked_dup_metrics.txt 
samtools index NGS0001_sorted_marked.bam
samtools view -F 1796 -q 20 -o NGS0001_sorted_filtered.bam NGS0001_sorted_marked.bam
samtools index NGS0001_sorted_filtered.bam
ls

######## Standard Alignment Statistics ########

samtools flagstat NGS0001_sorted_marked.bam
samtools flagstat NGS0001_sorted_filtered.bam
samtools idxstats NGS0001_sorted_filtered.bam
picard CollectInsertSizeMetrics \I=NGS0001_sorted_filtered.bam O=insert_size_metrics.txt H=insert_size_histogram.pdf DEVIATIONS=10 M=0.05 
cd ~/dnaseq/data/
mv annotation.bed alignment/
cd ~/dnaseq/data/alignment/
mosdepth --by annotation.bed NGS0001 NGS0001_sorted_filtered.bam
cat NGS0001.mosdepth.summary.txt

######## Remove excess files ########

rm NGS0001_sorted_marked.*
rm NGS0001_sorted.*
rm NGS0001.bam
rm ~/dnaseq/data/trimmed/*.fastq
rm ~/dnaseq/data/untrimmed/*.fastq.gz
cd ~/

######## Variant Calling ########

zcat ~/dnaseq/data/reference/hg19.fa.gz > ~/dnaseq/data/reference/hg19.fa 
samtools faidx ~/dnaseq/data/reference/hg19.fa
freebayes --bam ~/dnaseq/data/alignment/NGS0001_sorted_filtered.bam --fasta-reference ~/dnaseq/data/reference/hg19.fa --targets ~/dnaseq/data/alignment/annotation.bed --vcf ~/dnaseq/results/NGS0001.vcf
bgzip ~/dnaseq/results/NGS0001.vcf
ls ~/dnaseq/results/
tabix -p vcf ~/dnaseq/results/NGS0001.vcf.gz
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1"  ~/dnaseq/results/NGS0001.vcf.gz > ~/dnaseq/results/NGS0001_filtered.vcf
ls -lrH ~/dnaseq/results/
bgzip ~/dnaseq/results/NGS0001_filtered.vcf
tabix -p vcf ~/dnaseq/results/NGS0001_filtered.vcf.gz
rm ~/dnaseq/results/NGS0001.vcf.gz
cd ~/

######## Unzip Annovar and Annotate Variants ########

mkdir annovar
cd annovar/
wget -O annovar.latest.tar.gz 'https://www.openbioinformatics.org/annovar/download/annovar.latest.tar.gz'
tar -zxvf annovar.latest.tar.gz
cd annovar/
chmod +x annotate_variation.pl
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp150 humandb/
./convert2annovar.pl -format vcf4 ~/dnaseq/results/NGS0001_filtered.vcf.gz > ~/dnaseq/results/NGS0001_filtered.avinput
./table_annovar.pl ~/dnaseq/results/NGS0001_filtered.avinput humandb/ -buildver hg19 -out ~/dnaseq/results/NGS0001_filtered -remove -protocol refGene,ensGene,clinvar_20180603,exac03,avsnp150 -operation g,g,f,f,f -otherinfo -nastring . -csvout
cd ~/
csvcut -c avsnp150 ~/dnaseq/results/NGS0001_filtered.hg19_multianno.csv | head
awk -F',' 'NR==1 || ($6 == "\"exonic\"" && $29 == ".")' ~/dnaseq/results/NGS0001_filtered.hg19_multianno.csv > ~/dnaseq/results/NGS0001_exonic_rare.csv
rm -r annovar/annovar/humandb/

######## Download, unzip, and annotate with SnpEFF ########

wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
cd snpEff/
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz.tbi
conda create -n bcftools_env -c bioconda bcftools
conda activate bcftools_env

bcftools annotate \
  -a 00-All.vcf.gz \
  -c CHROM,POS,ID,REF,ALT \
  -h <(echo "##INFO=<ID=DBSNP,Number=1,Type=String,Description=\"dbSNP rsID\">") \
  -o NGS0001_dbsnp.vcf \
  -O z \
  ~/dnaseq/results/NGS0001_filtered.vcf.gz

java -Xmx4g -jar ~/snpEff/snpEff.jar -c ~/snpEff/snpEff.config hg19 ~/snpEff/NGS0001_dbsnp.vcf > ~/dnaseq/results/SnpEffanno.vcf 
