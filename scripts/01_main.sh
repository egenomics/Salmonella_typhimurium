conda activate bactopia

cd /home/jlvillanueva/Documents/SALMONELLA/fastq
bactopia prepare --path /home/jlvillanueva/Documents/SALMONELLA/fastq/ --fastq-ext '_001.fastq.gz' > /home/jlvillanueva/Documents/SALMONELLA/scripts/samples.tsv
# Replace unkown species by Salmonella_enterica in samples.tsv

'''
sample	runtype	genome_size	species	r1	r2	extra
SAL1-14_12598	paired-end	0	Salmonella_enterica	/home/jlvillanueva/Documents/SALMONELLA/fastq/20591-SAL1-598_S15_L001_R1_001.fastq.gz	/home/jlvillanueva/Documents/SALMONELLA/fastq/20591-SAL1-598_S15_L001_R2_001.fastq.gz	
SAL2-13_13995	paired-end	0	Salmonella_enterica	/home/jlvillanueva/Documents/SALMONELLA/fastq/20591-SAL2-995_S16_L001_R1_001.fastq.gz	/home/jlvillanueva/Documents/SALMONELLA/fastq/20591-SAL2-995_S16_L001_R2_001.fastq.gz	
SAL3-14_11275	paired-end	0	Salmonella_enterica	/home/jlvillanueva/Documents/SALMONELLA/fastq/20591-SAL3-275_S21_L001_R1_001.fastq.gz	/home/jlvillanueva/Documents/SALMONELLA/fastq/20591-SAL3-275_S21_L001_R2_001.fastq.gz	
SAL4-14_28806	paired-end	0	Salmonella_enterica	/home/jlvillanueva/Documents/SALMONELLA/fastq/20591-SAL4-510_S22_L001_R1_001.fastq.gz	/home/jlvillanueva/Documents/SALMONELLA/fastq/20591-SAL4-510_S22_L001_R2_001.fastq.gz	
SAL5-14_28823	paired-end	0	Salmonella_enterica	/home/jlvillanueva/Documents/SALMONELLA/fastq/20591-SAL5-311_S17_L001_R1_001.fastq.gz	/home/jlvillanueva/Documents/SALMONELLA/fastq/20591-SAL5-311_S17_L001_R2_001.fastq.gz	
SAL6-16_18784	paired-end	0	Salmonella_enterica	/home/jlvillanueva/Documents/SALMONELLA/fastq/20591-SAL6-450_S23_L001_R1_001.fastq.gz	/home/jlvillanueva/Documents/SALMONELLA/fastq/20591-SAL6-450_S23_L001_R2_001.fastq.gz	
SAL7-16_18792	paired-end	0	Salmonella_enterica	/home/jlvillanueva/Documents/SALMONELLA/fastq/20591-SAL7-614_S18_L001_R1_001.fastq.gz	/home/jlvillanueva/Documents/SALMONELLA/fastq/20591-SAL7-614_S18_L001_R2_001.fastq.gz	
REF-GSJ_2017-Sal-008	paired-end	0	Salmonella_enterica	/home/jlvillanueva/Documents/SALMONELLA/fastq/GSJ_2017-Sal-008_R1_001.fastq.gz	/home/jlvillanueva/Documents/SALMONELLA/fastq/GSJ_2017-Sal-008_R2_001.fastq.gz
'''

cd /home/jlvillanueva/Documents/SALMONELLA/scripts
bactopia -profile singularity \
    --samples samples.tsv \
    --max_cpus 4 \
    -qs 2 \
    --outdir /playground/dataset_salmonella/

cd /playground/dataset_salmonella

# search for strains of salmnonella with similar characteristics (min read lenght 249), coverage of 50 to our dataset
bactopia search --query 90371 -gsize 5027649 -mc 50 -mrl 249

# Run Rscript to select relevant strains for phylogenetic context
### data_explore.Rmd

# run for all the accesion samples located using bactopia
cd /playground/dataset_salmonella
bactopia -profile singularity \
    --accessions /playground/dataset_salmonella/search/selected_accessions.txt \
    --max_cpus 4 \
    -qs 8 \
    --outdir /playground/dataset_salmonella/

workflows='genotyphi mashtree seqsero2 sistr busco checkm plasmidfinder quast rgi'

# busco: Assembly completeness based on evolutionarily informed expectations
# checkm: Assess the assembly quality of your samples
# plasmidfinder: Plasmid identification from assemblies
# quast: A module for assessing the quality of assembled contigs
# rgi: Predict antibiotic resistance from assemblies
# mashtree: Quickly create a tree using Mash distances
# snippy: Rapid variant calling from Illumina sequence reads with optional core-SNP phylogeny
  #gubbins: phylogenetic tree
  #iqtree: phylogenetic tree
# genotyphi: Salmonella Typhi genotyping with Mykrobe outputs
# seqsero2: Salmonella serotype prediction from reads or assemblies
# sistr: Serovar prediction of Salmonella assemblies

cd /home/jlvillanueva/Documents/SALMONELLA/scripts
cat samples.tsv | cut -f1 | tail -8 > includes.txt

for wf in $workflows; do
bactopia -profile singularity --wf $wf \
  --bactopia /playground/dataset_salmonella/ \
  --include includes.txt --max_cpus 4 \
  -qs 8 -resume
done
s
# Reference from paper
# Merge reference strain with 4 plasmids
https://github.com/kblin/merge-gbk-records
https://www.nature.com/articles/s41564-020-00836-1

merge-gbk-records D23580_liv_o.gb D23580_liv_pBT1.gb D23580_liv_pBT2.gb D23580_liv_pBT3.gb D23580_liv_pSLT-BT.gb > D23580_4plasmids.gb

bactopia -profile singularity --wf snippy \
  --bactopia /playground/dataset_salmonella/ \
  --include includes.txt --max_cpus 8 \
  -qs 4 --reference /home/jlvillanueva/Documents/SALMONELLA/reference/D23580_4plasmids.gb -resume

# Other reference https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000006945.2/
bactopia -profile singularity --wf snippy \
  --bactopia /playground/dataset_salmonella/ \
  --include includes.txt --max_cpus 8 \
  -qs 4 --reference /home/jlvillanueva/Documents/SALMONELLA/reference/GCA_000006945.2/ncbi_dataset/data/GCA_000006945.2/genomic.gbff


bactopia -profile singularity --wf kraken2 \
  --bactopia /playground/dataset_salmonella/ \
  --kraken2_db /home/jlvillanueva/bactopia/kraken_db/k2_pluspf_20221209.tar.gz \
  --include includes.txt --max_cpus 8 \
  -qs 4 \
  --max_memory 20

bactopia-summary -profile singularity --bactopia /home/jlvillanueva/Documents/SALMONELLA/bactopia_v3_results

# reference strain (GSJ/2017-Sal-008)
bactopia -profile singularity --accession SRX6075121 --max_cpus 4 --outdir sra_test

fasterq-dump SRR9307305
mv SRR9307305_1.fastq.gz GSJ_2017-Sal-008_R1_001.fastq.gz
mv SRR9307305_2.fastq.gz GSJ_2017-Sal-008_R2_001.fastq.gz

cd /playground/dataset_salmonella/bactopia-runs/snippy-20240321-155218
trimal -in core-snp-clean.full.aln -out trimal_out.aln -automated1

raxml-ng --all --msa trimal_out_seq90_res07.aln --model GTR+G --prefix T5 --threads 24 --seed 2

raxml-ng --all --msa trimal_out_seq90_res07.aln --model GTR+G --prefix ST --tree pars{25},rand{25} --bs-trees 200 --threads 32 --seed 2

#testing models
iqtree2 -s trimal_out_seq90_res07.aln -alrt 1000 -B 1000 -T 32

#GTR+F+R4
iqtree2 -s trimal_out_seq90_res07.aln -alrt 1000 -B 1000 -T 2 -m GTR+F+R4

# TIM+F+R4
iqtree2 -s trimal_out_seq90_res07.aln -alrt 1000 -B 1000 -T 2 -m TIM+F+R4

# K3Pu+F+R4
iqtree2 -s trimal_out_seq90_res07.aln -alrt 1000 -B 1000 -T 32 -m K3Pu+F+R4
