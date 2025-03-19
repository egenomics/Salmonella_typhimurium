## Monophasic Salmonella Typhimurium circulating in seagulls from Barcelona, Spain is closely related to strains involved in human infections. 

The exploratory analysis to select and download similar NGS data from other sources is located in the folder:

**search_strains**

The code for the main analysis using bactopia and the generation of the tree figure is located in the folder:

**scripts**


Modules from bactopia used with their version and basic description:
- busco v5.7.0: Assembly completeness based on evolutionarily informed expectations
- checkm v1.2.2: Assess the assembly quality of samples
- plasmidfinder v2.1.6: Plasmid identification from assemblies
- quast v5.2.0: A module for assessing the quality of assembled contigs
- rgi v6.0.3: Predict antibiotic resistance from assemblies
- mashtree v1.4.6: Quickly create a tree using Mash distances
- snippy v4.6.0: Rapid variant calling from Illumina sequence reads with optional core-SNP phylogeny
- gubbins v3.3.4 phylogenetic tree
- iqtree v2.2.2.7: phylogenetic tree
- genotyphi v2.0.0: Salmonella Typhi genotyping with Mykrobe outputs
- seqsero2 v1.2.1: Salmonella serotype prediction from reads or assemblies
- sistr v1.1.2: Serovar prediction of Salmonella assemblies
