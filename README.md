# VariantCallingSimulations
Scripts and code used to perform variant calling simulations for non-model organisms.

Most scripts (not awk though) that I've written have a `-h` option, so take a look, and let me know if anything is unclear!

## Categories
1. [Genome simulation](README.md#simulation-related-scripts)
1. [Variant calling evaluation](README.md#evaluation-related-scripts)

### Compilation of C++ scripts:

C++ programs here will be compiled by a simple call to `make`, although they won't be installed to a location in your PATH.

## Simulation-related scripts:

### `simulateDivergedHaplotype.pl`

Given a haploid reference genome in unwrapped FASTA format, it simulates SNPs and indels with uniform spatial distributions along scaffolds, and in quantities determined by the percentage passed as the only positional argument. The indel rate is fixed at 0.04*SNP rate, and the length distribution of indels follows a modified geometric distribution with adjustable parameter (passed in by the `-g` option). The geometric distribution is modified such that the 0 class has density 0 (since a 0-length indel means nothing). The underlying implementation is simply that any time a 0 is drawn, the geometric is resampled until a non-zero length is generated. Insertions and deletions should occur at equal frequency (and the choice between the two is determined by rounding a random number between 0 and 1). The user can choose to omit or include indels (include indels by passing the `-n` flag).

The script outputs the simulated FASTA (unwrapped) to a file specified by the -o option. It also outputs a "SNP log" and an "indel log".

The SNP log is a TSV file with 4 columns:

1. Scaffold
1. Position in input FASTA coordinate space on the scaffold
1. Original allele
1. New allele

The indel log is a TSV file with 5 columns:

1. Scaffold
1. Position in input FASTA coordinate space on the scaffold
1. "ins" for insertion, or "del" for deletion
1. Indel length
1. String of bases that were inserted or deleted

Note that for deletions, column 2 denotes the first base that was deleted, not the base before the deletion. On the other hand, for insertions, column 2 denotes the base prior to the inserted bases.

Example call simulating 0.5% divergence, geometric parameter = 0.35, and including indels:

`simulateDivergedHaplotype.pl -i my_reference_unwrapped.fasta -o my_diverged_genome.fasta -n -g 0.35 0.5`

The SNP log for this example would be called `my_diverged_genome_SNPs.log` and the indel log would be called `my_diverged_genome_indels.log`.

### `mergeSNPlogs`

Mutations that occurred along the reference-ancestor branch need to be combined with mutations that occurred along the ancestor-haploid branch. On top of that, when indels are simulated along the reference-ancestor branch, this changes the coordinate space of the ancestor's FASTA, so we cannot simply perform a set join of the two SNP logs, we need to adjust the positions of ancestor-haploid branch SNPs back into the coordinate space of the reference. In order to perform this position adjustment, we need to know the positions and sizes of indels along the reference-ancestor branch, which we pass in via the ref-anc branch indel log (the `-i` option). Of course, we then need the ref-anc branch and anc-haploid branch SNP logs, which are passed in via the `-b` and `-c` options, respectively. The merged and adjusted SNP log is output to `STDOUT`, which we redirect to a file in the example call.

Example call:

`mergeSNPlogs -i my_ref_anc_indels.log -b my_ref_anc_SNPs.log -c my_anc_hap1_SNPs.log > haploid1_merged_SNPs.log`

### `diploidizeSNPlog`

Here we create the SNP log appropriate for a diploid formed by the two simulated haploids, in the coordinate space of the reference. In order to guarantee the correct scaffold ordering of the output, we pass in a FASTA index (.fai file, generated by `samtools faidx [reference FASTA]`) using the `-i` option. Then we pass in the two merged SNP logs using the `-a` and `-b` options (the order doesn't matter). Iterating over the scaffolds in the order presented in the .fai file, this program identifies SNPs present only in haploid A, SNPs present only in haploid B, and SNPs present in both haploids, and creates diploid records out of them. In particular, SNPs present only in one of the two haploids will be heterozygous positions in the diploid, consisting of the reference allele and the new allele, so they are output as the degenerate IUPAC code for that type of heterozygote. SNPs present in both haploids are output as the degeneration of the two alleles, so if the two alleles are the same, the site is output as homozygous for that allele, and if the alleles are different, the site is output as heterozygous for the two new alleles.

Note that this implies there may be sites with true genotype `1/2`, i.e. a triallelic site.

The output is a SNP log formatted as before, but containing IUPAC-encoded degenerate bases in the 4th column wherever heterozygous sites occur.

Example call:

`diploidizeSNPlog -i my_reference_unwrapped.fasta.fai -a haploid1_merged_SNPs.log -b haploid2_merged_SNPs.log > diploid_SNPs.log`

## Evaluation-related scripts

### `VCFtoUnfilteredINSNP_skipInsertions.awk`

Generates an unfiltered INSNP from the ERCGVCF.vcf produced directly by GATK HaplotypeCaller.

Example call:

`VCFtoUnfilteredINSNP_skipInsertions.awk Dyak_2Mreads/Dyak_2Mreads_realigned_HC_ERCGVCF.vcf > Dyak_2Mreads/Dyak_2Mreads_MD_IR_ERCGVCF_unfiltered_INSNP.tsv`

### `VCFtoUnfilteredINSNP_skipInsertions_GGVCFs.awk`

Generates an unfiltered INSNP from the GGVCFs.vcf produced by GATK GenotypeGVCFs.

Example call:

`VCFtoUnfilteredINSNP_skipInsertions_GGVCFs.awk Dyak_2Mreads/Dyak_2Mreads_realigned_HC_GGVCFs.vcf > Dyak_2Mreads/Dyak_2Mreads_MD_IR_GGVCFs_unfiltered_INSNP.tsv`

### `VCFtoUnfilteredINSNP_skipInsertions_samtools.awk`

Generates an unfiltered INSNP from the bgzipped VCF produced by BCFtools call.

Example call:

`zcat Dyak_2Mreads/Dyak_2Mreads_realigned_mpileupcall.vcf.gz | VCFtoUnfilteredINSNP_skipInsertions_samtools.awk > Dyak_2Mreads/Dyak_2Mreads_MD_IR_mpileup_unfiltered_INSNP.tsv`

### `compareSNPlogs`

Example call:

`compareSNPlogs -i [.fai FASTA index of the reference FASTA] -e [expected diploid SNP log] -o [observed INSNP file] -n [output FN log] -p [output FP log] -t [output TP log] -r [output ER log] -m [min depth]`

Here, the .fai file is used for calculating the number of true negatives, since we need the scaffold length to compute this from TP, FP, and FN, which we know from the logs.  Of course, you pass in the expected diploid SNP log generated by diploidizeSNPlog as above using the `-e` option, and you pass in the observed INSNP generated by applying one of the above awk scripts to the VCF output by the pseudoreference pipeline using the `-o` option.  The optional `-n`, `-p`, `-r`, and `-t` parameters allow you to get an INSNP-formatted log of all of the false negatives, false positives, true positives, and ambiguous miscalls (ERs), respectively, for further analysis of trends in these four categories.  For example, a significant fraction of the false negatives may be due to low sequencing depth, or low depth used by the variant caller, so you can process this INSNP with awk to generate a pattern file for fgrep, and then determine raw depth from an fgrep of the output of samtools depth, or determine the variant caller-used depth by using fgrep on the all-sites VCF.

The final parameter, `-m` or `--min_depth`, specifies a minimum depth, provided as a fifth column of the expected diploid SNP log, for a site to be considered "callable".  This allows us to calculate the statistics on only putatively "callable" sites -- sites with sufficient read coverage to have alleles detected.

Note that the INSNP-like log TSVs output can easily be converted to BED format, and the intervals of true negatives inferred from the set complement of the union of FP, FN, TP, and ER intervals. For example, to convert such an INSNP-like TSV to a BED:

`awk 'BEGIN{FS="\t";OFS="\t";}{print $1, $2-1, $2;}' Dyak_2Mreads/Dyak_2Mreads_MD_IR_mpileup_FPs.tsv | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i - > Dyak_2Mreads/Dyak_2Mreads_MD_IR_mpileup_FPs_merged.bed`

And to find the set of TN intervals:

`cat Dyak_2Mreads/Dyak_2Mreads_MD_IR_mpileup_{ER,FN,FP,TP}s_merged.bed | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i - | bedtools complement -i - -g <(cut -f1,2 Dyak_NY73_Quiver_Scaffolded_w60.fasta.fai) > Dyak_2Mreads/Dyak_2Mreads_MD_IR_mpileup_TNs_merged.bed`

### `subsetVCFstats.pl`

Example call:

`subsetVCFstats.pl -i Dyak_2Mreads/Dyak_2Mreads_realigned_HC_GGVCFs.vcf -b Dyak_2Mreads/Dyak_2Mreads_MD_IR_GATK_GGVCFs_FPs_merged.bed | less`

This script is fairly generic in that it will take any tab-separated file whose first two columns are Scaffold and Position, and will print out the lines of that file corresponding to the intervals provided in the BED file. The script was originally written to subset lines out of stats files made by `bcftools query` or GATK VariantsToTable, but can be applied to the VCF itself, or really any file fitting these criteria.