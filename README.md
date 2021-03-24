# Cross-Disease Fine-Mapping QC and Analysis Pipeline


## Contact
[Matthew R Lincoln](mailto:matthew.lincoln@yale.edu)
[Chris Cotsapas](mailto:cotsapas@broadinstitute.org)


## Introduction
This repository contains code used for quality control, association, joint likelihood mapping and fine-mapping in [Joint analysis reveals shared autoimmune disease associations and identifies common mechanisms](https://localhost).


## Disclaimer
Code is configured for use on the Yale HPC [farnam](https://docs.ycrc.yale.edu/clusters-at-yale/clusters/farnam/) cluster using the slurm workload manager. The code here is provided as-is. We have validated that it runs without errors on our system, but we cannot guarantee performance in other computing environments.


## Pipeline Overview
The analysis pipeline is described in detail in the supplementary materials. There are ten principal steps, as follows:
1. Setup.
2. Manifest tidying.
3. Liftover to hg19.
4. Quality control.
5. Recode individuals.
6. Association test.
7. Phasing and imputation.
8. Compile QC data.
9. Identify shared associations across diseases.
10. Cross-disease meta-analysis and fine-mapping.


## Software Dependencies
- [PLINK 1.90](https://www.cog-genomics.org/plink/1.9/)
- [PLINK 2.00](https://www.cog-genomics.org/plink/2.0/)
- [R](https://www.r-project.org) 3.5.0 or later
- [Python](https://www.python.org/download/releases/2.7/) 2.7
- [liftOver](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver)
- [flashpca](https://github.com/gabraham/flashpca)
- [SHAPEIT2](http://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html)
- [IMPUTE2](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html)
- [SNPTEST](https://mathgen.stats.ox.ac.uk/genetics_software/snptest/snptest.html)
- [VCFtools](https://vcftools.github.iohttps://vcftools.github.io)
- [samtools](https://github.com/samtools) tabix
- [Joint Likelihood Mapping (JLIM)](https://github.com/cotsapaslab/jlim)
- [FINEMAP](http://www.christianbenner.com)


## Data Dependencies
The pipeline requires raw ImmunoChip genotype data for each of the six diseases under study. Data sources and/or relevant contact details are provided in the supplementary information.

In addition, the following reference data is required:
- [1,000 Genomes (Phase 3) genotypes](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/)
- [1,000 Genomes (Phase 3) phased haplotypes](https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz)


## Running the Pipeline
The pipeline is invoked from src/immchip.master.sh, which calls relevant dependencies to execute the main analysis steps in order. Detailed descriptions are provided at the top of each script, and also as inline comments. The pipeline is not intended to be run as standalone, unsupervised code; user input is required at several stages.


## Files Included in This Distribution
| File | Description |
| ---- | --- |
| README.md | This file |
| src/  | Contains all scripts required in by the pipeline |
| src/immchip.master.sh | The main analysis script |
| src/resolve.manifests.1.R | Manifest resolution |
| src/resolve.manifests.2.R | Manifest resolution |
| src/snp.similarity.sh | Manifest resolution |
| src/dup.similarity.sh | Manifest resolution |
| src/liftover.py | Wrapper script for liftOver |
| src/immchip.ced.qc.sh | Celiac disease quality control|
| src/immchip.ibd.qc.sh | Inflammatory bowel disease quality control |
| src/immchip.ms.qc.sh | Multiple sclerosis quality control |
| src/immchip.sle_g.qc.sh | Systemic lupus erythematosus Genentecch quality control |
| src/immchip.sle_o.qc.sh | Systemic lupus erythematosus OMRF quality control |
| src/immchip.t1d.qc.sh | Type 1 diabetes GRID quality control |
| src/immchip.t1d_asp.qc.sh | Type 1 diabetes ASP quality control |
| src/immchip.ra.qc.sh | Rheumatoid arthritis quality control |
| src/resolve.ra.manifest.R | Manifest resolution |
| src/reassign.sex.R | Identify sex inconsistencies |
| src/plot.flashpca.R | Identify population outliers |
| src/remove.pop.outliers.R | Identify population outliers |
| src/cluster.pca.R | Identify population outliers |
| src/het.miss.R | Identify heterozygosity outliers |
| src/identify.dups.and.rels.to.remove.R | Identify duplicate and related subjects |
| src/rename.subjects.R | Harmonize sample nomenclature |
| src/immchip.assoc.sh | Initial association testing |
| src/immchip.impute.sh | Genotype imputation |
| src/identify.shapeit.flips.awk | Identify strand flips prior to imputation |
| src/identify.shapeit.monos.awk | Identify monomorphic SNPs prior to imputation |
| src/immchip.postimputation.qc.sh | Post-imputation quality control |
| src/account.snp.sample.fate.R | Account for SNPs and samples through QC |
| src/immchip.jlim.impute.sh | Identify shared traits |
| src/jlim.cond.select.imputed.snps.R | Identify conditioning SNPs for JLIM association analysis |
| src/jlim.cond.metafor.R | Meta-analysis of conditional associations |
| src/jlim.impute.pairs.R | Identify trait pairs to analyze with JLIM |
| src/jlim.indep.metafor.R | Meta-analysis of conditionally independent associations |
| src/jlim.permute.sample.pheno.R | Permute sample phenotypes for JLIM |
| src/jlim.impute.make.perm.matrix.R | Create permutation matrix for JLIM |
| src/jlim.pairs.to.repeat.R | Identify trait pairs to repeat JLIM analysis with differing intervals |
| src/compile.jlim.results.R | Compile JLIM results |
| src/immchip.finemap.impute.sh | Fine-map across shared trait clusters |
| src/jlim.impute.finemap.metafor.R | Meta-analysis for shared trait clusters |
