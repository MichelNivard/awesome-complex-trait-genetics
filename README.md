# awesome-complex-trait-genetics

[![Awesome](https://awesome.re/badge.svg)](https://awesome.re)

A list of awesome tools for human complex trait genetics, biased/slanted towards post-GWAS analysis.

How can you add your tool? via a pull request, if you dont know what that is, read the contributing guidlines: (https://github.com/MichelNivard/awesome-genetics/blob/main/contributing.md)

**This list may for the basis for a static review**, in which case serious contributors, so those who contribute to the list and to the preamle paragraphs for each section, will be ivited as collaborators.


## genetic architecture

**Sumstats based:**

- [LDSC](https://github.com/bulik/ldsc) ldsc is a python command line tool for estimating heritability and genetic correlation from GWAS summary statistics. ldsc also computes LD Scores. A [Python3 port](https://github.com/belowlab/ldsc) is available.

- [LDSR](https://github.com/Ararder/ldsR) LDSC implemented in R. Several quality of life features that simplifies the user interface to LDSC, and makes it dirt simple to run heritability estimates, genetic correlations and partitioned heritability from in-memory summary statistics.
  
- [GCTB](https://cnsgenomics.com/software/gctb/#Overview) GCTB is a software tool that comprises a family of Bayesian linear mixed models for complex trait analyses using genome-wide SNPs. It was developed to simultaneously estimate the joint effects of all SNPs and the genetic architecture parameters for a complex trait. There are now extensions to estimate the same bayesian linear model parameters based on summary data. 

**Raw data based:**

- [GCTA](https://yanglab.westlake.edu.cn/software/gcta/) GCTA (Genome-wide Complex Trait Analysis) is a software package initially developed to estimate the proportion of phenotypic variance explained by all genome-wide SNPs for a complex trait.

- [RHE-mc/GENIE](https://github.com/sriramlab/GENIE/tree/main) RHE-mc is a method to estimate the proportion of phenotypic variance explained by SNPs, and GENIE extends this to model GxE effects

- [BOLT-LMM/BOLT-REML](https://alkesgroup.broadinstitute.org/BOLT-LMM/BOLT-LMM_manual.html) The BOLT-LMM software package currently consists of two main algorithms, the BOLT-LMM algorithm for mixed model association testing, and the BOLT-REML algorithm for variance components analysis (i.e., partitioning of SNP-heritability and estimation of genetic correlations).

- [LDAK/SumHer/PCGC](https://dougspeed.com/) LDAK is a software package for analysing association study data. PCGC (phenotype-correlation genotype-correlation) Regression is an alternative to REML when estimating heritability for binary traits (i.e., diseases).

### Univariate models (heritability/poligenicity/stratified/geneset enrichment etc)

- [i-LDSC](https://github.com/lcrawlab/i-LDSC)   interaction-LD score (i-LDSC) regression: Model an additional score that measures the amount of non-additive genetic variation that is tagged by each variant in the data.
  
- [ACLR](https://github.com/arminschoech/ACLR) Autocorrelation LD regression: a tool to efficiently estimate the autocorrelation of latent effects in large genetic data sets. (WARNING Python 2.7 code)
  
- [HAMSTA](https://github.com/tszfungc/hamsta) HAMSTA is a python package that estimate heritability explained by local ancestry using summary statistics from admixture mapping studies.
  
- [MAGMA](https://cncr.nl/research/magma/) MAGMA: Generalized gene-set analysis of GWAS data.
  
- [MiXeR](https://github.com/precimed/mixer) Causal mixture model (MiXeR) is a tool designed to estimate the polygenic architecture of a single trait, including the total number of causal variants and the distribution of their effect sizes.



## Genetic correlation (LD score derivatives/extensions)


- [HDL](https://github.com/zhenin/HDL) High-Definition Likelihood (HDL) is a likelihood-based method for estimating genetic correlation using GWAS summary statistics. Compared to LD Score regression (LDSC), It reduces the variance of a genetic correlation estimate by about 60%.

### Stratified/local genetic correlatons

- [LAVA](https://github.com/josefin-werme/LAVA) `LAVA` (Local Analysis of [co]Variant Association) is a tool developed for local genetic correlation (rg) analysis.

### Ancestry aware Genetic correlations:

- [s-ldxr](https://github.com/huwenboshi/s-ldxr) `S-LDXR` is a method to stratify squared trans-ethnic genetic correlation by genomic annotations from GWAS summary statistics.

- [Popcorn](https://github.com/brielin/Popcorn) `Popcorn` is a program for estimaing the correlation of causal variant effect sizes across populations in GWAS.

- [mama](https://github.com/JonJala/mama) `mama` is a Python-based command line tool that meta-analyzes GWAS summary statistics generated from distinct ancestry groups.

## Model trait relationships beyond correlation

### Genetic SEM/Factor models

**sumstat based:**

- [GenomicSEM](https://github.com/GenomicSEM/GenomicSEM) R package for Genomic Structural Equation Models. user defined models of the relation between complex traits based on GWAS summary data.
  
- [GUIDE](https://github.com/daniel-lazarev/GUIDE) Genetic Unmixing by Independent Decomposition (`GUIDE`), uses ICA to estimate statistically independent latent factors that best express the patterns of association across many traits.
  
- [FactorGO](https://github.com/mancusolab/FactorGo) FactorGo is a scalable variational factor analysis model that learns pleiotropic factors using GWAS summary statistics!
  
- [GNA](https://github.com/GenomicNetworkAnalysis/GNA) GNA is an R package for performing network analysis of genetic overlap derived from GWAS summary statistics
  
- [partialLDSC](https://github.com/GEMINI-multimorbidity/partialLDSC) is an R-package to estimate partial genetic correlations from GWAS summary statistics, and compare them to their unadjusted counterparts, to quantify the contribution of a given confounder in explaining genetic similarity between conditions.

**these tools also give you the % of SNPs that is pleiotropic:**

- [bivarite MiXeR](https://github.com/precimed/mixer) Bivariate causal mixture model (MiXeR) is a tool designed to estimate the degree of genetic overlap and shared causal variants between two traits.
  
- [trivariate MiXeR](https://github.com/precimed/mix3r) Trivariate causal mixture model (MiXeR) is a tool designed to estimate the degree of genetic overlap and shared causal variants between three traits.


**Raw data based:**

- [Gsens](https://github.com/LeonardFrach/Gsens): Genetically informed sensitivity analysis to estimate role of genetic confounding in phenotypic associations and for causal inference from observation data.


### Two sample Mendelian Randomisation

- [TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR) An R package for performing Mendelian randomization using GWAS summary data.
  
- [SMR](https://yanglab.westlake.edu.cn/software/smr/#Overview)The SMR software tool was originally developed to implement the SMR & HEIDI methods to test for pleiotropic association between the expression level of a gene and a complex trait of interest using summary-level data from GWAS and expression quantitative trait loci (eQTL) studies

### MR/Genetic architecture hybrid models

- [lhcMR](https://github.com/LizaDarrous/lhcMR) lhcMR is an R package that performs bi-directional causal estimation between a pair of traits, while accounting for the presence of a potential heritable confounder acting on the pair.
  
- [CAUSE](https://github.com/jean997/cause) an MR method, Causal Analysis Using Summary Effect Estimates (CAUSE), that accounts for correlated and uncorrelated horizontal pleiotropic effects.
  
- [MR-cML](https://github.com/xue-hr/MRcML) a simple and effective approach to MR, based on constrained maximum likelihood and model averaging, called cML-MA, applicable to GWAS summary data.
  
- [LCV](https://github.com/lukejoconnor/LCV) LCV is an LD score based method for inferring genetically causal relationships using GWAS data.

### Mendelian randomization in _cis_

- [MR-link-2](https://github.com/adriaan-vd-graaf/mrlink2) _cis_ Mendelian randomization that is more robust to violations of the pleiotropy assumption from a single region

## Colocalisation/finemapping of causal variants

- [coloc](https://github.com/chr1swallace/coloc) R package to perform genetic colocalisation analysis, test whether two related phenotypes share common genetic causal variant(s) in a given region.
  
- [fastENLOC](https://github.com/xqwen/fastenloc) This repository contains the software implementation of FastENLOC, which enables integrative genetic association analysis of molecular QTL data and GWAS data.

- [SharePro](https://github.com/zhwm/SharePro_coloc) This repository contains the software implementation of SharePro, which uses an effect group-level approach to integrate LD modeling and colocalization assessment to account for multiple causal variants in colocalization analysis.

- [FINEMAP](http://www.christianbenner.com) FINEMAP: Efficient variable selection using summary data from genome-wide association studies.
  
- [polyfun](https://github.com/omerwe/polyfun) Package contains **PolyFun** for functionally-informed fine-mapping, **PolyLoc** for polygenic localization of complex trait heritability.
  
- [OPERA](https://github.com/wuyangf7/OPERA): (omics pleiotropic association) method tests combinatorial pleiotropic associations between multiple molecular phenotypes (eQTL, DNAm) with a complex trait using summary-level data.
  
- [SuSiEx](https://github.com/getian107/SuSiEx): Cross-population finemapping using summary statistics and LD reference panels.

## gene-level analysis (TWAS)

- [FUSION](https://github.com/gusevlab/fusion_twas) FUSION is a suite of tools for performing transcriptome-wide and regulome-wide association studies (TWAS and RWAS).
  
- [FOCUS](https://github.com/mancusolab/ma-focus) FOCUS (Fine-mapping Of CaUsal gene Sets) is software to fine-map transcriptome-wide association study statistics at genomic risk regions

- [cTWAS](https://github.com/xinhe-lab/ctwas) Causal-TWAS (cTWAS) jointly models genetically regulated expression and the direct effects of genetic variants on the phenotype


## Simulation

- [GWASBrewer](https://github.com/jean997/GWASBrewer) A flexible tool for simulating realistic GWAS summary statistics for one, or many, traits.

- [magenpy](https://github.com/shz9/magenpy): A `python`-based library that provides utilities for simulating complex traits with various genetic architectures.

## Genomic data wrangling

- [HAIL](https://github.com/hail-is/hail) Hail is an open-source, general-purpose, Python-based data analysis tool with additional data types and methods for working with genomic data.
  
- [bigsnpr](https://github.com/privefl/bigsnpr) R package for the analysis of massive SNP arrays, primarily designed for human genetics.
  
- [ukbrapR](https://github.com/lcpilling/ukbrapR) ukbrapR (phonetically: 'U-K-B-wrapper') is an R package for working in the UK Biobank Research Analysis Platform (RAP). The aim is to make it quicker, easier, and more reproducible.
  
- [MungeSumstats](https://github.com/Al-Murphy/MungeSumstats) R package designed to facilitate the standardisation of GWAS summary statistics.

- [tidyGWAS](https://github.com/Ararder/tidyGWAS) R package that harmonises GWAS data, repairs/imputes missing columns, removing a minimal amount of rows.  
  
- [gwasRtools](https://github.com/lcpilling/gwasRtools) R package to (1) identify loci and independent lead SNPs (using online or local reference panel) and (2) annotate variants with nearest gene from GENCODE database.
  
- [qgg](https://psoerensen.github.io/qgg/) & [gact](https://psoerensen.github.io/gact/) **qgg** provides an infrastructure for efficient processing of large-scale genetic and phenotypic data while **gact** is designed for establishing and populating a comprehensive database focused on genomic associations with complex traits, provies R implementations of popular follow up analysis (LDscore regresison, MAGMA, VEGAS, PoPS, etc). 

- [bcftools](https://samtools.github.io/bcftools/bcftools.html) bcftools is a tool for creating, editing, and manipulating VCF and BCF files

- [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) An R packages for working with genomic intervals

- [bedtools](https://bedtools.readthedocs.io/en/latest/) Collectively, the bedtools utilities are a swiss-army knife of tools for a wide-range of genomics analysis tasks. A very fast and easy way to intersect, merge, count, complement, and shuffle genomic intervals from multiple files in widely-used genomic file formats such as BAM, BED, GFF/GTF, VCF

- [giggle](https://github.com/ryanlayer/giggle) GIGGLE indexes several BED/VCFs to perform efficient genomic interval searches across all files in the index

- [magenpy](https://github.com/shz9/magenpy): A `python`-based library that provides utilities for interfacing with genotype data (BED format), computing large-scale LD matrices, harmonizing various GWAS data sources (genotypes, LD, sumstats, annotations), and simulating complex traits.

## Polygenic scores

- [PRSice](https://choishingwan.github.io/PRSice/) PRSice (pronounced 'precise') is a Polygenic Risk Score software for calculating, applying, evaluating and plotting the results of polygenic risk scores (PRS) analyses.
  
- [LDpred2](https://privefl.github.io/bigsnpr/articles/LDpred2.html) LDpred-2 is one of the dedicated PRS programs which is an R package that uses a Bayesian approach to polygenic risk scoring.
  
- [GCTB](https://cnsgenomics.com/software/gctb/#Overview) GCTB is a software tool that comprises a family of Bayesian linear mixed models for complex trait analyses using genome-wide SNPs. It was developed to simultaneously estimate the joint effects of all SNPs and the genetic architecture parameters for a complex trait.

- [LDpred-funct](https://github.com/carlaml/LDpred-funct) LDpred-funct is an extension of LDpred that uses functional annotations of the SNPs to modify the prior.

- [DPR](https://github.com/biostatpzeng/DPR) DPR is a Bayesian PRS method that uses a non-parametric dirichlet process prior
  
- [VIPRS](https://github.com/shz9/viprs) VIPRS is a fast Bayesian PRS method that uses Variational Inference techniques to efficiently approximate the posterior for the effect sizes. The python package provides utilities for performing inference as well as computing polygenic scores and common evaluation metrics used in PRS settings.

## Pharmacogenomics

- [PharmCAT](https://github.com/PharmGKB/PharmCAT) Pharmacogenomics clinical annotation tool, a caller for pharmacogene alleles and their corresponding metaboliser phenotypes. Originally intended for clinical applications and single samples, but can be [scaled](https://pharmcat.org/using/Multi-Sample-Analysis/) to datasets of an arbitrary size.

- [PGxPOP](https://github.com/PharmGKB/PGxPOP) Population-scale pharmacogenetic allele and phenotype caller. Allele definitions are based on a 2020 release of the [PharmVar](https://www.pharmvar.org/) database and might need updating.
  
- [PyPGx](https://github.com/sbslee/pypgx) Pharmacogenomic allele and phenotype caller that can be used on various genomic datasets, including next-generation sequencing (NGS), single nucleotide polymorphism (SNP) arrays, and long-read sequencing projects. It can detect and account structural variants if used on NGS data. It supports 87 [pharmacogenes](https://pypgx.readthedocs.io/en/latest/genes.html), but not all have genotype/phenotype mappings.

## GWAS result repositories (preferably with an API)

- [ieugwasr](https://github.com/MRCIEU/ieugwasr) The OpenGWAS database comprises over 50,000 curated, QC'd and harmonised complete GWAS summary datasets and can be queried using an API. See here for documentation on the API itself. This R package is a wrapper to make generic calls to the API.
  
- [GWAScatalog API](https://www.ebi.ac.uk/gwas/summary-statistics/docs/) As of 2024-11-20, the GWAS Catalog contains 7083 publications, 692444 top associations and 96947 full summary statistics.
GWAS Catalog data is currently mapped to Genome Assembly GRCh38.p14 and dbSNP Build 156.

- [GWAS atlas](https://atlas.ctglab.nl) This atlas is a database of publicly available GWAS summary statistics. Each GWAS can be browsed with the manhattan plot, risk loci, MAGMA (i.e. gene-based) results, SNP heritability and genetic correlations with other GWAS in the database.

- [S4 programs](https://github.com/jpt34/S4_programs) S4 programs to calculate PRSs and evaluate them.

## Mendelian randomization result repositories 
- [Multi-ancestry Proteome-Phenome MR atlas](https://broad.io/protein_mr_atlas) Multi-ancestry atlas of protein-phenotype associations in European, African, and East Asian ancestries using MR and colocalization.

- [Proteome-phenome atlas](https://proteome-phenome-atlas.com/) Atlas of protein-phenotype associations in the UK Biobank (2024).

- [Proteome MR](https://www.epigraphdb.org/pqtl/) Mendelian Randomization and sensitivity analyses results for 1,545 proteins on 8 complex diseases in European and African ancestries (2022).
  
- [Proteome MR atlas](https://www.epigraphdb.org/pqtl/) Mendelian Randomization and sensitivity analyses results for 989 proteins and 225 traits in Europeans (2020).

## Online tools

- [gnomAD browser](https://gnomad.broadinstitute.org/) The gnomAD (genome aggregation database) browser is an online tool for querying the gnomAD and ExAC population data

- [Open Targets Platform](https://platform.opentargets.org/) The Open Targets Platform is a comprehensive tool that supports systematic identification and prioritisation of potential therapeutic drug targets.

- [genebass](https://app.genebass.org/) Genebass is a resource of exome-based association statistics, made available to the public. The dataset encompasses 4,529 phenotypes with gene-based and single-variant testing across 394,841 individuals with exome sequence data from the UK Biobank.

- [All by All](https://allbyall.researchallofus.org/) The All by All browser maps known and novel associations between genotypes and phenotypes using data contributed by All of Us Research Program participants as of July 1, 2022. All by All encompasses about 3,400 phenotypes with gene-based and single-variant associations across nearly 250,000 whole genome sequences

- [Bravo](https://bravo.sph.umich.edu/) Variant browser from the Trans-Omics for Precision Medicine (TOPMed) consortium

- [Regeneron Exome browser](https://rgc-research.regeneron.com/me/home) Variant browser from Regeneron from ~983K individuals

 
