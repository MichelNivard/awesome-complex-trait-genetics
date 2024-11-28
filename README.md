# awesome-complex-trait-genetics

[![Awesome](https://awesome.re/badge.svg)](https://awesome.re)

A list of awesome tools for human complex trait genetics, biased/slanted towards post-GWAS analysis.

How can you add your tool? via a pull request, if you dont know what that is, read the contributing guidlines: (https://github.com/MichelNivard/awesome-genetics/blob/main/contributing.md)

**This list may for the basis for a static review**, in which case serious contributors, so those who contribute to the list and to the preamle paragraphs for each section, will be ivited as collaborators.


## genetic architecture

**Sumstats based:**

- [LDSC](https://github.com/bulik/ldsc) ldsc is a python command line tool for estimating heritability and genetic correlation from GWAS summary statistics. ldsc also computes LD Scores. A [Python3 port](https://github.com/belowlab/ldsc) is available.

**Raw data based:**

- [GCTA](https://yanglab.westlake.edu.cn/software/gcta/) GCTA (Genome-wide Complex Trait Analysis) is a software package initially developed to estimate the proportion of phenotypic variance explained by all genome-wide SNPs for a complex trait.

### Univariate models (heritability/poligenicity/stratified/geneset enrichment etc)

- [i-LDSC](https://github.com/lcrawlab/i-LDSC)   interaction-LD score (i-LDSC) regression: Model an additional score that measures the amount of non-additive genetic variation that is tagged by each variant in the data.
- [ACLR](https://github.com/arminschoech/ACLR) Autocorrelation LD regression: a tool to efficiently estimate the autocorrelation of latent effects in large genetic data sets. (WARNING Python 2.7 code)
- [HAMSTA](https://github.com/tszfungc/hamsta) HAMSTA is a python package that estimate heritability explained by local ancestry using summary statistics from admixture mapping studies.
- [MAGMA](https://cncr.nl/research/magma/) MAGMA: Generalized gene-set analysis of GWAS data. 


## Genetic correlation (LD score derivatives/extensions)


- [HDL](https://github.com/zhenin/HDL) High-Definition Likelihood (HDL) is a likelihood-based method for estimating genetic correlation using GWAS summary statistics. Compared to LD Score regression (LDSC), It reduces the variance of a genetic correlation estimate by about 60%.

### Stratified/local genetic correlatons

- [LAVA](https://github.com/josefin-werme/LAVA) LAVA (Local Analysis of [co]Variant Association) is a tool developed for local genetic correlation (rg) analysis.

### Ancestry aware Genetic correlations:

- [s-ldxr](https://github.com/huwenboshi/s-ldxr) S-LDXR is a method to stratify squared trans-ethnic genetic correlation by genomic annotations from GWAS summary statistics.

## Model trait relationships

### Genetic SEM/Factor models

**sumstat based:**

- [GenomicSEM](https://github.com/GenomicSEM/GenomicSEM) R package for Genomic Structural Equation Models. user defined models of the relation between complex traits based on GWAS summary data.
- [GUIDE](https://github.com/daniel-lazarev/GUIDE) Genetic Unmixing by Independent Decomposition (GUIDE), uses ICA to estimate statistically independent latent factors that best express the patterns of association across many traits.
- [FactorGO](https://github.com/mancusolab/FactorGo) FactorGo is a scalable variational factor analysis model that learns pleiotropic factors using GWAS summary statistics!
- [GNA](https://github.com/GenomicNetworkAnalysis/GNA) GNA is an R package for performing network analysis of genetic overlap derived from GWAS summary statistics
- [partialLDSC](https://github.com/GEMINI-multimorbidity/partialLDSC) is an R-package to estimate partial genetic correlations from GWAS summary statistics, and compare them to their unadjusted counterparts, to quantify the contribution of a given confounder in explaining genetic similarity between conditions.

**Raw data based:**

- [Gsens](https://github.com/LeonardFrach/Gsens): Genetically informed sensitivity analysis to estimate role of genetic confounding in phenotypic associations and for causal inference from observation data.


### Two sample Mendelian Randomisation

- [TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR) An R package for performing Mendelian randomization using GWAS summary data.

### MR/Genetic architecture hybrid models

- [lhcMR](https://github.com/LizaDarrous/lhcMR) lhcMR is an R package that performs bi-directional causal estimation between a pair of traits, while accounting for the presence of a potential heritable confounder acting on the pair.
- [CAUSE](https://github.com/jean997/cause) an MR method, Causal Analysis Using Summary Effect Estimates (CAUSE), that accounts for correlated and uncorrelated horizontal pleiotropic effects. 
- [LCV](https://github.com/lukejoconnor/LCV) LCV is an LD score based method for inferring genetically causal relationships using GWAS data.

## Colocalisation/finemapping of causal variants

- [coloc](https://github.com/chr1swallace/coloc) R package to perform genetic colocalisation analysis, test whether two related phenotypes share common genetic causal variant(s) in a given region.
- [fastENLOC](https://github.com/xqwen/fastenloc) This repository contains the software implementation of FastENLOC, which enables integrative genetic association analysis of molecular QTL data and GWAS data.
- [FINEMAP](http://www.christianbenner.com) FINEMAP: Efficient variable selection using summary data from genome-wide association studies. 
- [polyfun](https://github.com/omerwe/polyfun) Package contains **PolyFun** for functionally-informed fine-mapping, **PolyLoc** for polygenic localization of complex trait heritability.
- [OPERA](https://github.com/wuyangf7/OPERA): (omics pleiotropic association) method tests combinatorial pleiotropic associations between multiple molecular phenotypes (eQTL, DNAm) with a complex trait using summary-level data.
- [SuSiEx](https://github.com/getian107/SuSiEx): Cross-population finemapping using summary statistics and LD reference panels.

## gene-level analysis (TWAS)

- [FUSION](https://github.com/gusevlab/fusion_twas) FUSION is a suite of tools for performing transcriptome-wide and regulome-wide association studies (TWAS and RWAS). 
- [FOCUS](https://github.com/mancusolab/ma-focus) FOCUS (Fine-mapping Of CaUsal gene Sets) is software to fine-map transcriptome-wide association study statistics at genomic risk regions



## Simulation

- [GWASBrewer](https://github.com/jean997/GWASBrewer) A flexible tool for simulating realistic GWAS summary statistics for one, or many, traits.

## Genomic data wrangling

- [HAIL](https://github.com/hail-is/hail) Hail is an open-source, general-purpose, Python-based data analysis tool with additional data types and methods for working with genomic data.
- [bigsnpr](https://github.com/privefl/bigsnpr) R package for the analysis of massive SNP arrays, primarily designed for human genetics.
- [ukbrapR](https://github.com/lcpilling/ukbrapR) ukbrapR (phonetically: 'U-K-B-wrapper') is an R package for working in the UK Biobank Research Analysis Platform (RAP). The aim is to make it quicker, easier, and more reproducible.
- [MungeSumstats](https://github.com/Al-Murphy/MungeSumstats) R package designed to facilitate the standardisation of GWAS summary statistics.
- [gwasRtools](https://github.com/lcpilling/gwasRtools) R package to (1) identify loci and independent lead SNPs (using online or local reference panel) and (2) annotate variants with nearest gene from GENCODE database.



