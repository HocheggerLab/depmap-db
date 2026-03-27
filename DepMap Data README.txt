# DepMap Public 25Q2

##############
## Overview ##
##############

This DepMap Release contains new cell models and data from Whole Genome/Exome Sequencing (Copy Number and Mutation), RNA Sequencing (Expression and Fusions), Genome-wide CRISPR knockout screens. Also included are updated metadata and mapping files for information about cell models and data relationships, respectively. Each release may contain improvements to our pipelines that generate this data so you may notice changes from the last release.

###############
## Pipelines ##
###############

### Achilles

This Achilles dataset contains the results of genome-scale CRISPR knockout screens for Achilles (using Avana Cas9 and Humagne-CD Cas12 libraries) and Achilles combined with Sanger's Project SCORE (KY Cas9 library) screens.

The dataset was processed using the following steps:
- Sum raw readcounts by SequenceID and collapse pDNAs by median (removing sequences with fingerprinting issues, total reads below 500000 or less then 25 mean reads per guide)
- Normalize readcounts such that pDNA sequences have the same mode, and all other sequences have median of nonessentials that matches the median of nonessentials in the corresponding pDNA batch
- Remove intergenic/non-targeting controls, Cas9 guides with multiple alignments to a gene, guides with >5 genomic alignments, and sgRNAs with inconsistent signal (see DropReason of guide maps).
- NA sgRNAs with pDNA counts less than one millionth of the pDNA pool
- Calculate the mean reads per guide for each sequence. Sequences with more than 185 mean reads per guide are considered passing.
- Compute the log2-fold-change (LFC) from pDNA counts for each sequence . Collapse to a gene-level by taking the median of guides.
- Calculate the Null-Normalized Median Difference (NNMD) from the LFC using the following equation: (median(essentials) - median(nonessentials)) / MAD(nonessentials).  Sequences must be below a threshold of -1.25 to be considered passing.
- Calculate the fraction of reads originating from the other library than what is annotated. The fraction of reads must be below a threshold of 0.1 to be considered passing.
- Calculate residual LFC for each library and replicate by creating a linear model of the replicate's LFC as a function of the mean LFC for each gene across replicates. Only genes that show high variance in gene effect across screens and libraries are considered when calculating residual LFC.  Specifically, this set of controls was identified by taking the intersection of the top 1250 genes with the highest variance gene effect from each library. This removes the influence of common essential dropout on replicate correlations.
- Remove sequences that do not have a Pearson coefficient > .19 with at least one other replicate sequence for the same screen using the residual LFC. Remove sequences that have a Pearson coefficient > .8 with any other sequence that is not part of the same screen and is not the same lineage. This checks for inconsistent biology between replicates, or unexpected similarity to another screen.
- Calculate the NNMD for each screen after averaging passing sequences. Screen must be below a threshold of -1.25 to be considered passing.
- In each passing replicate, NA vectors if any contained sgRNA overlaps with a SNP in that cell line, either somatic or germline. Recompute LFC after correction.
- Compute the naive gene score by collapsing the passing LFC data to a screen x gene matrix using the median of all replicates for each screen and the median of sgRNAs for each gene.
- Prior to running Chronos, NA apparent outgrowths in readcounts.
- Run Chronos per library-screen type batch to generate screen-level gene effect scores, where screen type batch 3D includes all organoid screens (2DO, 3DO) and screen type batch 2D includes all standard, neurosphere, and anchor screens (2DS, 2DD, 2DN, 3DN, 2DND, 3DND).
-  Apply copy number correction, then scale such that median of common essentials is at -1.0 and median of nonessentials is at 0.
- To correct proximity bias, the median gene effect of each chromosome arm is aligned to be the same across all screens. See Venceti et al. 2024 (https://doi.org/10.1186/s13059-024-03336-1) for more details on this correction.
- Concatenate gene effects from all libraries into a single ScreenGeneEffect matrix.
- Run Chronos jointly on all libraries-screen types to generate model-level gene effect scores which are provided as CRISPRGeneEffectUncorrected, then perform the same corrections as described above. This produces the integrated CRISPRGeneEffect matrix using Chronos's innate batch correction.
- Removal of Confounding Principal Components (RCPC): Preprocess CRISPRGeneEffect matrix by subtracting column mean and standardizing column variance. Perform PCA on preprocessed CRISPRGeneEffect matrix and correlate component loadings with median of common essential genes. Remove component loadings with Pearson correlation higher than 0.2 with common essential dropout. Multiply the truncated component loadings and principal components to revert back to the CRISPRGeneEffect matrix of original dimension. Add column mean and multiply column variance to CRISPRGeneEffect to revert back to the original scale prior to preprocessing. To learn more refer to (www.biorxiv.org/content/10.1101/720243v1).
- Using the CRISPRGeneEffect, identify pan-dependent genes as those for whom 90% of cell lines rank the gene above a given dependency cutoff. The cutoff is determined from the central minimum in a histogram of gene ranks in their 90th percentile least dependent line.
- For each Chronos gene effect score from both the ScreenGeneEffect and CRISPRGeneEffect, infer the probability that the score represents a true dependency. This is done using an EM step until convergence independently in each screen/model. The dependent distribution is given by the list of essential genes. The null distribution is determined from unexpressed gene scores in those cell lines that have expression data available, and from the nonessential gene list in the remainder.
The essential  and nonessential controls used throughout the analysis are the Hart reference nonessentials and the intersection of the Hart and Blomen essentials. See Hart et al., Mol. Syst. Biol, 2014 and Blomen et al., Science, 2015. They are provided with this release as AchillesCommonEssentialControls.csv and AchillesNonessentialControls.csv.
Genes that are exclusively targeted by KY or Humagne library have been removed from CRISPRGeneEffect/Uncorrected, CRISPRGeneDependency, and CRISPRInferredCommonEssential. These genes have shown generalized correlation with other genes driven by batch effect between sequencing technologies, which we are working to address.


### Expression

DepMap expression data for genes and transcripts(TPMs) is quantified from RNAseq files using Salmon version v1.10.0 (Patro et al April 2017). Raw counts are quantified by STAR v2.7.11b. For code and configuration details, see https://github.com/broadinstitute/depmap-omics-rna/tree/main/workflows/quantify_sr_rna. Based on Gencode v38.

### Copy number

DepMap WGS/WES relative copy number data is generated by running the GATK copy number pipeline aligned to hg38. Tutorials and descriptions of this method can be found here https://software.broadinstitute.org/gatk/documentation/article?id=11682, https://software.broadinstitute.org/gatk/documentation/article?id=11683.
Absolute copy number data from WGS/WES is generated using PureCN (https://scfbm.biomedcentral.com/articles/10.1186/s13029-016-0060-z).

### Mutations

DepMap mutation calls are generated using Mutect2 and annotated and filtered downstream. Variants are aligned to hg38. Detailed documentation can be found here https://storage.googleapis.com/shared-portal-files/Tools/25Q2_Mutation_Pipeline_Documentation.pdf.

### fusionss

DepMap generates RNAseq based fusion calls using 2.5.0 (Uhrig et al March 2021). For code and configuration details, see https://github.com/broadinstitute/depmap-omics-rna/tree/main/workflows/call_sr_rna_fusions. Based on Gencode v38.

###########
## Files ##
###########

### ScreenSequenceMap.csv

Pipeline: Achilles

**_Pre-Chronos_**

Mapping of SequenceIDs to ScreenID and related info.

**Columns:**

- SequenceID
- ScreenID
- ModelConditionID
- ModelID
- ScreenType: 2DS = 2D standard
- Library
- Days
- pDNABatch
- PassesQC
- ExcludeFromCRISPRCombined

### README.txt

README containing descriptions of each file

### ScreenGeneEffect.csv

Pipeline: Achilles

**_Post-Chronos_**

Gene effect estimates for all screens Chronos processed by library, copy number corrected, scaled, screen quality corrected then concatenated.

**Rows:**

- ScreenID

**Columns:**

- Gene

### ScreenGeneEffectUncorrected.csv

Pipeline: Achilles

**_Post-Chronos_**

Gene effect estimates for all screens Chronos processed by library, then concatenated. No copy number correction or scaling.

**Rows:**

- ScreenID

**Columns:**

- Gene

### ScreenNaiveGeneScore.csv

Pipeline: Achilles

**_Post-Chronos_**

LFC collapsed by mean of sequences and median of guides, computed per library-screen type then concatenated.

**Rows:**

- ScreenID

**Columns:**

- Gene

### ScreenGeneDependency.csv

Pipeline: Achilles

**_Post-Chronos_**

Gene dependency probability estimates for all screens Chronos processed by library-screen type, then concatenated.

**Rows:**

- ScreenID

**Columns:**

- Gene

### CRISPRScreenMap.csv

Pipeline: Achilles

**_Post-Chronos_**

Map from ModelID to all ScreenIDs combined to make up a given model's data in the CRISPRGeneEffect matrix.

**Columns:**

- ModelID
- ScreenID

### CRISPRGeneEffect.csv

Pipeline: Achilles

**_Post-Chronos_**

Gene effect estimates for all models, integrated using Chronos.  Copy number corrected, scaled, and screen quality corrected.

**Rows:**

- ModelID

**Columns:**

- Gene

### CRISPRGeneEffectUncorrected.csv

Pipeline: Achilles

**_Post-Chronos_**

Gene effect estimates for all models, integrated using Chronos. No copy number correction or scaling.

**Rows:**

- ModelID

**Columns:**

- Gene

### CRISPRInferredLibraryEffect.csv

Pipeline: Achilles

**_Post-Chronos_**

The estimates for the library batch effects identified by Chronos.

**Rows:**

- LibraryBatch

**Columns:**

- Gene

### CRISPRInitialOffset.csv

Pipeline: Achilles

**_Post-Chronos_**

Estimated log fold pDNA error for each sgrna in each library, as identified by Chronos.

**Rows:**

- sgRNA

**Columns:**

- pDNABatch


### CRISPRGeneDependency.csv

Pipeline: Achilles

**_Post-Chronos_**

Gene dependency probability estimates for all models in the integrated gene effect.

**Rows:**

- ModelID

**Columns:**

- Gene


### CRISPRInferredGuideEfficacy.csv

Pipeline: Achilles

**_Post-Chronos_**

The estimates for the efficacies of all reagents in the different libraries-screen types, computed from the Chronos runs.

**Columns:**

- sgRNA
- Efficacy

### CRISPRInferredModelGrowthRate.csv

Pipeline: Achilles

**_Post-Chronos_**

The estimates for the growth rate of all models in the different libraries-screen types, computed from the Chronos runs.

**Columns:**

- ScreenID
- Achilles-Avana-2D
- Achilles-Humagne-CD-2D
- Achilles-Humagne-CD-3D
- Project-Score-KY


### CRISPRInferredModelEfficacy.csv

Pipeline: Achilles

**_Post-Chronos_**

The estimates for the efficacy of all models in the different libraries-screen types, computed from the Chronos runs.

**Columns:**

- ModelID
- Achilles-Avana-2D
- Achilles-Humagne-CD-2D
- Achilles-Humagne-CD-3D
- Project-Score-KY


### CRISPRInferredSequenceOverdispersion.csv

Pipeline: Achilles

**_Post-Chronos_**

List of genes identified as dependencies across all lines. Each entry is separated by a newline.

### CRISPRInferredCommonEssentials.csv

Pipeline: Achilles

**_Post-Chronos_**

List of genes identified as dependencies across all lines. Each entry is separated by a newline.

### AchillesCommonEssentialControls.csv

Pipeline: Achilles

**_Pre-Chronos_**

List of genes used as positive controls, intersection of Biomen (2014) and Hart (2015) essentials. Each entry is separated by a newline.

The scores of these genes are used as the dependent distribution for inferring dependency probability

### AchillesNonessentialControls.csv

Pipeline: Achilles

**_Pre-Chronos_**

List of genes used as negative controls (Hart (2014) nonessentials). Each entry is separated by a newline.

### AvanaRawReadcounts.csv

Pipeline: Achilles

**_Pre-Chronos_**

Summed guide-level read counts for each sequence screened with the Avana Cas9 library.

**Rows:**

- sgRNA

**Columns:**

- SequenceID

### HumagneRawReadcounts.csv

Pipeline: Achilles

**_Pre-Chronos_**

Summed guide-level read counts for each sequence screened with the Humagne-CD Cas12 library.

**Rows:**

- sgRNA

**Columns:**

- SequenceID

### KYRawReadcounts.csv

Pipeline: Achilles

**_Pre-Chronos_**

Summed guide-level read counts for each sequence screened with the Sanger's KY Cas9 library.

**Rows:**

- sgRNA

**Columns:**

- SequenceID

### AvanaLogfoldChange.csv

Pipeline: Achilles

**_Pre-Chronos_**

Log2-fold-change from pDNA counts for each sequence screened with the Avana Cas9 Library.

**Rows:**

- sgRNA

**Columns:**

- SequenceID

### HumagneLogfoldChange.csv

Pipeline: Achilles

**_Pre-Chronos_**

Log2-fold-change from pDNA counts for each sequence screened with the Humagne-CD Library.

**Rows:**

- sgRNA

**Columns:**

- SequenceID

### KYLogfoldChange.csv

Pipeline: Achilles

**_Pre-Chronos_**

Log2-fold-change from pDNA counts for each sequence screened with the Sanger's KY Cas9 Library.

**Rows:**

- sgRNA

**Columns:**

- SequenceID

### AvanaGuideMap.csv

Pipeline: Achilles

**_Pre-Chronos_**

Mapping of sgRNAs to Genes in the Avana Cas9 library.

**Columns:**

- sgRNA: guide in vector
- GenomeAlignment: alignment to hg38
- Gene: HUGO (entrez)
- nAlignments: total number of alignments for a given sgRNA
- UsedByChronos: boolean indicating if sgRNA was included in Chronos analysis
- DropReason: why a guide was removed prior to Chronos analysis


### HumagneGuideMap.csv

Pipeline: Achilles

**_Pre-Chronos_**

Mapping of sgRNAs to Genes in the Humagne-CD Cas12 library.

**Columns:**

- sgRNA: first guide in vector
- GenomeAlignment: alignment of first guide to hg38
- sgRNA2: second guide in vector
- GenomeAlignment2: alignment of second guide to hg38
- Set: which half of the library the vector is part of (C, D, or Both)
- Gene: HUGO (entrez)
- nAlignments: total number of alignments for a given sgRNA
- UsedByChronos: boolean indicating if sgRNA was included in Chronos analysis
- DropReason: why a guide was removed prior to Chronos analysis


### KYGuideMap.csv

Pipeline: Achilles

**_Pre-Chronos_**

Mapping of sgRNAs to Genes in the Sanger's KY Cas9 library.

**Columns:**

- sgRNA: guide in vector
- GenomeAlignment: alignment to hg38
- Gene: HUGO (entrez)
- nAlignments: total number of alignments for a given sgRNA
- UsedByChronos: boolean indicating if sgRNA was included in Chronos analysis
- DropReason: why a guide was removed prior to Chronos analysis


### AchillesScreenQCReport.csv

Pipeline: Achilles

**_Pre-Chronos_**

Screen-level quality control metrics.

**Columns:**

- ScreenID
- ScreenNNMD: null-normalized median difference (threshold of -1.25)
- ScreenROCAUC: area under the Receiver Operating Characteristic curve
for essential (positive) vs nonessentials (negative) controls
- ScreenFPR: false positive rate, computed as fraction of nonessentials in
the 15th percentile most depleted genes
- ScreenMeanEssentialDepletion: mean LFC of essentials
- ScreenMeanNonessentialDepletion: mean LFC of nonessentials
- ScreenSTDEssentials: standard deviation of essential gene LFCs
- ScreenSTDNonessentials: standard deviation of nonessential gene LFCs
- PassesQC: boolean indicating if screen passes QC thresholds
- CanInclude: boolean indicating if screen can be included in dataset
- QCStatus: string describing QC status
- HasCopyNumber: boolean indicating if screen has copy number data
- ModelConditionID
- ModelID
- Library
- ScreenType
- CasActivity: percentage of cells remaining GFP negative on days 12-14 of cas9
  activity assay as measured by FACs
- ScreenDoublingTime: hours for cells to double
- nPassingSequences: number of sequences that pass sequence-level QC
- nIncludedSequences: number of sequences that can be included


### AchillesSequenceQCReport.csv

Pipeline: Achilles

**_Pre-Chronos_**

Sequence-level quality control metrics.

**Columns:**

- SequenceID
- SequenceTotalReads: total number of raw readcounts
- SequenceMeanReadsPerGuide: mean reads per guide following normalization (threshold of 185)
- SequenceNNMD: null-normalized median difference (threshold of -1.25)
- PassesQC: boolean indicating if sequence passes QC thresholds
- CanInclude: boolean indicating if sequence can be included in dataset
- QCStatus: string describing QC status
- SequenceMaxCorr: max correlation of residual LFC with other
sequences from same screen (threshold of 0.19)
- SequenceFracReadsFromOtherLibrarySubset: the fraction of reads originating from a library other than what is annotated (threshold of 0.1)
- UnexpectedHighCorrelation: max correlation of residual LFC with other
sequences that are not from the same screen (threshold of 0.8)
- UnexpectedHighCorrelationPartners: list of sequences that are highly correlated
with sequences that are not from the same screen (threshold of 0.8)
- ConfirmedNoSwap: confirms a sample swap did not occur and overrides QC failure due to high correlation with random sequences reported in column “UnexpectedHighCorrelationPartners"
- ScreenID
- ScreenPassesQC: boolean indicating if corresponding screen passes QC thresholds


### AchillesHighVarianceGeneControls.csv

Pipeline: Achilles

**_Pre-Chronos_**

List of genes with variable gene effects across models and libraries. Used for sequence correlation in QC.

### OmicsCNSegmentsProfileWGS.csv

Pipeline: Copy number

Segment level copy number data derived from WGS.

More information on the DepMap Omics copy number processing pipeline is available at <https://github.com/broadinstitute/depmap_omics> and <https://github.com/broadinstitute/depmap-omics-wgs/tree/main/workflows/call_cnvs>.

**Columns:**

- Profile ID
- Chromosome
- Start (bp start of the segment)
- End (bp end of the segment)
- Num_Probes (the number of targeting probes that make up this segment)
- Segment_Mean (relative copy ratio for that segment)
- Amplification status (+,-,0,.)

### OmicsCNGeneWGS.csv

Pipeline: Copy number

Gene-level copy number data inferred from WGS data only

Additional copy number datasets are available for download as part of the full DepMap Data Release.

More information on the DepMap Omics copy number processing pipeline is available at <https://github.com/broadinstitute/depmap_omics> and <https://github.com/broadinstitute/depmap-omics-wgs/tree/main/workflows/call_cnvs>.

### PortalOmicsCNGeneLog2.csv

A log2(x+1) transformed version of OmicsCNGene.csv for use in the portal's Data Explorer and custom analysis tools.

### OmicsSomaticMutations.csv

Pipeline: Mutations

MAF-like file containing information on all the somatic point mutations and indels called in the DepMap cell lines.
The calls are generated from Mutect2.

Additional processed mutation matrices containing genotyped mutation calls are available for download as part of the full DepMap Data Release.

**Columns:**

- Chrom
- Pos
- Ref
- Alt
- AF
- DP
- RefCount
- AltCount
- GT
- PS
- VariantType
- VariantInfo
- DNAChange
- ProteinChange
- HugoSymbol
- EnsemblGeneID
- EnsemblFeatureID
- HgncName
- HgncFamily
- UniprotID
- DbsnpRsID
- GcContent
- LofGeneName
- LofGeneId
- LofNumberOfTranscriptsInGene
- LofPercentOfTranscriptsAffected
- NMD
- MolecularConsequence
- VepImpact
- VepBiotype
- VepHgncID
- VepExistingVariation
- VepManeSelect
- VepENSP
- VepSwissprot
- Sift
- Polyphen
- GnomadeAF
- GnomadgAF
- VepClinSig
- VepSomatic
- VepPliGeneValue
- VepLofTool
- OncogeneHighImpact
- TumorSuppressorHighImpact
- TranscriptLikelyLof
- Brca1FuncScore
- CivicID
- CivicDescription
- CivicScore
- LikelyLoF
- HessDriver
- HessSignature
- RevelScore
- PharmgkbId
- DidaID
- DidaName
- GwasDisease
- GwasPmID
- GtexGene
- ProveanPrediction
- AMClass
- AMPathogenicity
- Hotspot
- Rescue
- EntrezGeneID
- ModelID

For details, see https://storage.googleapis.com/shared-portal-files/Tools/25Q2_Mutation_Pipeline_Documentation.pdf


### OmicsSomaticMutationsProfile.csv

Pipeline: Mutations

MAF-like formatted file containing information on all the somatic point mutations and indels called in the DepMap cell lines.

Calls are generated from Mutect2.

Additional processed mutation matrices containing genotyped mutation calls are available for download as part of the full DepMap Data Release.

Samples are identified by Profile IDs.

**Columns:**

- Chrom
- Pos
- Ref
- Alt
- AF
- DP
- RefCount
- AltCount
- GT
- PS
- VariantType
- VariantInfo
- DNAChange
- ProteinChange
- HugoSymbol
- EnsemblGeneID
- EnsemblFeatureID
- HgncName
- HgncFamily
- UniprotID
- DbsnpRsID
- GcContent
- LofGeneName
- LofGeneId
- LofNumberOfTranscriptsInGene
- LofPercentOfTranscriptsAffected
- NMD
- MolecularConsequence
- VepImpact
- VepBiotype
- VepHgncID
- VepExistingVariation
- VepManeSelect
- VepENSP
- VepSwissprot
- Sift
- Polyphen
- GnomadeAF
- GnomadgAF
- VepClinSig
- VepSomatic
- VepPliGeneValue
- VepLofTool
- OncogeneHighImpact
- TumorSuppressorHighImpact
- TranscriptLikelyLof
- Brca1FuncScore
- CivicID
- CivicDescription
- CivicScore
- LikelyLoF
- HessDriver
- HessSignature
- RevelScore
- PharmgkbId
- DidaID
- DidaName
- GwasDisease
- GwasPmID
- GtexGene
- ProveanPrediction
- AMClass
- AMPathogenicity
- Hotspot
- Rescue
- EntrezGeneID
- ProfileID

For details, see https://storage.googleapis.com/shared-portal-files/Tools/25Q2_Mutation_Pipeline_Documentation.pdf


### OmicsSomaticMutationsMAFProfile.maf

Pipeline: Mutations

MAF file containing information on all the somatic point mutations and indels called in the DepMap cell lines.

Calls are generated from Mutect2.

Additional processed mutation matrices containing genotyped mutation calls are available for download as part of the full DepMap Data Release.

This file contains the same variants as OmicsSomaticMutationsProfile.csv, but follows the standard MAF format suitable for downstream analysis tools such as maftools.

Samples are identified by Profile IDs.

### OmicsSomaticMutationsMatrixHotspot.csv

Pipeline: Mutations

Genotyped matrix determining for each cell line whether each gene has at least one hot spot mutation.

A variant is considered a hot spot if it's present in one of the following: Hess et al. 2019 paper, OncoKB hotspot, COSMIC mutation significance tier 1, TERT promoter mutations (C228T or C250T), mutations in the polypyrimidine track in intron 13 of MET.

(0 == no mutation; If there is one or more hot spot mutations in the same gene for the same cell line, the allele frequencies are summed, and if the sum is greater than 0.95, a value of 2 is assigned and if not, a value of 1 is assigned.)

### OmicsSomaticMutationsMatrixDamaging.csv

Pipeline: Mutations

Genotyped matrix determining for each cell line whether each gene has at least one damaging mutation.

A variant is considered a damaging mutation if LikelyLoF == True.
(0 == no mutation; If there is one or more damaging mutations in the same gene for the same cell line, the allele frequencies are summed, and if the sum is greater than 0.95, a value of 2 is assigned and if not, a value of 1 is assigned.)

### OmicsProfiles.csv

Omics metadata and ID mapping information for files indexed by Profile ID.

**Columns:**

- ModelID: "ACH-*" Model IDs
- StrippedCellLineName
- DepMapCode
- Lineage
- ModelConditionID
- SourceModelCondition
- CellFormat
- GrowthMedia
- GrowthPattern
- ProfileID
- Product
- SequencingID: unique identifier of this table
- DataType
- Stranded
- SharedToDbGaP
- is_default_entry: whether or not each sequencing is selected to represent the model in model-level datasets


### OmicsDefaultModelConditionProfiles.csv

**THIS MAPPING FILE WILL NO LONGER BE RELEASED AFTER 25Q2**

Indicates which profiles are selected for each model condition in Achilles postprocessing.

**Columns:**

- ModelConditionID
- ProfileID
- ProfileType (DNA/RNA)

### OmicsGuideMutationsBinaryKY.csv

Binary matrix indicating whether there are mutations in guide locations from the KY library.

KY guide library (same as the one used in project Score) can be accessed from https://score.depmap.sanger.ac.uk/downloads.

**Columns:**

- Chrom
- Start
- End
- sgRNA
- ModelID

### OmicsGuideMutationsBinaryHumagne.csv

Binary matrix indicating whether there are mutations in guide locations from the Humangne library.

Humagne guide library can be accessed from Addgene (https://www.addgene.org/pooled-library/broadgpp-human-knockout-humagne/).

**Columns:**

- Chrom
- Start
- End
- sgRNA
- ModelID

### OmicsGuideMutationsBinaryAvana.csv

Binary matrix indicating whether there are mutations in guide locations from the Avana library.

Avana guide library can be accessed from AvanaGuideMap.csv.

**Columns:**

- Chrom
- Start
- End
- sgRNA
- ModelID

### OmicsGlobalSignatures.csv

Model-level genomic signatures extracted from WES and WGS data.

Signatures include Ploidy, Chromosomal Instability (CIN), Whole Genome Doubling (WGD), Loss of Heterozygosity Fraction (LoH Fraction), and Microsatellite Instability (MSI Score).

Ploidy, CIN, WGD, and LoH Fraction are generated using PureCN ( https://github.com/lima1/PureCN ).

MSI Score is generated by MSIsensor2 ( https://github.com/niu-lab/msisensor2 ).

Aneuploidy score is implemented based on Ben-David 2021 paper ( https://www.nature.com/articles/s41586-020-03114-6 )

**Rows:**

- ModelID

**Columns:**

- MSIScore (float)
- Ploidy (float)
- CIN (float)
- WGD (binary, 1 indicates presence of WGD, 0 otherwise)
- LoHFraction (float)
- Aneuploidy (int)

### OmicsGlobalSignaturesProfile.csv

Profile-level genomic signatures extracted from WES and WGS data.

Signatures include Ploidy, Chromosomal Instability (CIN), Whole Genome Doubling (WGD), Loss of Heterozygosity Fraction (LoH Fraction), and Microsatellite Instability (MSI Score).

Ploidy, CIN, WGD, and LoH Fraction are generated using PureCN (https://github.com/lima1/PureCN).

MSI Score is generated by MSIsensor2 (https://github.com/niu-lab/msisensor2).

Aneuploidy score is implemented based on Ben-David 2021 paper (https://www.nature.com/articles/s41586-020-03114-6)

**Rows:**

- ProfileID

**Columns:**

- MSIScore (float)
- Ploidy (float)
- CIN (float)
- WGD (binary, 1 indicates presence of WGS, 0 otherwise)
- LoHFraction (float)
- Aneuploidy (int)

### Model.csv

Metadata describing all cancer models/cell lines which are referenced by a dataset contained within the DepMap portal.

- ModelID: Unique identifier for the model

- PatientID: Unique identifier for models derived from the same tissue sample

- CellLineName: Commonly used cell line name

- StrippedCellLineName: Commonly used cell line name without characters or spaces

- DepmapModelType: Abbreviated ID for model type. For most cancer models, this field is from Oncotree, information for other disease types are generated by DepMap

- OncotreeLineage: Lineage of model. For cancer models, this field is from Oncotree, information for other disease types are generated by DepMap

- OncotreePrimaryDisease: Primary disease of model. For cancer models, this field is from Oncotree, information for other disease types are generated by DepMap

- OncotreeSubtype: Subtype of model. For cancer models, this field is from Oncotree, information for other disease types are generated by DepMap

- OncotreeCode: For cancer models, this field is based on Oncotree. For some models for which no corresponding code exists, this field is left blank

- PatientSubtypeFeatures: Aggregated features known for the patient tumor

- RRID: Cellosaurus ID

- Age: Age at time of sampling

- AgeCategory: Age category at time of sampling (Adult, Pediatric, Fetus, Unknown)

- Sex: Sex at time of sampling (Female, Male, Unknown)

- PatientRace: Patient/clinical indicated race (not derived)

- PrimaryOrMetastasis: Site of the primary tumor where cancer originated from (Primary, Metastatic, Recurrence, Other, Unknown)

- SampleCollectionSite: Site of tissue sample collection

- SourceType: Indicates where model was onboarded from (Commercial, Academic lab, Other)

- SourceDetail: Details on where model was onboarded from

- CatalogNumber: Catalog number of cell model, if commercial

- ModelType: Type of model at onboarding (Cell Line or Organoid)

- TissueOrigin: Indicates tissue model was derived from (Human, Mouse, Other)

- ModelDerivationMaterial: Indicates what material a model was derived from (Fresh tissue, PDX, Other)

- ModelTreatment: Indicates a virus used to transform a cell line, if any

- PatientTreatmentStatus: Indicates if sample was collected before, during, or after the patient's cancer treatment (Pre-treatment, Active treatment, Post-treatment, Unknown)

- PatientTreatmentType: Type of treatment patient received prior to, or at the time of, sampling (e.g. chemotherapy, immunotherapy, etc.), if known

- PatientTreatmentDetails: Details about patient treatment

- Stage: Stage of patient tumor

- StagingSystem: Classification system used to categorize disease stage (e.g. AJCC Pathologic Stage), if known

- PatientTumorGrade: Grade (or other marker of proliferation) of the patient tumor, if known

- PatientTreatmentResponse: Any response to treatment, if known

- GrowthPattern: Format model onboarded in (Adherent, Suspension, Dome, Spheroid, Mixed, Unknown)

- OnboardedMedia: Description of onboarding media

- FormulationID: The unique identifier of the onboarding media

- SerumFreeMedia: Indicates a non-serum based media (<1% serum)

- PlateCoating: Coating on plate model onboarded in (Laminin, Matrigel, Collagen, None)

- EngineeredModel: Indicates if model was engineered (genetic knockout, genetic knock down, cultured to resistance, other)

- EngineeredModelDetails: Detailed information for genetic knockdown/out models

- CulturedResistanceDrug: Drug of resistance used for cultured to resistance models

- PublicComments: Comments released to portals

- CCLEName: CCLE name for the cell line

- HCMIID: Identifier models available through the Human Cancer Models Initiative (HCMI)

- PediatricModelType: Indicates if this model represents a pediatric cancer subtype

- ModelAvailableInDbgap: Indicates whether genomics data for a model is approved for data sharing to a specific repository. Refer to OmicsProfiles.csv for specific Omics Profile data available

- ModelSubtypeFeatures: Curated list of confirmed molecular features seen in the model

- WTSIMasterCellID: WTSI ID

- SangerModelID: Sanger ID

- COSMICID: Cosmic ID

- ModelIDAlias: Prior known Model IDs for a given cell line



### ModelCondition.csv

The condition(s) under which the model was assayed.

- ModelConditionID: Unique identifier for each model condition.

- ModelID: Unique identifier for each model (same ID as in model table)

- ParentModelConditionID: ID of parental model condition for new model conditions derived from other model conditions

- DataSource: Site where source data was generated (e.g. Broad, etc.)

- CellFormat: Format the cell line is being grown in (Adherent, Suspension, Dome, Spheroid, Mixed, Unknown). This can differ from OnboardedGrowthPattern.

- Morphology: Description of morphological features of the model in a particular growth condition

- PassageNumber: Number of cell line passages (<5, 6-10, 10+)

- GrowthMedia: Media condition was grown in at the model condition level

- FormulationID: Media formulation

- PlateCoating: Substrate used to coat plates (Laminin, Matrigel, None, Unknown)

- SerumFreeMedia: Indicates a non-serum based media (<1% serum)

- PrescreenTreatmentDrug: Drug used in pretreatment prior to screening

- PrescreenTreatmentDrugDays: Duration of drug treatment prior to screening

- AnchorDrug: Name of drug treatment used in anchor screen

- AnchorDrugConcentration: Concentration of drug treatment used in anchor screen (with units)

- AnchorDaysWithDrug: Range of days the model was treated with drug in anchor screens



### OmicsMicrosatelliteRepeatsProfile.csv

Weighted mean of the number of motif repeats predicted by [MSIsensor2](https://github.com/niu-lab/msisensor2) at each microsatellite site. This data is upstream of MSIscore in OmicsSignatures.

Indexed by ProfileID

### Gene.csv

This file contains the metadata about all genes which are loaded into the DepMap portal's database. Most pipelines will now use this table for mapping between various gene and protein identifiers as we work to get them all switched over. This file will be used to ensure all files have consistent gene symbols and other gene identifiers.

This table was taken from the quarterly export the gene information published by HGNC and downloaded from https://www.genenames.org/download/archive/quarterly/tsv/


### PortalCompounds.csv

This file contains the metadata on all compounds which is loaded into the DepMap portal along with the set of IDs used in the original drug screens data files which transformed and loaded into the portal. Several columns may be blank for compounds used by older drug screens. As new drugs are added we're curating metadata, but still need to backfill these columns for past drug screens.

Columns:
* CompoundID: The internal DepMap portal compound ID
* CompoundName: The canonical name of the compound which will be shown as the default name for the compound
* GeneSymbolOfTargets: A ";" delimited list of gene symbols which are targets of this compound
* TargetOrMechanism: A brief description of the targets and/or mechanism
* Synonyms: A ";" delimited list of alternative names of this compound
* SampleIDs: A ";" delimited list of IDs used to reference this compound in different drug datasets. Each ID is prefixed with an abbreviation representing which data source used the ID
* ChEMBLID: ChEMBLID for this compound taken from PubChem (only populated for those compounds which have PubChemCID populated)
* SMILES: Canonical SMILES string describing molecular structure of the compound taken from PubChem (only populated for those compounds which have PubChemCID populated)
* InChIKey: InChIKey taken from PubChem (only populated for those compounds which have PubChemCID populated)
* DoseUnit: Either "uM" or "mg/mL"
* PubChemCID: PubChem Compound ID



### OmicsCNGeneWES.csv

Gene-level copy number data inferred from WES using the GATK CNV pipeline, for cell lines that do not have WGS profiles. DepMap is no longer profiling cell lines using WES, so this dataset is not actively maintained.

### OmicsCNSegmentsProfileWES.csv

Segment level copy number data derived from WES using the GATK CNV pipeline, for cell lines that do not have WGS profiles. DepMap is no longer profiling cell lines using WES, so this dataset is not actively maintained.

**Columns:**

- Profile ID
- Chromosome
- Start (bp start of the segment)
- End (bp end of the segment)
- Num_Probes (the number of targeting probes that make up this segment)
- Segment_Mean (relative copy ratio for that segment)
- Amplification status (+,-,0,.)

### OmicsCNGeneProfileWGS.csv

Gene-level copy number data inferred from WGS, indexed by ProfileIDs.


Additional copy number datasets are available for download as part of the full DepMap Data Release.

More information on the DepMap Omics copy number processing pipeline is available at <https://github.com/broadinstitute/depmap_omics> and <https://github.com/broadinstitute/depmap-omics-wgs/tree/main/workflows/call_cnvs>.

### OmicsCNGeneProfileWES.csv

Gene-level copy number data inferred from WES using the GATK CNV pipeline, indexed by ProfileIDs. DepMap is no longer profiling cell lines using WES, so this dataset is not actively maintained.

### OmicsExpressionEffectiveLengthHumanAllGenes.csv

Effective length at gene level derived from unstranded RNA-seq for all genes in human

One row per profile
Columns are ProfileID, is_default_entry, ModelID and then a column per gene


### OmicsExpressionEffectiveLengthHumanAllGenesStranded.csv

Effective length at gene level derived from strand-specific RNA-seq for all genes in human

One row per profile
Columns are ProfileID, is_default_entry, ModelID and then a column per gene


### OmicsExpressionEffectiveLengthHumanProteinCodingGenes.csv

Effective length at gene level derived from unstranded RNA-seq for protein-coding genes in human

One row per profile
Columns are ProfileID, is_default_entry, ModelID and then a column per gene


### OmicsExpressionEffectiveLengthHumanProteinCodingGenesStranded.csv

Effective length at gene level derived from strand-specific RNA-seq for protein-coding genes in human

One row per profile
Columns are ProfileID, is_default_entry, ModelID and then a column per gene


### OmicsExpressionExpectedCountHumanAllGenes.csv

Expected count at gene level derived from unstranded RNA-seq for all genes in human

One row per profile
Columns are ProfileID, is_default_entry, ModelID and then a column per gene


### OmicsExpressionExpectedCountHumanAllGenesStranded.csv

Expected count at gene level derived from strand-specific RNA-seq for all genes in human

One row per profile
Columns are ProfileID, is_default_entry, ModelID and then a column per gene


### OmicsExpressionExpectedCountHumanProteinCodingGenes.csv

Expected count at gene level derived from unstranded RNA-seq for protein-coding genes in human

One row per profile
Columns are ProfileID, is_default_entry, ModelID and then a column per gene


### OmicsExpressionExpectedCountHumanProteinCodingGenesStranded.csv

Expected count at gene level derived from strand-specific RNA-seq for protein-coding genes in human

One row per profile
Columns are ProfileID, is_default_entry, ModelID and then a column per gene


### OmicsExpressionProteinCodingGenesTPMLogp1.csv

Model-level TPMs derived from Salmon v1.10.0 (Patro et al 2017)
Rows: Model IDs
Columns: Gene names

### OmicsExpressionProteinCodingGenesTPMLogp1Stranded.csv

Model level TPMs derived from Salmon in a strand-aware manner where applicable
Rows are model IDs
Columns are gene names

### OmicsExpressionRawReadCountHumanAllGenes.csv

Raw read count at gene level derived from unstranded RNA-seq for all genes in human

One row per profile




### OmicsExpressionRawReadCountHumanAllGenesStranded.csv

Raw read count at gene level derived from strand-specific RNA-seq for all genes in human

One row per profile
Columns are ProfileID, is_default_entry, ModelID and then a column per gene


### OmicsExpressionRawReadCountHumanProteinCodingGenes.csv

Raw read count at gene level derived from unstranded RNA-seq for protein-coding genes in human

One row per profile
Columns are ProfileID, is_default_entry, ModelID and then a column per gene


### OmicsExpressionRawReadCountHumanProteinCodingGenesStranded.csv

Raw read count at gene level derived from strand-specific RNA-seq for protein-coding genes in human

One row per profile
Columns are ProfileID, is_default_entry, ModelID and then a column per gene


### OmicsExpressionTPMLogp1HumanAllGenes.csv

Log-transformed TPM (Transcripts Per Million) values at gene level derived from unstranded RNA-seq for all genes in human

One row per profile
Columns are ProfileID, is_default_entry, ModelID and then a column per gene


### OmicsExpressionTPMLogp1HumanAllGenesStranded.csv

Log-transformed TPM (Transcripts Per Million) values at gene level derived from strand-specific RNA-seq for all genes in human

One row per profile
Columns are ProfileID, is_default_entry, ModelID and then a column per gene


### OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv

Log-transformed TPM (Transcripts Per Million) values at gene level derived from unstranded RNA-seq for protein-coding genes in human

One row per profile
Columns are ProfileID, is_default_entry, ModelID and then a column per gene


### OmicsExpressionTPMLogp1HumanProteinCodingGenesStranded.csv

Log-transformed TPM (Transcripts Per Million) values at gene level derived from strand-specific RNA-seq for protein-coding genes in human

One row per profile
Columns are ProfileID, is_default_entry, ModelID and then a column per gene


### OmicsExpressionTranscriptEffectiveLengthHumanAllGenes.csv

Effective length at transcript level derived from unstranded RNA-seq for all genes in human

One row per profile
Columns are ProfileID, is_default_entry, ModelID and then a column per transcript


### OmicsExpressionTranscriptEffectiveLengthHumanAllGenesStranded.csv

Effective length at transcript level derived from strand-specific RNA-seq for all genes in human

One row per profile
Columns are ProfileID, is_default_entry, ModelID and then a column per transcript


### OmicsExpressionTranscriptExpectedCountHumanAllGenes.csv

Expected count at transcript level derived from unstranded RNA-seq for all genes in human

One row per profile
Columns are ProfileID, is_default_entry, ModelID and then a column per transcript


### OmicsExpressionTranscriptExpectedCountHumanAllGenesStranded.csv

Expected count at transcript level derived from strand-specific RNA-seq for all genes in human

One row per profile
Columns are ProfileID, is_default_entry, ModelID and then a column per transcript


### OmicsExpressionTranscriptTPMLogp1HumanAllGenes.csv

Log-transformed TPM (Transcripts Per Million) values at transcript level derived from unstranded RNA-seq for all genes in human

One row per profile
Columns are ProfileID, is_default_entry, ModelID and then a column per transcript


### OmicsExpressionTranscriptTPMLogp1HumanAllGenesStranded.csv

Log-transformed TPM (Transcripts Per Million) values at transcript level derived from strand-specific RNA-seq for all genes in human

One row per profile
Columns are ProfileID, is_default_entry, ModelID and then a column per transcript


### OmicsFusionFiltered.csv

Gene fusions derived from RNA-seq. Fusions are collapsed at the gene level, averaging over different breakpoints within each gene. Only high and medium confidence fusions are output in this file.

### OmicsFusionFilteredSupplementary.csv

Gene fusions derived from RNA-seq. This is the full data set of gene fusions where each row represents an individual breakpoint.
