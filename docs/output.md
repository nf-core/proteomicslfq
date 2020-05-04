# nf-core/proteomicslfq: Output

This document describes the output produced by the pipeline.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

1. (optional) Conversion of spectra data to indexedMzML: Using ThermoRawFileParser if Thermo Raw or using OpenMS' FileConverter if just an index is missing
1. (optional) Decoy database generation for the provided DB (fasta) with OpenMS
1. Database search with either MSGF+ or Comet through OpenMS adapters
1. Re-mapping potentially identified peptides to the database for consistency and error-checking (using OpenMS' PeptideIndexer)
1. (Intermediate score switching steps to use appropriate scores for the next step)
1. PSM rescoring using PSMFeatureExtractor and Percolator or a PeptideProphet-like distribution fitting approach in OpenMS
1. (Intermediate score switching steps to use appropriate scores for the next step)
1. PSM/Peptide-level FDR filtering
1. Protein inference and labelfree quantification based on MS1 feature detection, alignment and integration with OpenMS' ProteomicsLFQ

## Output

Output is by default written to the $NXF_WORKSPACE/results folder. You can change that with TODO
The output consists of the following folders:

results

* ids
  * [${infile}\*.idXML](#identifications)
* logs
  * ...
* msstats
  * [ComparisonPlot.pdf](#msstats-plots)
  * [VolcanoPlot.pdf](#msstats-plots)
  * [Heatmap.pdf](#msstats-plots)
  * [msstats\_results.csv](#msstats-table)
* pipeline\_info
  * [...](#nextflow-pipeline-info)
* proteomics\_lfq
  * [debug\_\*.idXML](#debug-output)
  * [out.consensusXML](#consenusxml)
  * [out.csv](#msstats-ready-quantity-table)
  * [out.mzTab](#mztab)
* ptxqc
  * [report\_v1.0.2\_out.yaml](#ptxqc-yaml-config)
  * [report\_v1.0.2\_out\_${hash}.html](#ptxqc-report)
  * [report\_v1.0.2\_out\_${hash}.pdf](#ptxqc-report)

### Nextflow pipeline info

### ProteomicsLFQ main output

#### ConsensusXML

#### MSstats-redy quantity table

#### mzTab

### MSstats output

#### MSstats table

#### MSstats plots

### PTXQC output

#### PTXQC report

#### PTXQC yaml config
