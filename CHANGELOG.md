# nf-core/proteomicslfq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0 - Lovely Logan [18.10.2020]

Initial release of nf-core/proteomicslfq, created with the [nf-core](https://nf-co.re/) template.

### `Added`

The initial version of the pipeline features the following steps:

    - (optional) Conversion of spectra data to indexedMzML: Using ThermoRawFileParser if Thermo Raw or using OpenMS' FileConverter if just an index is missing
    - (optional) Decoy database generation for the provided DB (fasta) with OpenMS
    - Database search with either MSGF+ and/or Comet through OpenMS adapters
    - Re-mapping potentially identified peptides to the input database for consistency and error-checking (using OpenMS' PeptideIndexer)
    - PSM rescoring using PSMFeatureExtractor and Percolator or a PeptideProphet-like distribution fitting approach in OpenMS
    - If multiple search engines were chosen, the results are combined with OpenMS' ConsensusID
    - If multiple search engines were chosen, a combined FDR is calculated
    - Single run PSM/Peptide-level FDR filtering
    - If localization of modifications was requested, Luciphor2 is applied via the OpenMS adapter
    - Protein inference and labelfree quantification based on spectral counting or MS1 feature detection, alignment and integration with OpenMS' ProteomicsLFQ. Performs an additional experiment-wide FDR filter on protein level (and if requested peptide/PSM level).

### `Known issues`

If you experience nextflow running forever after a failed step, try settings errorStrategy = terminate. See the corresponding [nextflow issue](https://github.com/nextflow-io/nextflow/issues/1457).

### `Fixed`

### `Dependencies`

### `Deprecated`
