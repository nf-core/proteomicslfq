# nf-core/proteomicslfq: Usage

## Table of contents

* [Table of contents](#table-of-contents)
* [Introduction](#introduction)
* [Running the pipeline](#running-the-pipeline)
  * [Updating the pipeline](#updating-the-pipeline)
  * [Reproducibility](#reproducibility)
* [Main arguments](#main-arguments)

  Either (using a PRIDE Sample to data relation format file):
  * [`--sdrf`](#--sdrf)
  * [`--root_folder`](#--root_folder)

  Or (using spectrum files and an OpenMS style experimental design):
  * [`--spectra`](#--spectra)
  * [`--exp_design`](#--exp_design)

  And:
  * [`--database`](#--database)
  * [`-profile`](#-profile)
* [Decoy database generation](#decoy-database-generation)
  * [`--add_decoys`](#--add_decoys)
  * [`--decoy_affix`](#-decoy_affix)
  * [`--affix_type`](#-profile)
* [Database search](#database-search)
  * [`--search_engine`](#--search_engine)
  * [`--enzyme`](#--enzyme)
  * [`--num_enzyme_termini`](#--num_enzyme_termini)
  * [`--num_hits`](#--num_hits)
  * [`--fixed_mods`](#--fixed_mods)
  * [`--variable_mods`](#--variable_mods)
  * [`--precursor_mass_tolerance`](#--precursor_mass_tolerance)
  * [`--precursor_mass_tolerance_unit`](#--precursor_mass_tolerance_unit)
  * [`--fragment_mass_tolerance`](#--fragment_mass_tolerance)
  * [`--fragment_mass_tolerance_unit`](#--fragment_mass_tolerance_unit)
  * [`--allowed_missed_cleavages`](#--allowed_missed_cleavages)
  * [`--psm_level_fdr_cutoff`](#--psm_level_fdr_cutoff)
  * [`--min_precursor_charge`](#--min_precursor_charge)
  * [`--max_precursor_charge`](#--max_precursor_charge)
  * [`--min_peptide_length`](#--min_peptide_length)
  * [`--max_peptide_length`](#--max_peptide_length)
  * [`--instrument`](#--instrument)
  * [`--protocol`](#--protocol)
  * [`--fragment_method`](#--fragment_method)
  * [`--isotope_error_range`](#--isotope_error_range)
  * [`--max_mods`](#--max_mods)
  * [`--db_debug`](#--db_debug)
* [Peptide reindexing](#peptide-reindexing)
  * [`--IL_equivalent`](#--IL_equivalent)
  * [`--allow_unmatched`](#--allow_unmatched)
* [PSM rescoring](#psm-rescoring)
  * [`--posterior_probabilities`](#--posterior_probabilities)
  * [`--rescoring_debug`](#--rescoring_debug)
  * [`--psm_pep_fdr_cutoff`](#--psm_pep_fdr_cutoff)
  * [Percolator specific](#percolator-specific)
    * [`--train_FDR`](#--train_FDR)
    * [`--test_FDR`](#--test_FDR)
    * [`--percolator_fdr_level`](#--percolator_fdr_level)
    * [`--description_correct_features`](#--description_correct_features)
    * [`--generic-feature-set`](#--feature)
    * [`--subset-max-train`](#--subset-max-train)
    * [`--klammer`](#--klammer)
  * [Distribution specific](#distribution-specific)
    * [`--outlier_handling`](#--outlier_handling)
    * [`--top_hits_only`](#--top_hits_only)
* [Inference and Quantification](#inference-and-quantification)
  * [`--inf_quant_debug`](#--inf_quant_debug)
  * [Inference](#inference)
    * [`--protein_inference`](#--protein_inference)
    * [`--protein_level_fdr_cutoff`](#--protein_level_fdr_cutoff)
  * [Quantification](#quantification)
    * [`--transfer_ids`](#--transfer_ids)
    * [`--targeted_only`](#--targeted_only)
    * [`--mass_recalibration`](#--mass_recalibration)
    * [`--psm_pep_fdr_for_quant`](#--psm_pep_fdr_for_quant)
    * [`--protein_quantification`](#--protein_quantification)
* [Statistical post-processing](#statistical-post-processing)
  * [`--skip_post_msstats`](#--skip_post_msstats)
  * [`--ref_condition`](#--ref_condition)
  * [`--contrasts`](#--contrasts)
* [Quality control](#quality-control)
  * [`--ptxqc_report_layout`](#--ptxqc_report_layout)
* [Job resources](#job-resources)
  * [Automatic resubmission](#automatic-resubmission)
  * [Custom resource requests](#custom-resource-requests)
* [AWS Batch specific parameters](#aws-batch-specific-parameters)
  * [`--awsqueue`](#--awsqueue)
  * [`--awsregion`](#--awsregion)
  * [`--awscli`](#--awscli)
* [Other command line parameters](#other-command-line-parameters)
  * [`--outdir`](#--outdir)
  * [`--email`](#--email)
  * [`--email_on_fail`](#--email_on_fail)
  * [`-name`](#-name)
  * [`-resume`](#-resume)
  * [`-c`](#-c)
  * [`--custom_config_version`](#--custom_config_version)
  * [`--custom_config_base`](#--custom_config_base)
  * [`--max_memory`](#--max_memory)
  * [`--max_time`](#--max_time)
  * [`--max_cpus`](#--max_cpus)
  * [`--plaintext_email`](#--plaintext_email)
  * [`--monochrome_logs`](#--monochrome_logs)

## Introduction

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline

The most simple command for running the pipeline is as follows:

```bash
nextflow run nf-core/proteomicslfq --spectra '*.mzML' --database '*.fasta' -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/proteomicslfq
```

### Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/proteomicslfq releases page](https://github.com/nf-core/proteomicslfq/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For running a (not necessarily stable) development version of the pipeline you can use the `dev` branch with `-r dev`.

## Main arguments

The input to the pipeline can be specified in two mutually exclusive ways:

-----

__a)__ Either by using a path or URI to a PRIDE Sample to data relation format file (SDRF), e.g. as part of a submitted and
annotated PRIDE experiment (see here for examples). For this case, use:

### `--sdrf`

For the URI or path to the SDRF file. Input files will be downloaded and cached from the URIs specified in the SDRF file.
An OpenMS-style experimental design will be generated based on the factor columns of the SDRF. The settings for the
following parameters will currently be overwritten by the ones specified in the SDRF:

* `fixed_mods`,
* `variable_mods`,
* `precursor_mass_tolerance`,
* `precursor_mass_tolerance_unit`,
* `fragment_mass_tolerance`,
* `fragment_mass_tolerance_unit`,
* `fragment_method`,
* `enzyme`

### `--root_folder`

This optional parameter can be used to specify a root folder in which the spectrum files specified in the SDRF are searched.
It is usually used if you have a local version of the experiment already. Note that this option does not support recursive
searching yet.

-----

__b)__ By specifying globbing patterns to the input spectrum files in Thermo RAW or mzML format and a manual OpenMS-style
experimental design file.

### `--spectra`

Use this to specify the location of your input mzML or Thermo RAW files:

```bash
--spectra 'path/to/data/*.mzML'
```

or

```bash
--spectra 'path/to/data/*.raw'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character

-----

### `--exp_design`

Path or URL to an experimental design file (if not given, it assumes unfractionated, unrelated samples). See an example here (TODO).

```bash
--exp_design '[path to experimental design file in OpenMS-style tab separated format]'
```

### `--database`

Since the database is not included in an SDRF, this parameter always needs to be given to specify the input protein database
when you run the pipeline. Remember to include contaminants (and decoys if not added in the pipeline with --add-decoys)

```bash
--database '[path to Fasta protein database]'
```

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. Note that multiple profiles can be loaded, for example: `-profile docker` - the order of arguments is important!

If `-profile` is not specified at all the pipeline will be run locally and expects all software to be installed and available on the `PATH`.

* `docker`
  * A generic configuration profile to be used with [Docker](http://docker.com/)
  * Pulls software from Docker Hub: [`nfcore/proteomicslfq`](http://hub.docker.com/r/nfcore/proteomicslfq/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
  * Pulls software from Docker Hub: [`nfcore/proteomicslfq`](http://hub.docker.com/r/nfcore/proteomicslfq/)
* `conda`
  * Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker or Singularity.
  * A generic configuration profile to be used with [Conda](https://conda.io/docs/)
  * Pulls most software from the [Bioconda](https://bioconda.github.io/) and [conda-forge](https://conda-forge.org/) channels.
* `test`
  * A profile with a complete configuration for automated testing
  * Includes links to test data hosted on GitHub and therefore doesn't need additional parameters
* `test_full`
  * A profile with a complete configuration for automated testing on AWS
  * Includes links to test data on GitHub and PRIDE and therefore doesn't need additional parameters
  * Warning: Downloads roughly 9GB of raw data from PRIDE and analyzes

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

## Decoy database generation

### `--add_decoys`

If decoys were not yet included in the input database, they have to be appended by OpenMS DecoyGenerator by adding this flag
(TODO allow specifying type).
Default: pseudo-reverse peptides

### `--decoy_affix`

Specify the string that was or will be added to the protein accession to label it

### `--affix_type`

Is the decoy label a prefix or suffix. Prefix is highly recommended as some tools (e.g. Percolator) might not work well with suffixes

## Database search

### `--precursor_mass_tolerance`

Precursor mass tolerance used for database search. For High-Resolution instruments a precursor mass tolerance value of 5 ppm is recommended (i.e. 5). See also [`--precursor_mass_tolerance_unit`](#--precursor_mass_tolerance_unit).

### `--precursor_mass_tolerance_unit`

Precursor mass tolerance unit used for database search. Possible values are "ppm" (default) and "Da".

### `--enzyme`

Specify which enzymatic restriction should be applied, e.g. 'unspecific cleavage', 'Trypsin' (default), see OpenMS
[enzymes](https://github.com/OpenMS/OpenMS/blob/develop/share/OpenMS/CHEMISTRY/Enzymes.xml). Note: MSGF does not support extended
cutting rules, as used by default with "Trypsin". I.e. if you specify "Trypsin" with MSGF, it will be automatically converted to
"Trypsin/P" = "Trypsin without proline rule".

### `--num_enzyme_termini`

Specify the number of termini matching the enzyme cutting rules for a peptide to be considered. Valid values are
"fully" (default), "semi", or "none".

### `--num_hits`

Specify the maximum number of top peptide candidates per spectrum to be reported by the search engine. Default: 1

### `--fixed_mods`

Specify which fixed modifications should be applied to the database search (eg. '' or 'Carbamidomethyl (C)', see Unimod modifications
in the style '({unimod name} ({optional term specificity} {optional origin})').
All possible enzymes can be found in the restrictions mentioned in the command line documentation of e.g. [CometAdapter](https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_CometAdapter.html) (scroll down a bit for the complete set).
Multiple fixed modifications can be specified comma separated (e.g. 'Carbamidomethyl (C),Oxidation (M)')

### `--variable_mods`

Specify which variable modifications should be applied to the database search (eg. 'Oxidation (M)', see Unimod modifications
in the style '({unimod name} ({optional term specificity} {optional origin})').
All possible enzymes can be found in the restrictions mentioned in the command line documentation of e.g. [CometAdapter](https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/TOPP_CometAdapter.html) (scroll down a bit for the complete set).
Multiple variable modifications can be specified comma separated (e.g. 'Carbamidomethyl (C),Oxidation (M)')

### `--allowed_missed_cleavages`

Specify the maximum number of allowed missed enzyme cleavages in a peptide. The parameter is not applied if "unspecific cleavage" is specified as enzyme.

### `--instrument`

Type of instrument that generated the data. 'low_res' or 'high_res' (default; refers to LCQ and LTQ instruments)

### `--protocol`

MSGF only: Labeling or enrichment protocol used, if any

### `--fragment_method`

Currently unsupported. Defaults to "ALL" for Comet and "from_spectrum" for MSGF. Should be a sensible default for 99% of the cases.

### `--isotope_error_range`

Range of allowed isotope peak errors (MS-GF+ parameter '-ti'). Takes into account the error introduced by choosing a non-monoisotopic peak for fragmentation. Combined with 'precursor_mass_tolerance'/'precursor_error_units', this determines the actual precursor mass tolerance. E.g. for experimental mass 'exp' and calculated mass 'calc', '-precursor_mass_tolerance 20 -precursor_error_units ppm -isotope_error_range -1,2' tests '|exp - calc - n * 1.00335 Da| < 20 ppm' for n = -1, 0, 1, 2.

### `--min_precursor_charge`

Minimum precursor ion charge

### `--max_precursor_charge`

Maximum precursor ion charge

### `--min_peptide_length`

Minimum peptide length to consider (works with MSGF and in newer Comet versions)

### `--max_peptide_length`

Maximum peptide length to consider (works with MSGF and in newer Comet versions)

### `--max_mods`

Maximum number of modifications per peptide. If this value is large, the search may take very long.

### `--db_debug`

Set debug level for the search engines (regulates if intermediate output is kept and if you are going to see the output
of the underlying search engine)

## Peptide reindexing

### `--IL_equivalent`

Should isoleucine and leucine be treated interchangeably? Default: true

### `--allow_unmatched`

Ignore unmatched peptides (Default: false; only activate if you double-checked all other settings)

## PSM Rescoring

### `--posterior_probabilities`

How to calculate posterior probabilities for PSMs:

* "percolator" = Re-score based on PSM-feature-based SVM and transform distance
    to hyperplane for posteriors
* "fit_distributions" = Fit positive and negative distributions to scores
    (similar to PeptideProphet)

### `--rescoring_debug`

Debug level during PSM rescoring for additional text output and keeping temporary files

### `--psm_pep_fdr_cutoff`

FDR cutoff on PSM level (or potential peptide level; see Percolator options)
before going into feature finding, map alignment and inference.

### Percolator specific

In the following you can find help for the Percolator specific options that are only used if [`--posterior_probabilities`](#--posterior_probabilities) was set to "percolator".
Note that there are currently some restrictions to the original options of Percolator:

* no Percolator protein FDR possible (currently OpenMS' FDR is used on protein level)
* no support for separate target and decoy databases (i.e. no min-max q-value calculation or target-decoy competition strategy)
* no support for combined or experiment-wide peptide re-scoring. Currently search results per input file are submitted to Percolator independently.

With time, some of the limitations might be removed. Pull requests are always welcome.

### `--train_FDR`

False discovery rate threshold to define positive examples in training. Set to testFDR if 0

### `--test_FDR`

False discovery rate threshold for evaluating best cross validation result and reported end result

### `--percolator_fdr_level`

Level of FDR calculation ('peptide-level-fdrs' or 'psm-level-fdrs')

### `--description_correct_features`

Percolator provides the possibility to use so called description of correct features, i.e. features for which desirable values are learnt from the previously identified target PSMs. The absolute value of the difference between desired value and observed value is the used as predictive features.

1 -> iso-electric point

2 -> mass calibration

4 -> retention time

8 -> delta_retention_time\*delta_mass_calibration

### `--generic-feature-set`

Use only generic (i.e. not search engine specific) features. Generating search engine specific
features for common search engines by PSMFeatureExtractor will typically boost the identification rate significantly.

### `--subset-max-train`

Only train an SVM on a subset of PSMs, and use the resulting score vector to evaluate the other
PSMs. Recommended when analyzing huge numbers (>1 million) of PSMs. When set to 0, all PSMs are used for training as normal.
Default: 300,000

### `--klammer`

Retention time features are calculated as in Klammer et al. instead of with Elude. Default: false

### Distribution-fitting (IDPEP) specific

Use this instead of Percolator if there are problems with Percolator (e.g. due to bad separation) or for performance
reasons.

### `--outlier_handling`

How to handle outliers during fitting:

* ignore_iqr_outliers (default): ignore outliers outside of 3\*IQR from Q1/Q3 for fitting
* set_iqr_to_closest_valid: set IQR-based outliers to the last valid value for fitting
* ignore_extreme_percentiles: ignore everything outside 99th and 1st percentile (also removes equal values like potential censored max values in XTandem)
* none: do nothing

### `--top_hits_only`

Use only the top peptide hits per spectrum for fitting. Default: true

## Inference and Quantification

### `--inf_quant_debug`

Debug level during inference and quantification. (WARNING: Higher than 666 may produce a lot of additional output files)

### Inference

### `--protein_inference`

Infer proteins through:

* "aggregation"  = aggregates all peptide scores across a protein (by calculating the maximum) (default)
* "bayesian"     = computes a posterior probability for every protein based on a Bayesian network
* ("percolator" not yet supported)

### `--protein_level_fdr_cutoff`

Protein level FDR cutoff (Note: this affects and chooses the peptides used for quantification). Default: 0.05

### Quantification

### `--transfer_ids`

Transfer IDs over aligned samples to increase # of quantifiable features (WARNING: increased memory consumption). (default: "false")

### `--targeted_only`

Only looks for quantifiable features at locations with an identified spectrum. (default: "true")

### `--mass_recalibration`

Recalibrates masses to correct for instrument biases. (default: "false")

### `--psm_pep_fdr_for_quant`

PSM/peptide level FDR used for quantification after inference (*in addition to protein-level filtering*)
If Bayesian inference was chosen, this will be a peptide-level FDR and only the best PSMs per
peptide will be reported. (default: off = 1.0)

### `--protein_quantification`

Quantify proteins based on:

* "unique_peptides" = use peptides mapping to single proteins or a group of indistinguishable proteins (according to the set of experimentally identified peptides)
* "strictly_unique_peptides" = use peptides mapping to a unique single protein only
* "shared_peptides" = use shared peptides, too, but only greedily for its best group (by inference score)

## Statistical post-processing

### `--skip_post_msstats`

Skip MSstats for statistical post-processing?

### `--ref_condition`

Instead of all pairwise contrasts (default), uses the given condition name/number (corresponding to your experimental design) as a reference and creates pairwise contrasts against it. (TODO not yet fully implemented)

### `--contrasts`

Specify a set of contrasts in a semicolon seperated list of R-compatible contrasts with the
condition names/numbers as variables (e.g. "1-2;1-3;2-3"). Overwrites "--reference" (TODO not yet fully implemented)

## Quality control

### `--skip_qc`

Skip generation of quality control report by PTXQC? default: "true" since it is still unstable

### `--ptxqc_report_layout`

Specify a yaml file for the report layout (see PTXQC documentation) (TODO not yet fully implemented)

Note that you can use the same configuration setup to save sets of reference files for your own use, even if they are not part of the iGenomes resource. See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for instructions on where to save such a file.

## Job resources

### Automatic resubmission

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests

Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files hosted at [`nf-core/configs`](https://github.com/nf-core/configs/tree/master/conf) for examples.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack).

## AWS Batch specific parameters

Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use [`-profile awsbatch`](https://github.com/nf-core/configs/blob/master/conf/awsbatch.config) and then specify all of the following parameters.

### `--awsqueue`

The JobQueue that you intend to use on AWS Batch.

### `--awsregion`

The AWS region in which to run your job. Default is set to `eu-west-1` but can be adjusted to your needs.

### `--awscli`

The [AWS CLI](https://www.nextflow.io/docs/latest/awscloud.html#aws-cli-installation) path in your custom AMI. Default: `/home/ec2-user/miniconda/bin/aws`.

Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.

## Other command line parameters

### `--outdir`

The output directory where the results will be saved.

### `--email`

Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

### `--email_on_fail`

This works exactly as with `--email`, except emails are only sent if the workflow is not successful.

### `-name`

Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`

Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### `--custom_config_version`

Provide git commit id for custom Institutional configs hosted at `nf-core/configs`. This was implemented for reproducibility purposes. Default: `master`.

```bash
## Download and use config file with following git commid id
--custom_config_version d52db660777c4bf36546ddb188ec530c3ada1b96
```

### `--custom_config_base`

If you're running offline, nextflow will not be able to fetch the institutional config files
from the internet. If you don't need them, then this is not a problem. If you do need them,
you should download the files from the repo and tell nextflow where to find them with the
`custom_config_base` option. For example:

```bash
## Download and unzip the config files
cd /path/to/my/configs
wget https://github.com/nf-core/configs/archive/master.zip
unzip master.zip

## Run the pipeline
cd /path/to/my/data
nextflow run /path/to/pipeline/ --custom_config_base /path/to/my/configs/configs-master/
```

> Note that the nf-core/tools helper package has a `download` command to download all required pipeline
> files + singularity containers + institutional configs in one go for you, to make this process easier.

### `--max_memory`

Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`

### `--max_time`

Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`

Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--plaintext_email`

Set to receive plain-text e-mails instead of HTML formatted.

### `--monochrome_logs`

Set to disable colourful command line output and live life in monochrome.
