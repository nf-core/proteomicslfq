#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/proteomicslfq
========================================================================================
 nf-core/proteomicslfq Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/proteomicslfq
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/proteomicslfq --spectra '*.mzML' --database '*.fasta' -profile docker

    Main arguments:
      Either:
      --sdrf                        Path to PRIDE Sample to data relation format file
      --root_folder                 (Optional) If given, looks for the filenames in the SDRF in this folder, locally
      Or:
      --spectra                     Path to input spectra as mzML or Thermo Raw
      --expdesign                   Path to optional experimental design file (if not given, it assumes unfractionated, unrelated samples)

      And:
      --database                    Path to input protein database as fasta

    Decoy database:
      --add_decoys                  Add decoys to the given fasta
      --decoy_affix                 The decoy prefix or suffix used or to be used (default: DECOY_)
      --affix_type                  Prefix (default) or suffix (WARNING: Percolator only supports prefices)

    Database Search:
      --search_engine               Which search engine: "comet" (default) or "msgf"
      --enzyme                      Enzymatic cleavage (e.g. 'unspecific cleavage' or 'Trypsin' [default], see OpenMS enzymes)
      --num_enzyme_termini          Specify the termini where the cleavage rule has to match (default:
                                         'fully' valid: 'semi', 'fully')
      --num_hits                    Number of peptide hits per spectrum (PSMs) in output file (default: '1')
      --fixed_mods                  Fixed modifications ('Carbamidomethyl (C)', see OpenMS modifications)
      --variable_mods               Variable modifications ('Oxidation (M)', see OpenMS modifications)
      --enable_mod_localization     Enable localization scoring with Luciphor
      --mod_localization            Specify the var. modifications whose localizations should be rescored with the luciphor algorithm
      --precursor_mass_tolerance    Mass tolerance of precursor mass
      --precursor_mass_tolerance_unit Da or ppm
      --fragment_mass_tolerance     Mass tolerance for fragment masses (currently only controls Comets fragment_bin_tol)
      --fragment_mass_tolerance_unit Da or ppm (currently always ppm)
      --allowed_missed_cleavages    Allowed missed cleavages
      --min_precursor_charge        Minimum precursor ion charge
      --max_precursor_charge        Maximum precursor ion charge
      --min_peptide_length          Minimum peptide length to consider
      --max_peptide_length          Maximum peptide length to consider
      --instrument                  Type of instrument that generated the data (currently only 'high_res' [default] and 'low_res' supported)
      --protocol                    Used labeling or enrichment protocol (if any)
      --fragment_method             Used fragmentation method (currently unused since we let the search engines consider all MS2 spectra and let                                     them determine from the spectrum metadata)
      --max_mods                    Maximum number of modifications per peptide. If this value is large, the search may take very long
      --db_debug                    Debug level during database search

      //TODO probably also still some options missing. Try to consolidate them whenever the two search engines share them

    PSM Rescoring:
      --posterior_probabilities     How to calculate posterior probabilities for PSMs:
                                    "percolator" = Re-score based on PSM-feature-based SVM and transform distance
                                        to hyperplane for posteriors
                                    "fit_distributions" = Fit positive and negative distributions to scores
                                        (similar to PeptideProphet)
      --rescoring_debug             Debug level during PSM rescoring
      --psm_pep_fdr_cutoff          FDR cutoff on PSM level (or potential peptide level; see Percolator options) before going into
                                    feature finding, map alignment and inference.

      Percolator specific:
      --train_FDR                   False discovery rate threshold to define positive examples in training. Set to testFDR if 0
      --test_FDR                    False discovery rate threshold for evaluating best cross validation result and reported end result
      --percolator_fdr_level        Level of FDR calculation ('peptide-level-fdrs' or 'psm-level-fdrs')
      --description_correct_features Description of correct features for Percolator (0, 1, 2, 4, 8, see Percolator retention time and calibration)
      --generic-feature-set         Use only generic (i.e. not search engine specific) features. Generating search engine specific
                                    features for common search engines by PSMFeatureExtractor will typically boost the identification rate significantly.
      --subset-max-train            Only train an SVM on a subset of PSMs, and use the resulting score vector to evaluate the other
                                    PSMs. Recommended when analyzing huge numbers (>1 million) of PSMs. When set to 0, all PSMs are used for training as normal.
      --klammer                     Retention time features are calculated as in Klammer et al. instead of with Elude

      Distribution specific:
      --outlier_handling            How to handle outliers during fitting:
                                    - ignore_iqr_outliers (default): ignore outliers outside of 3*IQR from Q1/Q3 for fitting
                                    - set_iqr_to_closest_valid: set IQR-based outliers to the last valid value for fitting
                                    - ignore_extreme_percentiles: ignore everything outside 99th and 1st percentile (also removes equal values like potential censored max values in XTandem)
                                    - none: do nothing
      --top_hits_only               Use only the top hits for fitting

      //TODO add more options for rescoring part

    Inference and Quantification:
      --inf_quant_debug             Debug level during inference and quantification. (WARNING: Higher than 666 may produce a lot
                                    of additional output files)
      Inference:
      --protein_inference           Infer proteins through:
                                    "aggregation"  = aggregates all peptide scores across a protein (by calculating the maximum)
                                    "bayesian"     = computes a posterior probability for every protein based on a Bayesian network
                                    ("percolator" not yet supported)
      --protein_level_fdr_cutoff    Protein level FDR cutoff (this affects and chooses the peptides used for quantification)

      Quantification:
      --transfer_ids                Transfer IDs over aligned samples to increase # of quantifiable features (WARNING:
                                    increased memory consumption). (default: false) TODO must specify true or false
      --targeted_only               Only ID based quantification. (default: true) TODO must specify true or false
      --mass_recalibration          Recalibrates masses to correct for instrument biases. (default: false) TODO must specify true
                                    or false

      //TODO the following need to be passed still
      --psm_pep_fdr_for_quant       PSM/peptide level FDR used for quantification (if filtering on protein level is not enough)
                                    If Bayesian inference was chosen, this will be a peptide-level FDR and only the best PSMs per
                                    peptide will be reported.
                                    (default: off = 1.0)
      --protein_quantification      Quantify proteins based on:
                                    "unique_peptides" = use peptides mapping to single proteins or a group of indistinguishable proteins (according to the set of experimentally identified peptides)
                                    "strictly_unique_peptides" = use peptides mapping to a unique single protein only
                                    "shared_peptides" = use shared peptides only for its best group (by inference score)

    Statistical post-processing:
      --skip_post_msstats           Skip MSstats for statistical post-processing?
      --ref_condition               Instead of all pairwise contrasts, uses the given condition number (corresponding to your experimental design) as a reference and
                                    creates pairwise contrasts against it (TODO fully implement)
      --contrasts                   Specify a set of contrasts in a semicolon seperated list of R-compatible contrasts with the
                                    condition numbers as variables (e.g. "1-2;1-3;2-3"). Overwrites "--reference" (TODO fully implement)

    Quality control:
      --ptxqc_report_layout         Specify a yaml file for the report layout (see PTXQC documentation) (TODO fully implement)

    Other nextflow options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the
                                    run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random
                                    mnemonic.
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity, awsbatch, test and more.

    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Stage config files
ch_output_docs = Channel.fromPath("$baseDir/docs/output.md")


// Validate inputs
if (!(params.spectra || params.sdrf) || (params.spectra && params.sdrf))
{
  log.error "EITHER spectra data OR SDRF needs to be provided. Make sure you have used either of those options."; exit 1
}

params.database = params.database ?: { log.error "No protein database provided. Make sure you have used the '--database' option."; exit 1 }()
params.outdir = params.outdir ?: { log.warn "No output directory provided. Will put the results into './results'"; return "./results" }()

/*
 * Create a channel for input files
 */

 //Filename        FixedModifications      VariableModifications   Label   PrecursorMassTolerance  PrecursorMassToleranceUnit      FragmentMassTolerance   DissociationMethod      Enzyme


if (!params.sdrf)
{
  ch_spectra = Channel.fromPath(params.spectra, checkIfExists: true)
  ch_spectra
  .multiMap{ it -> id = UUID.randomUUID().toString()
                    comet_settings: msgf_settings: tuple(id,
                                    params.fixed_mods,
                                    params.variable_mods,
                                    "", //labelling modifications currently not supported
                                    params.precursor_mass_tolerance,
                                    params.precursor_mass_tolerance_unit,
                                    params.fragment_mass_tolerance,
                                    params.fragment_mass_tolerance_unit,
                                    params.fragment_method,
                                    params.enzyme)
                    idx_settings: tuple(id,
                                    params.enzyme)
                    luciphor_settings: 
                                  tuple(id,
                                    params.fragment_method)
                    mzmls: tuple(id,it)}
  .set{ch_sdrf_config}
}
else
{
  ch_sdrf = Channel.fromPath(params.sdrf, checkIfExists: true)
  /*
   * STEP 0 - SDRF parsing
   */
  process sdrf_parsing {

      publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

      input:
       file sdrf from ch_sdrf

      output:
       file "experimental_design.tsv" into ch_expdesign
       file "openms.tsv" into ch_sdrf_config_file

      when:
        params.sdrf

      script:
       """
       ## -t2 since the one-table format parser is broken in OpenMS2.5
       ## -l for legacy behavior to always add sample columns
       parse_sdrf.py convert-openms -t2 -l -s ${sdrf} > sdrf_parsing.log
       """
  }

  //TODO use header and reference by col name instead of index
  ch_sdrf_config_file
  .splitCsv(skip: 1, sep: '\t')
  .multiMap{ row -> id = UUID.randomUUID().toString()
                    comet_settings: msgf_settings: tuple(id,
                                    row[2],
                                    row[3],
                                    row[4],
                                    row[5],
                                    row[6],
                                    row[7],
                                    row[8],
                                    row[9],
                                    row[10])
                    idx_settings: tuple(id,
                                    row[10])
                    luciphor_settings: 
                                  tuple(id,
                                    row[9])
                    mzmls: tuple(id, params.root_folder.length() == 0 ? row[0] : (params.root_folder + "/" + row[1]))}
  .set{ch_sdrf_config}
}

ch_db_for_decoy_creation = Channel.fromPath(params.database)

// overwrite experimental design if given additionally to SDRF
//TODO think about that
if (params.expdesign)
{
    Channel
        .fromPath(params.expdesign)
        .set { ch_expdesign }
}

ch_sdrf_config.mzmls
.branch {
        raw: hasExtension(it[1], 'raw')
        mzML: hasExtension(it[1], 'mzML')
}
.set {branched_input}


//TODO we could also check for outdated mzML versions and try to update them
branched_input.mzML
.branch {
    nonIndexedMzML: file(it[1]).withReader {
                        f = it;
                        1.upto(5) {
                            if (f.readLine().contains("indexedmzML")) return false;
                        }
                        return true;
                    }
    inputIndexedMzML: file(it[1]).withReader {
                        f = it;
                        1.upto(5) {
                            if (f.readLine().contains("indexedmzML")) return true;
                        }
                        return false;
                    }
}
.set {branched_input_mzMLs}

//Push raw files through process that does the conversion, everything else directly to downstream Channel with mzMLs


//This piece only runs on data that is a.) raw and b.) needs conversion
//mzML files will be mixed after this step to provide output for downstream processing - allowing you to even specify mzMLs and RAW files in a mixed mode as input :-)


//GENERAL TODOS
// - Check why we depend on full filepaths and if that is needed
/* Proposition from nextflow gitter https://gitter.im/nextflow-io/nextflow?at=5e25fabea259cb0f0607a1a1
*
* unless the specific filenames are important (depends on the tool you're using), I usually use the pattern outlined here:
* https://www.nextflow.io/docs/latest/process.html#multiple-input-files
* e.g: file "?????.mzML" from mzmls_plfq.toSortedList() and ProteomicsLFQ -in *.mzML -ids *.id
*/

/*
 * STEP 0.1 - Raw file conversion
 */
process raw_file_conversion {

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
     tuple mzml_id, path(rawfile) from branched_input.raw

    output:
     set mzml_id, file("*.mzML") into mzmls_converted

    script:
     """
     ThermoRawFileParser.sh -i=${rawfile} -f=2 -o=./ > ${rawfile}_conversion.log
     """
}

/*
 * STEP 0.2 - MzML indexing
 */
process mzml_indexing {

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
     set mzml_id, path(mzmlfile) from branched_input_mzMLs.nonIndexedMzML

    output:
     set mzml_id, file("out/*.mzML") into mzmls_indexed
     file "*.log"

    script:
     """
     mkdir out
     FileConverter -in ${mzmlfile} -out out/${mzmlfile.baseName}.mzML > ${mzmlfile.baseName}_mzmlindexing.log
     """
}

//Mix the converted raw data with the already supplied mzMLs and push these to the same channels as before

branched_input_mzMLs.inputIndexedMzML.mix(mzmls_converted).mix(mzmls_indexed).into{mzmls_comet; mzmls_msgf; mzmls_luciphor; mzmls_plfq}


//Fill the channels with empty Channels in case that we want to add decoys. Otherwise fill with output from database.
(searchengine_in_db_msgf, searchengine_in_db_comet, pepidx_in_db, plfq_in_db) = ( params.add_decoys
                    ? [ Channel.empty(), Channel.empty(), Channel.empty(), Channel.empty() ]
                    : [ Channel.fromPath(params.database), Channel.fromPath(params.database), Channel.fromPath(params.database), Channel.fromPath(params.database)  ] )

//Add decoys if params.add_decoys is set appropriately
process generate_decoy_database {

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
     file(mydatabase) from ch_db_for_decoy_creation

    output:
     file "${mydatabase.baseName}_decoy.fasta" into searchengine_in_db_decoy_msgf, searchengine_in_db_decoy_comet, pepidx_in_db_decoy, plfq_in_db_decoy
     file "*.log"

    when:
     params.add_decoys

    script:
     """
     DecoyDatabase  -in ${mydatabase} \\
                 -out ${mydatabase.baseName}_decoy.fasta \\
                 -decoy_string ${params.decoy_affix} \\
                 -decoy_string_position ${params.affix_type} \\
                 > ${mydatabase.baseName}_decoy_database.log
     """
}

// Doesnt work. Py script needs all the inputs to be together in a folder
// Wont work with nextflow. It needs to accept a list of paths for the inputs!!
//process generate_simple_exp_design_file {
//    publishDir "${params.outdir}", mode: 'copy'
//    input:
//      val mymzmls from mzmls.collect()

//    output:
//        file "expdesign.tsv" into expdesign
//    when:
//        !params.expdesign

//    script:
//     strng = new File(mymzmls[0].toString()).getParentFile()
//     """
//       create_trivial_design.py ${strng} 1 > expdesign.tsv
//     """
//}

if (params.enzyme == "unspecific cleavage")
{
  params.num_enzyme_termini == "none"
}

pepidx_num_enzyme_termini = params.num_enzyme_termini
if (params.num_enzyme_termini == "fully")
{
  pepidx_num_enzyme_termini = "full"
}

/// Search engine
if (params.search_engine == "msgf")
{
    search_engine_score = "SpecEValue"
} else { //comet
    search_engine_score = "expect"
}

process search_engine_msgf {

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    // ---------------------------------------------------------------------------------------------------------------------
    // ------------- WARNING: THIS IS A HACK. IT JUST DOES NOT WORK IF THIS PROCESS IS RETRIED -----------------------------
    // ---------------------------------------------------------------------------------------------------------------------
    // I actually dont know, where else this would be needed.
    errorStrategy 'terminate'

    input:
     tuple file(database), mzml_id, path(mzml_file), fixed, variable, label, prec_tol, prec_tol_unit, frag_tol, frag_tol_unit, diss_meth, enzyme from searchengine_in_db_msgf.mix(searchengine_in_db_decoy_msgf).combine(mzmls_msgf.join(ch_sdrf_config.msgf_settings))

     // This was another way of handling the combination
     //file database from searchengine_in_db.mix(searchengine_in_db_decoy)
     //each file(mzml_file) from mzmls
    when:
      params.search_engine == "msgf"

    output:
     set mzml_id, file("${mzml_file.baseName}.idXML") into id_files_msgf
     file "*.log"

    script:
      if (enzyme == 'Trypsin') enzyme = 'Trypsin/P'
      else if (enzyme == 'Arg-C') enzyme = 'Arg-C/P'
      else if (enzyme == 'Asp-N') enzyme = 'Asp-N/B'
      else if (enzyme == 'Chymotrypsin') enzyme = 'Chymotrypsin/P'
      else if (enzyme == 'Lys-C') enzyme = 'Lys-C/P'

     """
     MSGFPlusAdapter -in ${mzml_file} \\
                     -out ${mzml_file.baseName}.idXML \\
                     -threads ${task.cpus} \\
                     -database "${database}" \\
                     -instrument ${params.instrument} \\
                     -protocol "${params.protocol}" \\
                     -matches_per_spec ${params.num_hits} \\
                     -min_precursor_charge ${params.min_precursor_charge} \\
                     -max_precursor_charge ${params.max_precursor_charge} \\
                     -min_peptide_length ${params.min_peptide_length} \\
                     -max_peptide_length ${params.max_peptide_length} \\
                     -enzyme "${enzyme}" \\
                     -tryptic ${params.num_enzyme_termini} \\
                     -precursor_mass_tolerance ${prec_tol} \\
                     -precursor_error_units ${prec_tol_unit} \\
                     -fixed_modifications ${fixed.tokenize(',').collect { "'${it}'" }.join(" ") } \\
                     -variable_modifications ${variable.tokenize(',').collect { "'${it}'" }.join(" ") } \\
                     -max_mods ${params.max_mods} \\
                     -debug ${params.db_debug} \\
                     > ${mzml_file.baseName}_msgf.log
     """
}

process search_engine_comet {

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    // ---------------------------------------------------------------------------------------------------------------------
    // ------------- WARNING: THIS IS A HACK. IT JUST DOES NOT WORK IF THIS PROCESS IS RETRIED -----------------------------
    // ---------------------------------------------------------------------------------------------------------------------
    // I actually dont know, where else this would be needed.
    errorStrategy 'terminate'
    input:
     tuple file(database), mzml_id, path(mzml_file), fixed, variable, label, prec_tol, prec_tol_unit, frag_tol, frag_tol_unit, diss_meth, enzyme from searchengine_in_db_comet.mix(searchengine_in_db_decoy_comet).combine(mzmls_comet.join(ch_sdrf_config.comet_settings))

     //or
     //file database from searchengine_in_db_comet.mix(searchengine_in_db_decoy_comet)
     //each file(mzml_file) from mzmls

    when:
      params.search_engine == "comet"

    output:
     set mzml_id, file("${mzml_file.baseName}.idXML") into id_files_comet
     file "*.log"

    //TODO we currently ignore the activation_method param to leave the default "ALL" for max. compatibility
    script:
     """
     CometAdapter  -in ${mzml_file} \\
                   -out ${mzml_file.baseName}.idXML \\
                   -threads ${task.cpus} \\
                   -database "${database}" \\
                   -instrument ${params.instrument} \\
                   -allowed_missed_cleavages ${params.allowed_missed_cleavages} \\
                   -num_hits ${params.num_hits} \\
                   -num_enzyme_termini ${params.num_enzyme_termini} \\
                   -enzyme "${enzyme}" \\
                   -precursor_charge ${params.min_precursor_charge}:${params.max_precursor_charge} \\
                   -fixed_modifications ${fixed.tokenize(',').collect { "'${it}'" }.join(" ") } \\
                   -variable_modifications ${variable.tokenize(',').collect { "'${it}'" }.join(" ") } \\
                   -max_variable_mods_in_peptide ${params.max_mods} \\
                   -precursor_mass_tolerance ${prec_tol} \\
                   -precursor_error_units ${prec_tol_unit} \\
                   -fragment_bin_tolerance ${frag_tol} \\
                   -debug ${params.db_debug} \\
                   > ${mzml_file.baseName}_comet.log
     """
}


process index_peptides {

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
     tuple mzml_id, file(id_file), enzyme, file(database) from id_files_msgf.mix(id_files_comet).join(ch_sdrf_config.idx_settings).combine(pepidx_in_db.mix(pepidx_in_db_decoy))

     //each mzml_id, file(id_file) from id_files_msgf.mix(id_files_comet)
     //file database from pepidx_in_db.mix(pepidx_in_db_decoy)

    output:
     set mzml_id, file("${id_file.baseName}_idx.idXML") into id_files_idx_ForPerc, id_files_idx_ForIDPEP
     file "*.log"

    script:
     """
     PeptideIndexer -in ${id_file} \\
                    -out ${id_file.baseName}_idx.idXML \\
                    -threads ${task.cpus} \\
                    -fasta ${database} \\
                    -enzyme:name "${enzyme}" \\
                    -enzyme:specificity ${pepidx_num_enzyme_termini}
                    > ${id_file.baseName}_index_peptides.log
     """
}


// ---------------------------------------------------------------------
// Branch a) Q-values and PEP from Percolator


process extract_percolator_features {

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
     set mzml_id, file(id_file) from id_files_idx_ForPerc

    output:
     set mzml_id, file("${id_file.baseName}_feat.idXML") into id_files_idx_feat
     file "*.log"

    when:
     params.posterior_probabilities == "percolator"

    script:
     """
     PSMFeatureExtractor -in ${id_file} \\
                         -out ${id_file.baseName}_feat.idXML \\
                         -threads ${task.cpus} \\
                         > ${id_file.baseName}_extract_percolator_features.log
     """
}


//Note: from here, we do not need any settings anymore. so we can skip adding the mzml_id to the channels
//TODO parameterize and find a way to run across all runs merged
process percolator {

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
     set mzml_id, file(id_file) from id_files_idx_feat

    output:
     set mzml_id, file("${id_file.baseName}_perc.idXML") into id_files_idx_feat_perc
     file "*.log"

    when:
     params.posterior_probabilities == "percolator"

    // NICE-TO-HAVE: the decoy-pattern is automatically detected from PeptideIndexer.
    // Parse its output and put the correct one here.
    script:
        if (params.klammer && params.description_correct_features == 0) {
            log.warn('Klammer was specified, but description of correct features was still 0. Please provide a description of correct features greater than 0.')
            log.warn('Klammer will be implicitly off!')
        }

        // currently post-processing-tdc is always set since we do not support separate TD databases
        """
        ## Percolator does not have a threads parameter. Set it via OpenMP env variable,
        ## to honor threads on clusters
        OMP_NUMBER_THREADS=${task.cpus} PercolatorAdapter \\
                            -in ${id_file} \\
                            -out ${id_file.baseName}_perc.idXML \\
                            -threads ${task.cpus} \\
                            -subset-max-train ${params.subset_max_train} \\
                            -decoy-pattern ${params.decoy_affix} \\
                            -post-processing-tdc \\
                            > ${id_file.baseName}_percolator.log
        """
}

process idfilter {

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'
    publishDir "${params.outdir}/ids", mode: 'copy', pattern: '*.idXML'

    input:
     set mzml_id, file(id_file) from id_files_idx_feat_perc

    output:
     set mzml_id, file("${id_file.baseName}_filter.idXML") into id_files_idx_feat_perc_filter
     file "*.log"

    when:
     params.posterior_probabilities == "percolator"

    script:
     """
     IDFilter -in ${id_file} \\
                        -out ${id_file.baseName}_filter.idXML \\
                        -threads ${task.cpus} \\
                        -score:pep ${params.psm_pep_fdr_cutoff} \\
                        > ${id_file.baseName}_idfilter.log
     """
}

process idscoreswitcher {

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
     set mzml_id, file(id_file) from id_files_idx_feat_perc_filter

    output:
     set mzml_id, file("${id_file.baseName}_switched.idXML") into id_files_idx_feat_perc_fdr_filter_switched_luciphor, id_files_idx_feat_perc_fdr_filter_switched_plfq
     file "*.log"

    when:
     params.posterior_probabilities == "percolator"

    script:
     """
     IDScoreSwitcher    -in ${id_file} \\
                        -out ${id_file.baseName}_switched.idXML \\
                        -threads ${task.cpus} \\
                        -old_score q-value \\
                        -new_score MS:1001493 \\
                        -new_score_orientation lower_better \\
                        -new_score_type "Posterior Error Probability" \\
                        > ${id_file.baseName}_scoreswitcher.log
     """
}

// ---------------------------------------------------------------------
// Branch b) Q-values and PEP from OpenMS

// Note: for IDPEP we never need any file specific settings so we can stop adding the mzml_idto the channels
process fdr_idpep {

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
     set mzml_id, file(id_file) from id_files_idx_ForIDPEP

    output:
     set mzml_id, file("${id_file.baseName}_fdr.idXML") into id_files_idx_ForIDPEP_fdr
     file "*.log"

    when:
     params.posterior_probabilities != "percolator"

    script:
     """
     FalseDiscoveryRate -in ${id_file} \\
                        -out ${id_file.baseName}_fdr.idXML \\
                        -threads ${task.cpus} \\
                        -protein false \\
                        -algorithm:add_decoy_peptides \\
                        -algorithm:add_decoy_proteins \\
                        > ${id_file.baseName}_fdr.log
     """
}

process idscoreswitcher_idpep_pre {

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
     set mzml_id, file(id_file) from id_files_idx_ForIDPEP_fdr

    output:
     set mzml_id, file("${id_file.baseName}_switched.idXML") into id_files_idx_ForIDPEP_fdr_switch
     file "*.log"

    when:
     params.posterior_probabilities != "percolator"

    script:
     """
     IDScoreSwitcher    -in ${id_file} \\
                        -out ${id_file.baseName}_switched.idXML \\
                        -threads ${task.cpus} \\
                        -old_score q-value \\
                        -new_score ${search_engine_score}_score \\
                        -new_score_orientation lower_better \\
                        -new_score_type ${search_engine_score} \\
                        > ${id_file.baseName}_scoreswitcher1.log
     """
}

process idpep {

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
     set mzml_id, file(id_file) from id_files_idx_ForIDPEP_fdr_switch

    output:
     set mzml_id, file("${id_file.baseName}_idpep.idXML") into id_files_idx_ForIDPEP_fdr_switch_idpep
     file "*.log"

    when:
     params.posterior_probabilities != "percolator"

    script:
     """
     IDPosteriorErrorProbability    -in ${id_file} \\
                                    -out ${id_file.baseName}_idpep.idXML \\
                                    -threads ${task.cpus} \\
                                    > ${id_file.baseName}_idpep.log
     """
}

process idscoreswitcher_idpep_post {

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
     set mzml_id, file(id_file) from id_files_idx_ForIDPEP_fdr_switch_idpep

    output:
     set mzml_id, file("${id_file.baseName}_switched.idXML") into id_files_idx_ForIDPEP_fdr_switch_idpep_switch
     file "*.log"

    when:
     params.posterior_probabilities != "percolator"

    script:
     """
     IDScoreSwitcher    -in ${id_file} \\
                        -out ${id_file.baseName}_switched.idXML \\
                        -threads ${task.cpus} \\
                        -old_score "Posterior Error Probability" \\
                        -new_score q-value \\
                        -new_score_orientation lower_better \\
                        > ${id_file.baseName}_scoreswitcher2.log
     """
}

process idfilter_idpep {

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'
    publishDir "${params.outdir}/ids", mode: 'copy', pattern: '*.idXML'

    input:
     set mzml_id, file(id_file) from id_files_idx_ForIDPEP_fdr_switch_idpep_switch

    output:
     set mzml_id, file("${id_file.baseName}_filter.idXML") into id_files_idx_ForIDPEP_fdr_switch_idpep_switch_filter
     file "*.log"

    when:
     params.posterior_probabilities != "percolator"

    script:
     """
     IDFilter -in ${id_file} \\
              -out ${id_file.baseName}_filter.idXML \\
              -threads ${task.cpus} \\
              -score:pep ${params.psm_pep_fdr_cutoff} \\
              > ${id_file.baseName}_idfilter1.log
     """
}

process idscoreswitcher_idpep_postfilter {

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
     set mzml_id, file(id_file) from id_files_idx_ForIDPEP_fdr_switch_idpep_switch_filter

    output:
     set mzml_id, file("${id_file.baseName}_switched.idXML") into id_files_idx_ForIDPEP_fdr_switch_idpep_switch_filter_switch_plfq, id_files_idx_ForIDPEP_fdr_switch_idpep_switch_filter_switch_luciphor
     file "*.log"

    when:
     params.posterior_probabilities != "percolator"

    script:
     """
     IDScoreSwitcher    -in ${id_file} \\
                        -out ${id_file.baseName}_switched.idXML \\
                        -threads ${task.cpus} \\
                        -old_score q-value \\
                        -new_score "Posterior Error Probability" \\
                        -new_score_orientation lower_better \\
                        > ${id_file.baseName}_scoreswitcher3.log
     """
}

process luciphor {

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
     tuple mzml_id, file(mzml_file), file(id_file), frag_method from mzmls_luciphor.join(id_files_idx_feat_perc_fdr_filter_switched_luciphor.mix(id_files_idx_ForIDPEP_fdr_switch_idpep_switch_filter_switch_luciphor)).join(ch_sdrf_config.luciphor_settings)

    output:
     set mzml_id, file("${id_file.baseName}_luciphor.idXML") into id_files_luciphor
     file "*.log"

    when:
     params.enable_mod_localization

    script:
     id_files_idx_ForIDPEP_fdr_switch_idpep_switch_filter_switch_plfq = Channel.empty()
     id_files_idx_feat_perc_fdr_filter_switched_plfq = Channel.empty()
     def losses = params.luciphor_neutral_losses ? '-neutral_losses "${params.luciphor_neutral_losses}"' : ''
     def dec_mass = params.luciphor_decoy_mass ? '-decoy_mass "${params.luciphor_decoy_mass}"' : ''
     def dec_losses = params.luciphor_decoy_neutral_losses ? '-decoy_neutral_losses "${params.luciphor_decoy_neutral_losses}' : ''
     """
     LuciphorAdapter    -id ${id_file} \\
                        -in ${mzml_file} \\
                        -out ${id_file.baseName}_luciphor.idXML \\
                        -threads ${task.cpus} \\
                        -num_threads ${task.cpus} \\
                        -target_modifications ${params.mod_localization.tokenize(',').collect { "'${it}'" }.join(" ") } \\
                        -fragment_method ${frag_method} \\
                        ${losses} \\
                        ${dec_mass} \\
                        ${dec_losses} \\
                        -max_charge_state ${params.max_precursor_charge} \\
                        -max_peptide_length ${params.max_peptide_length} \\
                        > ${id_file.baseName}_scoreswitcher.log
     """
                     //        -fragment_mass_tolerance ${} \\
                     //   -fragment_error_units ${} \\
}

// ---------------------------------------------------------------------
// Main Branch

// Join mzmls and ids by UID specified per mzml file in the beginning.
// ID files can come directly from the Percolator branch, IDPEP branch or
// after optional processing with Luciphor
mzmls_plfq
  .join(id_files_luciphor
        .mix(id_files_idx_ForIDPEP_fdr_switch_idpep_switch_filter_switch_plfq)
        .mix(id_files_idx_feat_perc_fdr_filter_switched_plfq))
  .multiMap{ it ->
      mzmls: it[1]
      ids: it[2]
  }
  .set{ch_plfq}

process proteomicslfq {

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'
    publishDir "${params.outdir}/proteomics_lfq", mode: 'copy'

    ///.toSortedList({ a, b -> b.baseName <=> a.baseName })
    input:
     file(mzmls) from ch_plfq.mzmls.collect().view()
     file(id_files) from ch_plfq.ids.collect().view()
     file expdes from ch_expdesign
     file fasta from plfq_in_db.mix(plfq_in_db_decoy)

    output:
     file "out.mzTab" into out_mzTab
     file "out.consensusXML" into out_consensusXML
     file "out.csv" into out_msstats
     file "debug_mergedIDs.idXML" optional true
     file "debug_mergedIDs_inference.idXML" optional true
     file "debug_mergedIDsGreedyResolved.idXML" optional true
     file "debug_mergedIDsGreedyResolvedFDR.idXML" optional true
     file "debug_mergedIDsGreedyResolvedFDRFiltered.idXML" optional true
     file "debug_mergedIDsFDRFilteredStrictlyUniqueResolved.idXML" optional true
     file "*.log"

    script:
     """
     ProteomicsLFQ -in ${(mzmls as List).join(' ')} \\
                    -ids ${(id_files as List).join(' ')} \\
                    -design ${expdes} \\
                    -fasta ${fasta} \\
                    -protein_inference ${params.protein_inference} \\
                    -targeted_only ${params.targeted_only} \\
                    -mass_recalibration ${params.mass_recalibration} \\
                    -transfer_ids ${params.transfer_ids} \\
                    -out out.mzTab \\
                    -threads ${task.cpus} \\
                    -out_msstats out.csv \\
                    -out_cxml out.consensusXML \\
                    -proteinFDR ${params.protein_level_fdr_cutoff} \\
                    -debug ${params.inf_quant_debug} \\
                    > proteomicslfq.log
     """
}


// TODO the script supports a control condition as third argument
// TODO the second argument can be "pairwise" or TODO later a user defined contrast string

process msstats {

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'
    publishDir "${params.outdir}/msstats", mode: 'copy'

    when:
     !params.skip_post_msstats

    input:
     file csv from out_msstats

    output:
     file "*.pdf"
     file "*.csv"
     file "*.log"

    script:
     """
     msstats_plfq.R ${csv} > msstats.log || echo "Optional MSstats step failed. Please check logs and re-run or do a manual statistical analysis."
     """
}

//TODO allow user config yml (as second arg to the script

process ptxqc {

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'
    publishDir "${params.outdir}/ptxqc", mode: 'copy'

    when:
     params.enable_qc

    input:
     file mzTab from out_mzTab

    output:
     file "*.html"
     file "*.yaml"
     file "*.Rmd"
     file "*.pdf"
     file "*.txt"

    script:
     """
     ptxqc.R ${mzTab} > ptxqc.log
     """
}
  


//--------------------------------------------------------------- //
//---------------------- Nextflow specifics --------------------- //
//--------------------------------------------------------------- //


// Header log info
log.info nfcoreHeader()
def summary = [:]
summary['Run Name']         = custom_runName ?: workflow.runName
// TODO nf-core: Report custom parameters here
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if(workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if(workflow.profile == 'awsbatch'){
   summary['AWS Region']    = params.awsregion
   summary['AWS Queue']     = params.awsqueue
}
summary['Config Profile'] = workflow.profile
if(params.config_profile_description) summary['Config Description'] = params.config_profile_description
if(params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if(params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if(params.email) {
  summary['E-mail Address']  = params.email
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "\033[2m----------------------------------------------------\033[0m"

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-proteomicslfq-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/proteomicslfq Workflow Summary'
    section_href: 'https://github.com/nf-core/proteomicslfq'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }


/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy',
        saveAs: { filename ->
                      if (filename.indexOf(".csv") > 0) filename
                      else null
                }

    output:
    file 'software_versions_mqc.yaml' into ch_software_versions_yaml
    file "software_versions.csv"

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    ThermoRawFileParser.sh --version &> v_thermorawfileparser.txt
    echo \$(FileConverter 2>&1) > v_fileconverter.txt || true
    echo \$(DecoyDatabase 2>&1) > v_decoydatabase.txt || true
    echo \$(MSGFPlusAdapter 2>&1) > v_msgfplusadapter.txt || true
    echo \$(msgf_plus 2>&1) > v_msgfplus.txt || true
    echo \$(CometAdapter 2>&1) > v_cometadapter.txt || true
    echo \$(comet 2>&1) > v_comet.txt || true
    echo \$(PeptideIndexer 2>&1) > v_peptideindexer.txt || true
    echo \$(PSMFeatureExtractor 2>&1) > v_psmfeatureextractor.txt || true
    echo \$(PercolatorAdapter 2>&1) > v_percolatoradapter.txt || true
    percolator -h &> v_percolator.txt
    echo \$(IDFilter 2>&1) > v_idfilter.txt || true
    echo \$(IDScoreSwitcher 2>&1) > v_idscoreswitcher.txt || true
    echo \$(FalseDiscoveryRate 2>&1) > v_falsediscoveryrate.txt || true
    echo \$(IDPosteriorErrorProbability 2>&1) > v_idposteriorerrorprobability.txt || true
    echo \$(ProteomicsLFQ 2>&1) > v_proteomicslfq.txt || true
    echo $workflow.manifest.version &> v_msstats_plfq.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

/*
 * STEP 3 - Output Description HTML
 */
/* TODO Deactivated for now
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    file output_docs from ch_output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}
*/


/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/proteomicslfq] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[nf-core/proteomicslfq] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/proteomicslfq] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            [ 'mail', '-s', subject, email_address ].execute() << email_txt
            log.info "[nf-core/proteomicslfq] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/proteomicslfq]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/proteomicslfq]${c_red} Pipeline completed with errors${c_reset}-"
    }

}


def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/proteomicslfq v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}

// Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}
