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

    Mandatory arguments:
      --spectra                     Path to input spectra as mzML or Thermo Raw
      --database                    Path to input protein database as fasta


    Database Search:
      --search_engine               Which search engine: "comet" or "msgf"
      --enzyme                      Enzymatic cleavage ('unspecific cleavage', 'Trypsin', see OpenMS enzymes)
      --fixed_mods                  Fixed modifications ('Carbamidomethyl (C)', see OpenMS modifications)
      --variable_mods               Variable modifications ('Oxidation (M)', see OpenMS modifications)
      --precursor_mass_tolerance    Mass tolerance of precursor mass (ppm)
      --allowed_missed_cleavages    Allowed missed cleavages 
      --psm_level_fdr_cutoff        Identification PSM-level FDR cutoff
      --posterior_probabilities     How to calculate posterior probabilities for PSMs:
                                    "percolator" = Re-score based on PSM-feature-based SVM and transform distance
                                        to hyperplane for posteriors 
                                    "fit_distributions" = Fit positive and negative distributions to scores
                                        (similar to PeptideProphet)


    Inference:
      --protein_inference           Infer proteins through:
                                    "aggregation"  = aggregates all peptide scores across a protein (by calculating the maximum)
                                    "bayesian"     = computes a posterior probability for every protein based on a Bayesian network
      --protein_level_fdr_cutoff    Identification protein-level FDR cutoff

    Quantification:
      --transfer_ids                Transfer IDs over aligned samples to increase # of quantifiable features (WARNING:
                                    increased memory consumption)
      --targeted_only               Only ID based quantification
      --mass_recalibration          Recalibrates masses to correct for instrument biases
      --protein_quantification      Quantify proteins based on:
                                    "unique_peptides" = use peptides mapping to single proteins or a group of indistinguishable proteins (according to the set of experimentally identified peptides)
                                    "strictly_unique_peptides" = use peptides mapping to a unique single protein only
                                    "shared_peptides" = use shared peptides only for its best group (by inference score)

    General Options:
      --expdesign                   Path to experimental design file (if not given, it assumes unfractionated, unrelated samples)
      --add_decoys                  Add decoys to the given fasta

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
params.spectra = params.spectra ?: { log.error "No spectra data provided. Make sure you have used the '--spectra' option."; exit 1 }()
params.database = params.database ?: { log.error "No protein database provided. Make sure you have used the '--database' option."; exit 1 }()
// params.expdesign = params.expdesign ?: { log.error "No read data privided. Make sure you have used the '--design' option."; exit 1 }()
params.outdir = params.outdir ?: { log.warn "No output directory provided. Will put the results into './results'"; return "./results" }()

/*
 * Create a channel for input read files
 */

ch_spectra = Channel.fromPath(params.spectra, checkIfExists: true)
ch_database = Channel.fromPath(params.database).set{ db_for_decoy_creation }
// ch_expdesign = Channel.fromPath(params.design, checkIfExists: true)

//use a branch operator for this sort of thing and access the files accordingly!

ch_spectra
.branch {
        raw: hasExtension(it, 'raw')
        mzML: hasExtension(it, 'mzML')
}
.set {branched_input}


//TODO we could also check for outdated mzML versions and try to update them
branched_input.mzML
.branch {
        nonIndexedMzML: file(it).withReader {
                            f = it;
                            1.upto(5) {
                              if (f.readLine().contains("indexedmzML")) return false;
                            }
                            return true;
                        }
        inputIndexedMzML: file(it).withReader {
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
// - Check how to avoid copying of the database for example (currently we get one copy for each SE run). Is it the
//   "each file()" pattern I used?

/*
 * STEP 0.1 - Raw file conversion
 */
process raw_file_conversion {

    input:
        file rawfile from branched_input.raw

    output:
        file "*.mzML" into mzmls_converted
    
    // TODO use actual ThermoRawfileConverter!!
    script:
        """
        mv ${rawfile} ${rawfile.baseName}.mzML
        """
}

/*
 * STEP 0.2 - MzML indexing
 */
process mzml_indexing {

    input:
        file mzmlfile from branched_input_mzMLs.nonIndexedMzML

    output:
        file "out/*.mzML" into mzmls_indexed
    
    script:
        """
        mkdir out
        FileConverter -in ${mzmlfile} -out out/${mzmlfile.baseName}.mzML
        """
}

//Mix the converted raw data with the already supplied mzMLs and push these to the same channels as before

branched_input_mzMLs.inputIndexedMzML.mix(mzmls_converted).mix(mzmls_indexed).into{mzmls; mzmls_plfq}


if (params.expdesign)
{
    Channel
        .fromPath(params.expdesign)
        .ifEmpty { exit 1, "params.expdesign was empty - no input files supplied" }
        .set { expdesign }
}


//Fill the channels with empty Channels in case that we want to add decoys. Otherwise fill with output from database.
(searchengine_in_db, pepidx_in_db, plfq_in_db) = ( params.add_decoys
                    ? [ Channel.empty(), Channel.empty(), Channel.empty() ]
                    : [ Channel.fromPath(params.database),Channel.fromPath(params.database), Channel.fromPath(params.database)  ] )   

//Add decoys if params.add_decoys is set appropriately
process generate_decoy_database {

    input:
        file(mydatabase) from db_for_decoy_creation

    output:
        file "${database.baseName}_decoy.fasta" into searchengine_in_db_decoy, pepidx_in_db_decoy, plfq_in_db_decoy
        //TODO need to add these channel with .mix(searchengine_in_db_decoy) for example to all subsequent processes that need this...

    when: params.add_decoys

    script:
        """
        DecoyDatabase  -in ${mydatabase} \\
                    -out ${mydatabase.baseName}_decoy.fasta \\
                    -decoy_string DECOY_ \\
                    -decoy_string_position prefix
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

/// Search engine
// TODO parameterize more
if (params.search_engine == "msgf")
{
   search_engine_score = "SpecEValue"

    process search_engine_msgf {
        echo true
        input:
         tuple file(database), file(mzml_file) from searchengine_in_db.mix(searchengine_in_db_decoy).combine(mzmls)
         
         // This was another way of handling the combination
         //file database from searchengine_in_db.mix(searchengine_in_db_decoy)
         //each file(mzml_file) from mzmls


        output:
         file "${mzml_file.baseName}.idXML" into id_files
     
        script:
         """
         MSGFPlusAdapter  -in ${mzml_file} \\
                       -out ${mzml_file.baseName}.idXML \\
                       -threads ${task.cpus} \\
                       -database ${database}
         """
     }

} else {

    search_engine_score = "expect"

    process search_engine_comet {
        echo true
        input:
         file database from searchengine_in_db.mix(searchengine_in_db_decoy)
         each file(mzml_file) from mzmls

        output:
         file "${mzml_file.baseName}.idXML" into id_files
     
        script:
         """
         CometAdapter  -in ${mzml_file} \\
                       -out ${mzml_file.baseName}.idXML \\
                       -threads ${task.cpus} \\
                       -database ${database}
         """
     }
}



process index_peptides {
    echo true
    input:
     each file(id_file) from id_files
     file database from pepidx_in_db.mix(pepidx_in_db_decoy)
     
    output:
     file "${id_file.baseName}_idx.idXML" into id_files_idx_ForPerc, id_files_idx_ForIDPEP

    script:
     """
     PeptideIndexer -in ${id_file} \\
                    -out ${id_file.baseName}_idx.idXML \\
                    -threads ${task.cpus} \\
                    -fasta ${database}
     """

}


// ---------------------------------------------------------------------
// Branch a) Q-values and PEP from Percolator


process extract_perc_features {
 
    input:
     file id_file from id_files_idx_ForPerc

    output:
     file "${id_file.baseName}_feat.idXML" into id_files_idx_feat

    when:
     params.posterior_probabilities == "percolator"

    script:
     """
     PSMFeatureExtractor -in ${id_file} \\
                        -out ${id_file.baseName}_feat.idXML \\
                        -threads ${task.cpus}
     """

}

//TODO parameterize and find a way to run across all runs merged
process percolator {
 
    input:
     file id_file from id_files_idx_feat

    output:
     file "${id_file.baseName}_perc.idXML" into id_files_idx_feat_perc

    when:
     params.posterior_probabilities == "percolator"

    script:
     """
     PercolatorAdapter -in ${id_file} \\
                        -out ${id_file.baseName}_perc.idXML \\
                        -threads ${task.cpus} \\
                        -post-processing-tdc -subset-max-train 100000 -decoy-pattern "rev"
     """

}

process idfilter {
 
    publishDir "${params.outdir}/ids", mode: 'copy'

    input:
     file id_file from id_files_idx_feat_perc

    output:
     file "${id_file.baseName}_filter.idXML" into id_files_idx_feat_perc_filter

    when:
     params.posterior_probabilities == "percolator"

    script:
     """
     IDFilter -in ${id_file} \\
                        -out ${id_file.baseName}_filter.idXML \\
                        -threads ${task.cpus} \\
                        -score:pep ${params.psm_level_fdr_cutoff}
     """

}

process idscoreswitcher {
 
    input:
     file id_file from id_files_idx_feat_perc_filter

    output:
     file "${id_file.baseName}_switched.idXML" into id_files_idx_feat_perc_fdr_filter_switched

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
                        -new_score_type "Posterior Error Probability"
     """

}



// ---------------------------------------------------------------------
// Branch b) Q-values and PEP from OpenMS

process fdr {
 
    input:
     file id_file from id_files_idx_ForIDPEP

    output:
     file "${id_file.baseName}_fdr.idXML" into id_files_idx_ForIDPEP_fdr

    when:
     params.posterior_probabilities != "percolator"

    script:
     """
     FalseDiscoveryRate -in ${id_file} \\
                        -out ${id_file.baseName}_fdr.idXML \\
                        -threads ${task.cpus} \\
                        -protein false -algorithm:add_decoy_peptides -algorithm:add_decoy_proteins
     """

}

process idscoreswitcher1 {
 
    input:
     file id_file from id_files_idx_ForIDPEP_fdr

    output:
     file "${id_file.baseName}_switched.idXML" into id_files_idx_ForIDPEP_fdr_switch

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
                        -new_score_type ${search_engine_score}
     """

}

//TODO probably not needed when using Percolator. You can use the qval from there
process idpep {
 
    input:
     file id_file from id_files_idx_ForIDPEP_fdr_switch

    output:
     file "${id_file.baseName}_idpep.idXML" into id_files_idx_ForIDPEP_fdr_switch_idpep

    when:
     params.posterior_probabilities != "percolator"

    script:
     """
     IDPosteriorErrorProbability    -in ${id_file} \\
                                    -out ${id_file.baseName}_idpep.idXML \\
                                    -threads ${task.cpus}
     """

}

process idscoreswitcher2 {
 
    input:
     file id_file from id_files_idx_ForIDPEP_fdr_switch_idpep

    output:
     file "${id_file.baseName}_switched.idXML" into id_files_idx_ForIDPEP_fdr_switch_idpep_switch

    when:
     params.posterior_probabilities != "percolator"

    script:
     """
     IDScoreSwitcher    -in ${id_file} \\
                        -out ${id_file.baseName}_switched.idXML \\
                        -threads ${task.cpus} \\
                        -old_score "Posterior Error Probability" \\
                        -new_score q-value \\
                        -new_score_orientation lower_better
     """

}

process idfilter2 {
 
    publishDir "${params.outdir}/ids", mode: 'copy'

    input:
     file id_file from id_files_idx_ForIDPEP_fdr_switch_idpep_switch

    output:
     file "${id_file.baseName}_filter.idXML" into id_files_idx_ForIDPEP_fdr_switch_idpep_switch_filter

    when:
     params.posterior_probabilities != "percolator"
      
    script:
     """
     IDFilter -in ${id_file} \\
                        -out ${id_file.baseName}_filter.idXML \\
                        -threads ${task.cpus} \\
                        -score:pep ${params.psm_level_fdr_cutoff}
     """

}

process idscoreswitcher3 {
 
    input:
     file id_file from id_files_idx_ForIDPEP_fdr_switch_idpep_switch_filter

    output:
     file "${id_file.baseName}_switched.idXML" into id_files_idx_ForIDPEP_fdr_switch_idpep_switch_filter_switch

    when:
     params.posterior_probabilities != "percolator"

    script:
     """
     IDScoreSwitcher    -in ${id_file} \\
                        -out ${id_file.baseName}_switched.idXML \\
                        -threads ${task.cpus} \\
                        -old_score q-value \\
                        -new_score "Posterior Error Probability" \\
                        -new_score_orientation lower_better
     """

}


// ---------------------------------------------------------------------
// Main Branch

process proteomicslfq {
 
    publishDir "${params.outdir}/proteomics_lfq", mode: 'copy'
    
    input:
     file mzmls from mzmls_plfq.toSortedList({ a, b -> b.baseName <=> a.baseName }).view()
     file id_files from id_files_idx_feat_perc_fdr_filter_switched
         .mix(id_files_idx_ForIDPEP_fdr_switch_idpep_switch_filter_switch)
         .toSortedList({ a, b -> b.baseName <=> a.baseName })
         .view()
     file expdes from expdesign
     file fasta from plfq_in_db.mix(plfq_in_db_decoy)

    output:
     file "out.mzTab" into out_mzTab
     file "out.consensusXML" into out_consensusXML
     file "out.csv" into out_msstats
     file "debug_mergedIDs.idXML" into debug_id
     file "debug_mergedIDs_inference.idXML" into debug_id_inf
     file "debug_mergedIDsGreedyResolved.idXML" into debug_id_resolve

    script:
     """
     ProteomicsLFQ -in ${(mzmls as List).join(' ')} \\
                    -ids ${(id_files as List).join(' ')} \\
                    -design ${expdes} \\
                    -fasta ${fasta} \\
                    -targeted_only "true" \\
                    -mass_recalibration "false" \\
                    -transfer_ids "false" \\
                    -out out.mzTab \\
                    -threads ${task.cpus} \\
                    -out_msstats out.csv \\
                    -out_cxml out.consensusXML \\
                    -debug 667

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
  summary['MultiQC maxsize'] = params.maxMultiqcEmailFileSize
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "\033[2m----------------------------------------------------\033[0m"

// Check the hostnames against configured profiles
checkHostname()

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-proteomicslfq-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/proteomicslfq Workflow Summary'
    section_href: 'https://github.com/nf-core/proteomicslfq'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}


/*
 * Parse software version numbers
 */
process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    // TODO nf-core: Get all tools to print their version number here
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    echo "foo" > software_versions_mqc.yaml
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
    markdown_to_html.r $output_docs results_description.html
    """
}
*/


/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/proteomicslfq] Successful: $workflow.runName"
    if(!workflow.success){
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
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if(workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // TODO nf-core: If not using MultiQC, strip out this code (including params.maxMultiqcEmailFileSize)
    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList){
                log.warn "[nf-core/proteomicslfq] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/proteomicslfq] Could not attach MultiQC report to summary email"
    }

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
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.maxMultiqcEmailFileSize.toBytes() ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[nf-core/proteomicslfq] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/proteomicslfq] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    if(workflow.success){
        log.info "${c_purple}[nf-core/proteomicslfq]${c_green} Pipeline complete${c_reset}"
    } else {
        checkHostname()
        log.info "${c_purple}[nf-core/proteomicslfq]${c_red} Pipeline completed with errors${c_reset}"
    }

}


def nfcoreHeader(){
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";

    return """    ${c_dim}----------------------------------------------------${c_reset}
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/proteomicslfq v${workflow.manifest.version}${c_reset}
    ${c_dim}----------------------------------------------------${c_reset}
    """.stripIndent()
}

def checkHostname(){
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if(params.hostnames){
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if(hostname.contains(hname) && !workflow.profile.contains(prof)){
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
