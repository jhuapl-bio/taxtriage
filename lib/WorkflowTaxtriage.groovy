//
// This file holds several functions specific to the workflow/taxtriage.nf in the nf-core/taxtriage pipeline
//

class WorkflowTaxtriage {

    //
    // SRA / ENA accession detection.
    //
    // Single source of truth for "is this fastq_1 value an accession rather than a
    // path?", used by workflows/taxtriage.nf (to skip its file-existence pre-flight
    // check) and subworkflows/local/input_check.nf (to route the row through the
    // download path). Mirrors the same patterns in bin/resolve_sra.py and
    // bin/check_samplesheet.py.
    //
    // Anything containing a path separator, dot or ';' is treated as a path, never
    // an accession, so a local file named "SRR13191702.fastq.gz" still resolves
    // normally.
    //
    public static boolean isAccession(value) {
        if (!value) {
            return false
        }
        def v = value.toString().trim()
        if (v.contains('/') || v.contains('\\') || v.contains('.') || v.contains(';')) {
            return false
        }
        return v ==~ /(?i)^(SRR|ERR|DRR|SRX|ERX|DRX|SRS|ERS|DRS|SRP|ERP|DRP)\d{5,}$/ ||
               v ==~ /(?i)^SAM(N|EA|EG|D)\d+$/ ||
               v ==~ /(?i)^PRJ(NA|EB|EA|DA|DB)\d+$/
    }

    //
    // Check and validate parameters
    //
    //
    // ── Where auto-downloaded databases are cached ────────────────────────────
    //
    // The caches are `storeDir` targets: the download happens once and every
    // later run that resolves to the same directory skips the process outright.
    // That only works if the directory is (a) writable by the executor and
    // (b) the SAME path on the next run.
    //
    // `${projectDir}` satisfies neither on Seqera / cloud executors. When the
    // pipeline is pulled rather than checked out, projectDir is the local clone
    // of one commit:
    //
    //     /.nextflow/assets/.repos/jhuapl-bio/taxtriage/clones/<sha>/dbs/kaiju
    //
    // which lives inside the head-job container (so a task running on another
    // node cannot write it, and nothing there survives the job), and is keyed by
    // COMMIT, so bumping the revision silently re-downloads everything. Worse,
    // on a cloud executor it is a local POSIX path where every other path is a
    // bucket URI, so the storeDir is simply lost.
    //
    // Resolution order:
    //   1. the per-backend override (--novelty_kaiju_db_cache, ...)   - explicit wins
    //   2. --db_cache_dir <base>/<kind>                               - one base for all
    //   3. <workDir>/dbs/<kind>   when the work dir is remote (s3://, gs://,
    //      az://) OR the pipeline is running from a pulled clone       - bucket-native,
    //      shared by every task, and stable across runs and revisions because the
    //      work dir is what the user (or Seqera) pins
    //   4. <projectDir>/dbs/<kind>                                    - a plain local
    //      checkout, i.e. the historical behaviour, unchanged
    //
    public static String dbCacheDir(params, workflow, String kind) {
        def override = null
        switch (kind) {
            case 'mmseqs':  override = params.novelty_db_cache;         break
            case 'kaiju':   override = params.novelty_kaiju_db_cache;   break
            case 'kraken2': override = params.novelty_kraken2_db_cache; break
        }
        if (override) {
            return override.toString()
        }
        if (params.db_cache_dir) {
            return "${params.db_cache_dir}/${kind}".toString()
        }

        def workDir    = workflow.workDir.toString()
        def projectDir = workflow.projectDir.toString()
        boolean isRemote = (workDir =~ /^[a-zA-Z0-9+.-]+:\/\//).find()
        // Nextflow stages a pulled pipeline under <home>/.nextflow/assets/... ,
        // and Seqera adds a per-commit `clones/<sha>` level. Either way the path
        // is not a stable, shared, writable place to park tens of GB.
        def isPulledClone = projectDir.contains('/.nextflow/assets/') || projectDir.contains('/clones/')

        return (isRemote || isPulledClone)
            ? "${workDir}/dbs/${kind}".toString()
            : "${projectDir}/dbs/${kind}".toString()
    }

    public static void initialise(params, log) {
        genomeExistsError(params, log)
        mergeHostTaxids(params, log)

        // if (!params.fasta) {
        //     log.error "Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file."
        //     System.exit(1)
        // }
    }

    //
    // Named host targets (conf/hosts.config) declare the taxids that go with the
    // genomes they de-host against.  Read removal is never perfect, so those
    // taxids are also merged into:
    //   --remove_taxids                 -> dropped from the classification report
    //   --report_flag_exclude_taxids    -> never counted as a detection
    // so a handful of surviving host reads cannot show up as an organism call.
    //
    // A value the user explicitly set to an empty string is left alone: that is
    // how you opt out ("count host like any other organism").
    //
    public static void mergeHostTaxids(params, log) {
        if (!params.genome || !params.genomes || !params.genomes.containsKey(params.genome)) {
            return
        }
        def entry = params.genomes[params.genome]
        if (!entry || !entry.taxids) {
            return
        }
        def host_taxids = entry.taxids.toString().split(/[\s,]+/).findAll { it }
        if (!host_taxids) {
            return
        }

        def merge = { current ->
            def existing = (current == null) ? [] : current.toString().split(/[\s,]+/).findAll { it }
            return (existing + host_taxids).unique().join(' ')
        }

        if (!(params.remove_taxids instanceof String && params.remove_taxids.trim() == '')) {
            params.remove_taxids = merge(params.remove_taxids)
        }
        if (!(params.report_flag_exclude_taxids instanceof String && params.report_flag_exclude_taxids.trim() == '')) {
            // null here means "the report's own default (9606)"; make that explicit
            // before adding the host taxids so human is not silently dropped.
            def base = (params.report_flag_exclude_taxids == null) ? '9606' : params.report_flag_exclude_taxids
            params.report_flag_exclude_taxids = merge(base)
        }

        log.info "Host target '${params.genome}': taxids ${host_taxids.join(', ')} merged into " +
                 "--remove_taxids ('${params.remove_taxids}') and " +
                 "--report_flag_exclude_taxids ('${params.report_flag_exclude_taxids}')."
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    //
    // Exit pipeline if incorrect --genome key provided
    //
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }
    }
}
