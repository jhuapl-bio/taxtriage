//
// This file holds several functions specific to the workflow/taxtriage.nf in the nf-core/taxtriage pipeline
//

import nextflow.Nextflow

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
    // ── Empty read files ──────────────────────────────────────────────────────
    //
    // A sample whose reads were ALL host (or all filtered out by the kraken2
    // de-hosting filter) comes out of HOST_REMOVAL as an EMPTY fastq. Every read
    // QC tool downstream -- trimgalore, porechop, fastp -- exits non-zero on an
    // empty input, and with the default `finish` error strategy that one dead-end
    // sample used to take the whole run down with it. Such samples are dropped
    // before trimming instead (see workflows/taxtriage.nf), so the rest of the
    // run completes.
    //
    // Size, not record count, on purpose: the channel holds Paths that may live in
    // a bucket, so `.size()` is a cheap HEAD request while counting records would
    // mean pulling every fastq onto the head node. An empty gzip member is 20-30
    // bytes (gzip writes 20-25, bgzip's EOF block is 28); a single 20 bp read
    // already gzips to ~59 bytes and a usable sample is orders of magnitude bigger,
    // so the threshold below separates "no records at all" from "some data" without
    // ever having to decompress anything.
    //
    public static final int EMPTY_READS_BYTES = 100

    //
    // True when every read file in the entry holds at least one record. A pair with
    // one empty mate counts as empty: the mates would no longer be in sync, and the
    // paired tools fail on that just as hard as on a fully empty input.
    //
    public static boolean hasReads(reads) {
        if (!reads) { return false }
        def files = (reads instanceof List) ? reads : [reads]
        if (!files) { return false }
        return files.every { f ->
            if (!f) { return false }
            def sz = -1
            try { sz = f.size() } catch (Exception e) { return true }  // unreadable here -> let the task decide
            return f.toString().endsWith('.gz') ? sz > EMPTY_READS_BYTES : sz > 0
        }
    }

    //
    // ── Where auto-downloaded databases are cached ────────────────────────────
    //
    // The caches are `storeDir` targets: the download happens once and every
    // later run that resolves to the same directory skips the process outright.
    // That only works if the directory is (a) writable by the executor, (b) the
    // SAME path on the next run, and (c) a path the task launcher will accept.
    //
    // Nothing auto-derived satisfies all three off a developer laptop:
    //   * `${projectDir}` is a throwaway per-revision clone inside the head job
    //     when the pipeline is PULLED rather than checked out, so it is head-node
    //     only and keyed by commit.
    //   * `${workDir}` is node-local scratch on Seqera (/nftass/scratch/<id>),
    //     which the Fusion script launcher rejects outright:
    //       Unexpected path for Fusion script launcher: /nftass/scratch/.../dbs/kaiju/viral
    //
    // So the default is NO storeDir -- the same contract as the main --db
    // (DOWNLOAD_DB has no storeDir either, which is exactly why `--db <alias>`
    // has always worked everywhere). Caching is opt-in.
    //
    // Resolution order:
    //   1. the per-backend override (--novelty_kaiju_db_cache, ...)   - explicit wins
    //   2. --db_cache_dir <base>/<kind>                               - one base for all
    //   3. <projectDir>/dbs/<kind>   ONLY for a genuine local checkout - historical
    //      behaviour, so an existing local dbs/ folder is still picked up
    //   4. null -> no storeDir; the db lands in the task work dir, like --db
    //
    // Any candidate that is a plain local path while the work dir is a bucket URI
    // is dropped (see usableStoreDir): cloud workers, and Fusion in particular,
    // cannot reach it.
    //
    public static String dbCacheDir(params, workflow, String kind) {
        def override = null
        switch (kind) {
            case 'mmseqs':  override = params.novelty_db_cache;         break
            case 'kaiju':   override = params.novelty_kaiju_db_cache;   break
            case 'kraken2': override = params.novelty_kraken2_db_cache; break
        }
        if (override) {
            return usableStoreDir(override.toString(), workflow)
        }
        if (params.db_cache_dir) {
            return usableStoreDir("${params.db_cache_dir}/${kind}".toString(), workflow)
        }

        // NO IMPLICIT CACHE. This is the --db behaviour, deliberately: DOWNLOAD_DB has
        // no storeDir at all, the db lands in the task work dir, and that is why
        // `--db <alias>` works on every executor (Seqera/Fusion, AWS Batch, local)
        // while an auto-derived `--novelty_db` cache did not.
        //
        // The old default was `<workDir>/dbs/<kind>` when projectDir was not a local
        // checkout. On Seqera the head job's work dir is a node-local scratch path
        // (e.g. /nftass/scratch/<id>), so that produced a storeDir like
        // `/nftass/scratch/<id>/dbs/kaiju/viral` -- an absolute POSIX path that the
        // Fusion script launcher rejects outright:
        //     Unexpected path for Fusion script launcher: /nftass/scratch/.../dbs/kaiju/viral
        // Fusion only accepts paths under its own mount (bucket-backed), so any local
        // path handed to it kills the task before it starts.
        //
        // A cache across runs is now strictly opt-in: --db_cache_dir (or the
        // per-backend --novelty_{,kaiju_,kraken2_}db_cache), pointed at something the
        // executor can actually reach -- an s3://... prefix on cloud, a shared mount
        // on a cluster, a plain folder locally.
        //
        // The one exception is a genuine local checkout, where `<projectDir>/dbs/<kind>`
        // has always worked and an existing dbs/ folder should keep being picked up.
        def projectDir = workflow.projectDir.toString()
        if (isLocalCheckout(projectDir)) {
            return usableStoreDir("${projectDir}/dbs/${kind}".toString(), workflow)
        }
        return null
    }

    //
    // A storeDir is only usable if the TASKS can reach it. When the work dir is a
    // bucket URI (s3://, gs://, az://) the tasks run on cloud workers -- possibly
    // behind Fusion, which refuses any path outside its own mount -- so a plain
    // local POSIX path is not a cache, it is a hard failure. Drop it and let the
    // download land in the work dir, exactly like DOWNLOAD_DB (--db) does.
    //
    private static String usableStoreDir(String dir, workflow) {
        if (!dir) { return null }
        def isRemote = { String p -> p ==~ /^[a-zA-Z0-9+.-]+:\/\/.*/ && !p.startsWith('file://') }
        if (isRemote(workflow.workDir.toString()) && !isRemote(dir)) {
            return null
        }
        return dir
    }

    //
    // True only for a genuine local clone of the pipeline the user controls -- never
    // for the copy Nextflow (or Seqera) stages when the pipeline is pulled by name.
    //
    private static boolean isLocalCheckout(String projectDir) {
        if (!projectDir || projectDir =~ /^[a-zA-Z0-9+.-]+:\/\//) {
            return false
        }
        def nxfHome = System.getenv('NXF_HOME') ?: "${System.getProperty('user.home')}/.nextflow".toString()
        def dir = new File(projectDir)
        def canon = null
        try { canon = dir.canonicalPath } catch (Exception e) { canon = projectDir }
        if (canon.startsWith(new File(nxfHome).absolutePath) ||
            canon.contains('/assets/') || canon.contains('/clones/') || canon.contains('/.nextflow/')) {
            return false
        }
        return dir.isDirectory() && dir.canWrite() &&
               (new File(dir, '.git').exists() || new File(dir, 'nextflow.config').canWrite())
    }

    //
    // The storeDir value for an auto-downloaded db, or null when no safe store
    // directory exists (see dbCacheDir). Modules call this directly so a null
    // cache dir drops the directive instead of producing the string "null/<name>".
    //
    public static String dbStoreDir(params, workflow, String kind, db_name) {
        def base = dbCacheDir(params, workflow, kind)
        if (!base) { return null }
        return "${base}/${db_name.toString().replaceAll('[^A-Za-z0-9._-]', '_')}".toString()
    }

    //
    // Human-readable description of where a db download will land, for the
    // up-front console messages.
    //
    public static String dbCacheDescription(params, workflow, String kind) {
        def base = dbCacheDir(params, workflow, kind)
        return base ? "cached at ${base}" : 'staged in the work directory (set --db_cache_dir to cache it across runs)'
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

    //
    // ── Initialisation helpers ────────────────────────────────────────────────
    //
    // These used to live inline at the top of `workflow TAXTRIAGE`. They were moved
    // here because Nextflow captures each workflow body's source verbatim as a
    // string constant in the compiled class, and the JVM class-file format caps a
    // string at 65,535 bytes -- the body had grown past that and the pipeline no
    // longer compiled ("Module compilation error ... String too long"). Keeping
    // pure param validation/resolution in this class keeps the workflow body well
    // under the cap. Behaviour is unchanged.
    //

    // Flags that may arrive as strings on NF v26 and must be coerced before
    // schema validation.
    private static final List<String> BOOLEAN_PARAMS = [
        'annotate', 'centrifuge', 'trim', 'downsample', 'low_memory',
        'download_taxdump', 'download_db', 'add_irregular_top_hits',
        'save_output_fastqs', 'save_unaligned', 'remove_commensal',
        'save_k2_read_assignment', 'save_classified_fastq',
        'include_singletons_hostremoval', 'include_singletons_removal',
        'split_prefix', 'use_megahit_longreads', 'use_bt2', 'use_hisat2',
        'use_diamond', 'use_denovo', 'skip_report', 'skip_consensus',
        'skip_variants', 'skip_realignment', 'skip_confidence',
        'enable_genbank', 'get_pathogens', 'conf_sens', 'disable_auto_weights',
        'auto_score_power', 'fuzzy', 'refresh_download', 'igenomes_ignore',
        'recursive_reference', 'decompress_pre_megahit', 'skip_plots',
        'skip_stats', 'skip_fastp', 'skip_kraken2', 'skip_refpull',
        'skip_krona', 'skip_features', 'skip_pathogens', 'unknown_sample',
        'ignore_missing', 'reference_assembly', 'pathogenicity', 'get_features',
        'get_variants', 'compress_species', 'fast', 'enable_matrix', 'no_subkey',
        'sort_alphabetical', 'show_potentials', 'show_opportunistics',
        'show_commensals', 'show_unidentified', 'integrate_strain_table',
        'skip_multiqc', 'email_on_fail', 'plaintext_email', 'monochrome_logs',
        'help', 'validate_params', 'show_hidden_params', 'enable_conda'
    ]

    public static void coerceBooleanParams(params) {
        BOOLEAN_PARAMS.each { p ->
            if (params[p] instanceof String) { params[p] = params[p].toBoolean() }
        }
    }

    //
    // Container engine / classifier / input pre-flight checks.
    //
    public static void validateInputs(params, workflow, log) {
        if (workflow.containerEngine != 'singularity' && workflow.containerEngine != 'docker') {
            Nextflow.error("Neither Docker or Singularity was selected as the container engine. Please specify with `-profile docker` or `-profile singularity`. Exiting...")
        }

        if (params.classifier != 'kraken2' && params.classifier != 'centrifuge' && params.classifier != 'metaphlan') {
            Nextflow.error("Classifier must be either kraken2, centrifuge or metaphlan")
        }

        println "Working Directory: ${workflow.workDir}"

        if (params.bam) {
            if (!Nextflow.file(params.bam).exists()) {
                Nextflow.error("ERROR: bam file does not exist: ${params.bam}")
            }
        } else if (params.fastq_1) {
            // An SRA/ENA accession is not a path -- nothing exists on disk yet, it is
            // downloaded by INPUT_CHECK. Skip the existence pre-flight for those.
            if (isAccession(params.fastq_1)) {
                println "Detected SRA/ENA accession for --fastq_1: ${params.fastq_1} (reads will be downloaded)"
                if (params.fastq_2) {
                    log.warn "--fastq_2 is ignored when --fastq_1 is an accession; paired-end layout is detected from the archive."
                }
            } else {
                if (!Nextflow.file(params.fastq_1).exists()) {
                    Nextflow.error("ERROR: fastq_1 file does not exist: ${params.fastq_1}")
                }
                if (params.fastq_2 && !Nextflow.file(params.fastq_2).exists()) {
                    Nextflow.error("ERROR: fastq_2 file does not exist: ${params.fastq_2}")
                }
            }
        } else if (!params.input) {
            Nextflow.error('ERROR: Please specify an input samplesheet (--input), a fastq_1 file (--fastq_1) or an alignment (--bam)!')
        }
    }

    //
    // PRE-ALIGNED (BAM) INPUT
    // Peek at the samplesheet up front so a BAM-only run can relax the checks that
    // exist purely for the read-based path (Kraken2 DB, QC, reference download) and
    // warn about flags that cannot apply without raw reads.
    //
    public static Map scanBamSamplesheet(infile) {
        def scan = [any: false, all: false]
        if (!infile) { return scan }
        try {
            def f = Nextflow.file(infile)
            if (!f.exists()) { return scan }
            def lines = f.readLines().findAll { it != null && it.trim() }
            if (lines.size() < 2) { return scan }
            def sep = lines[0].contains('\t') ? '\t' : ','
            // strip a UTF-8 BOM if the samplesheet was saved from Excel
            def hdr = lines[0].split(sep, -1).collect { it.trim().replaceAll('\\uFEFF', '') }
            def bidx = hdr.indexOf('bam')
            if (bidx < 0) { return scan }
            def flags = lines[1..-1].collect { l ->
                def cols = l.split(sep, -1)
                (bidx < cols.size() && cols[bidx] != null && cols[bidx].trim()) ? true : false
            }
            scan.any = flags.any { it }
            scan.all = flags.every { it }
        /* groovylint-disable-next-line CatchException */
        } catch (Exception e) {
            println "WARNING: could not pre-scan ${infile} for a bam column: ${e.message}"
        }
        return scan
    }

    //
    // match_paths.py needs the reference(s) the BAM was aligned against: the FASTA
    // for sourmash / shared-window / ANI comparison and the derived accession->taxid
    // map for -m.  A BAM header carries reference NAMES and LENGTHS but no bases, so
    // when no reference is supplied we reconstruct one by calling consensus off the
    // alignment itself (BAM_CONSENSUS).  --bam_consensus false opts out and instead
    // runs without the minhash / conflict component.
    //
    public static boolean resolveBamConsensusMode(params, boolean has_bam_samples) {
        if (!has_bam_samples) { return false }
        def consensus_opt_out = (params.bam_consensus != null && !params.bam_consensus)
        if (!params.reference_fasta && !params.get_pathogens) {
            if (consensus_opt_out) {
                println 'WARNING: pre-aligned input without --reference_fasta and --bam_consensus false: ' +
                        'no reference sequence is available, so sourmash/ANI comparison and conflict-based ' +
                        'read removal are disabled. Set --minhash_weight 0 to rebalance the TASS weights.'
                return false
            }
            println 'NOTE: pre-aligned input without --reference_fasta -> reconstructing reference ' +
                    'sequence from the alignment (samtools consensus).'
            println 'WARNING: consensus-derived references only cover positions with aligned reads, and ' +
                    'multi-mapping reads contribute to every reference they were placed on, which ' +
                    'overstates similarity between related organisms and makes conflict-driven read ' +
                    'removal more aggressive. Pass --reference_fasta whenever the true reference is available.'
            return true
        }
        if (params.bam_consensus) {
            println 'NOTE: --bam_consensus set explicitly; consensus sequence will be derived from the ' +
                    'alignment in addition to the supplied reference.'
            return true
        }
        return false
    }

    //
    // Flags that need raw reads / de novo contigs are turned off for a BAM-only run.
    //
    public static void applyBamOnlyOverrides(params, boolean bam_only_run, boolean has_bam_samples) {
        if (bam_only_run) {
            println 'BAM-only run detected: skipping read QC, trimming, host removal, classification and reference download.'
            // Nothing to classify and nothing to select references from.
            params.skip_kraken2 = true
            params.skip_refpull = true
            [
                'use_denovo', 'use_diamond', 'annotate', 'microbert', 'novelty',
                'generate_iss', 'generate_nanosim', 'reference_assembly', 'get_variants'
            ].each { flag ->
                if (params[flag]) {
                    println "WARNING: --${flag} is not supported for pre-aligned (BAM) input and has been disabled."
                    params[flag] = false
                }
            }
        } else if (has_bam_samples) {
            println 'Mixed FASTQ/BAM samplesheet detected: pre-aligned samples bypass QC, ' +
                    'classification, reference download, de novo assembly, MicrobeRT and novelty.'
        }
    }

    //
    // Database requirements. Pre-aligned samples are exempt: the references are
    // already fixed by the BAM, and their sequence comes either from
    // --reference_fasta or from BAM_CONSENSUS.
    //
    public static void requireDatabases(params, boolean has_bam_samples) {
        if (!params.skip_kraken2 && !params.db && !params.download_db) {
            Nextflow.error("If --skip_kraken2 is false, you must provide --db or --download_db")
        }
        if (params.skip_kraken2 && !has_bam_samples && !params.reference_fasta && !params.get_pathogens && !params.organisms && !params.organisms_file) {
            Nextflow.error("If you are skipping kraken2, you must provide a reference fasta, --get_pathogens to pull the pathogens file, organisms, or organisms_file")
        }
    }

    //
    // --pathogens override -> the sheet to read, else the bundled default.
    //
    public static String resolvePathogensSheet(params, projectDir) {
        if (!params.pathogens) {
            return "${projectDir}/assets/pathogen_sheet.csv".toString()
        }
        if (!(params.pathogens.endsWith('.csv') || params.pathogens.endsWith('.txt'))) {
            Nextflow.error("Pathogens file must end with .csv or .txt i.e. it is a .csv (comma-delimited) file!")
        }
        return params.pathogens.toString()
    }

    //
    // --assembly (+ optional --assembly_summary_genbank) -> a Path, a List of Paths,
    // or null when the summaries should be downloaded instead.
    //
    // Local RefSeq summary: --assembly or its alias --assembly_summary_refseq.
    public static String localRefseqSummary(params) {
        return params.assembly ?: params.assembly_summary_refseq ?: null
    }

    // Returns the local summary file(s) when NOTHING needs downloading, else null
    // (GET_ASSEMBLIES then downloads whatever is missing and the workflow mixes
    // in any local files). RefSeq is always element [0], GenBank [1].
    public static Object resolveAssemblyFiles(params) {
        def refseq  = localRefseqSummary(params)
        def genbank = params.assembly_summary_genbank
        def need_genbank_download = params.enable_genbank && !genbank
        if (!refseq || need_genbank_download) {
            println 'Assembly summaries: ' +
                (refseq  ? "RefSeq local (${refseq})" : 'RefSeq download') + ', ' +
                (genbank ? "GenBank local (${genbank})" : (params.enable_genbank ? 'GenBank download' : 'GenBank disabled (enable with --enable_genbank or --assembly_summary_genbank)'))
            return null
        }
        println "Assembly file present, using it to pull genomes from... ${refseq}"
        def files = [Nextflow.file(refseq, checkIfExists: true)]
        if (genbank) {
            println "GenBank assembly file also provided: ${genbank}"
            files << Nextflow.file(genbank, checkIfExists: true)
        }
        return files.size() == 1 ? files[0] : files
    }
}
