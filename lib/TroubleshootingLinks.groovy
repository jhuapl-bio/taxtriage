// ##############################################################################################
// # Copyright 2022 The Johns Hopkins University Applied Physics Laboratory LLC
// # All rights reserved.
// # Permission is hereby granted, free of charge, to any person obtaining a copy of this
// # software and associated documentation files (the "Software"), to deal in the Software
// # without restriction, including without limitation the rights to use, copy, modify,
// # merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
// # permit persons to whom the Software is furnished to do so.
// #
// # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
// # INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
// # PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// # LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
// # TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
// # OR OTHER DEALINGS IN THE SOFTWARE.
// #

//
// End-of-run troubleshooting links.
//
// Called from workflow.onComplete. Reads this run's execution trace, finds every task
// whose final attempt FAILED (both the one that stopped the run and the many optional
// steps that run with errorStrategy 'ignore'), reads the tail of each task's
// .command.err, and maps it to the matching section of docs/troubleshooting.md:
//
//   * printed to the terminal (plain https URLs: Cmd/Ctrl-click in most terminals)
//   * pipeline_info/troubleshooting_<timestamp>.html  (clickable, shareable)
//   * pipeline_info/troubleshooting_<timestamp>.tsv   (machine readable)
//
// Rule order matters: the error text is the most specific signal, then the exit code,
// then the process name. Anchors are the GitHub/MkDocs slugs of the doc headings
// (kept in sync by tests: every anchor below must exist in docs/troubleshooting.md).
//

import nextflow.Nextflow
import java.nio.file.Files
import java.nio.file.Path

class TroubleshootingLinks {

    static final String DEFAULT_URL  = 'https://jhuapl-bio.github.io/taxtriage/stable/troubleshooting/'
    static final int    MAX_DETAILED = 50     // tasks whose .command.err is read / listed in full
    static final int    TAIL_BYTES   = 65536  // only the end of each .command.err is read

    static final String A_COMPRESSED  = 'compressed-gz-reference-fastas'
    static final String A_VARIANTS    = 'variant-calling-and-reference-assembly-bcftools--bedtools_maskfasta'
    static final String A_HOST        = 'host-removal-minimap2--deacon'
    static final String A_ALIGN       = 'bowtie2--minimap2-alignment'
    static final String A_OOM         = 'out-of-memory-or-killed-exit-137--143'
    static final String A_REQ         = 'process-requirement-exceeds-available-cpus-or-memory'
    static final String A_CMD         = 'command-not-found-exit-127-or-container-errors'
    static final String A_DISK        = 'no-space-left-on-device'
    static final String A_NETWORK     = 'assemblies-not-downloading'
    static final String A_SHEET       = 'my-pipeline-fails-at-samplesheet_check'
    static final String A_UNCOMMITTED = 'uncommitted-changes-warning'
    static final String A_KRAKEN      = 'kraken2_kraken2'
    static final String A_FASTP       = 'fastp'
    static final String A_HIST        = 'samtools_hist_coverage'
    static final String A_ASSEMBLY    = 'spades--flye-de-novo-assembly'
    static final String A_SORT        = 'samtools_sort'
    static final String A_VIEW        = 'samtools_view'
    static final String A_NANOPLOT    = 'nanoplot-failing-on-empty-data'
    static final String A_REPORT      = 'a-sample-is-missing-from-the-interactive-report'
    static final String A_IGNORED     = 'steps-that-failed-but-the-run-still-finished'
    static final String A_DEBUG       = 'debugging-a-specific-step'
    static final String A_HELP        = 'getting-help'

    // [regex on the error text, anchor, short title]
    static final List TEXT_RULES = [
        [/(?i)cannot index files compressed with gzip|not in gzip format|error when opening file `?\S+\.b?gz|_M_range_check/, A_COMPRESSED, 'Compressed reference FASTA'],
        [/(?i)cannot open '\S+ (bowtie2|hisat2)'/,                                          A_VARIANTS,    'Variant calling given the index directory'],
        [/(?i)no space left on device|disk quota exceeded/,                                    A_DISK,        'Disk full'],
        [/(?i)process requirement exceeds available (cpus|memory)/,                            A_REQ,         'Asked for more CPUs/memory than available'],
        [/(?i)command not found|cannot connect to the docker daemon|docker: not found|unable to find image|manifest unknown|failed to pull image|singularity: not found|apptainer: not found/, A_CMD, 'Tool or container missing'],
        [/(?i)out of memory|outofmemoryerror|std::bad_alloc|cannot allocate memory|memory exhausted|\bkilled\b|oom-kill/, A_OOM, 'Out of memory / killed'],
        [/(?i)could not resolve host|name or service not known|connection (timed out|refused|reset)|curl: \((6|7|28|35|56|60)\)|ssl: certificate|certificate verify failed|urlopen error|http error (4|5)\d\d|temporary failure in name resolution/, A_NETWORK, 'Network / download failure'],
        [/(?i)\[critical\].*(fastq|samplesheet|line \d)/,                                     A_SHEET,       'Samplesheet problem'],
        [/(?i)uncommitted changes/,                                                            A_UNCOMMITTED, 'Uncommitted changes in the cached pipeline'],
        [/(?i)kraken2 db could not be loaded|\bhash\.k2d\b|\btaxo\.k2d\b|\bopts\.k2d\b/,       A_KRAKEN,      'Kraken2 database'],
        [/(?i)expected -c or -m option/,                                                       A_VARIANTS,    'bcftools call options missing'],
    ]

    // Only for failures outside any task (pipeline setup): missing inputs, bad params.
    static final List SETUP_RULES = [
        [/(?i)does not exist|no such file|cannot find|can't find|not a valid path|unable to read|not found:/, 'input-issues', 'Input file / parameter problem'],
        [/(?i)does not have a fasta extension|accepted: \.fa/,                                           A_COMPRESSED, 'Reference FASTA name not accepted'],
    ]

    // exit code -> [anchor, title]
    static final Map EXIT_RULES = [
        137: [A_OOM, 'Out of memory / killed'],
        143: [A_OOM, 'Killed by the scheduler (memory or time)'],
        127: [A_CMD, 'Command not found'],
    ]

    // [regex on the process name (last path component), anchor, short title]
    static final List PROCESS_RULES = [
        [/SAMPLESHEET_CHECK|INPUT_CHECK/,                                                   A_SHEET,    'Samplesheet check'],
        [/DEACON|FILTER_MINIMAP2|REMOVE_HOSTREADS|FETCH_HOST_REFS|CHECK_GZIPPED_READS|FILTERED_STATS/, A_HOST, 'Host removal'],
        [/FASTP/,                                                                           A_FASTP,    'fastp'],
        [/KRAKEN2/,                                                                         A_KRAKEN,   'Kraken2'],
        [/BCFTOOLS|BEDTOOLS_MASKFASTA|BEDTOOLS_BAMTOBED|SPLIT_VCF/,                          A_VARIANTS, 'Variant calling / reference assembly'],
        [/DECOMPRESS_|PULL_FASTA|MAP_LOCAL_ASSEMBLY|MAP_SAMPLE_FASTA/,                       A_COMPRESSED, 'Reference FASTA handling'],
        [/MINIMAP2|BOWTIE2|HISAT2|BWA_/,                                                     A_ALIGN,    'Alignment'],
        [/DOWNLOAD|FEATURES_DOWNLOAD|GET_ASSEMBLY|FCS_DOWNLOAD|MAP_TAXID_ASSEMBLY|REMOTE_BLAST|FETCH_/, A_NETWORK, 'Download from NCBI / remote'],
        [/SAMTOOLS_HIST_COVERAGE/,                                                          A_HIST,     'Coverage histogram'],
        [/SPADES|FLYE|MEGAHIT/,                                                             A_ASSEMBLY, 'De novo assembly'],
        [/SAMTOOLS_SORT/,                                                                   A_SORT,     'samtools sort'],
        [/SAMTOOLS_VIEW/,                                                                   A_VIEW,     'samtools view'],
        [/NANOPLOT/,                                                                        A_NANOPLOT, 'NanoPlot'],
        [/ALIGNMENT_PER_SAMPLE|ORGANISM_MERGE_REPORT|REPORT|CONFIDENCE|COMBINE_SAMPLES_JSON/, A_REPORT,  'Report / scoring'],
    ]

    private static boolean reported = false

    // ── public entry points ──────────────────────────────────────────────────────

    /**
     * workflow.onError hook. Errors raised while the workflow is being built (bad
     * parameter, missing input file, channel/Groovy error) fire onError but NOT
     * onComplete, and have no completion time yet -- report those here. Task
     * failures fire both; they are reported from onComplete, once the trace is final.
     */
    public static void onError(workflow, params, log) {
        if (workflow.complete == null) report(workflow, params, log)
    }

    /**
     * Called at the top of the workflow body (kept to one line there: the TAXTRIAGE
     * body is close to the JVM's 64 KB method-constant limit).
     */
    public static void registerOnError(workflow, params, log) {
        def wf = workflow, p = params, lg = log
        workflow.onError { onError(wf, p, lg) }
    }

    /** workflow.onComplete hook (also safe to call twice: prints once per run). */
    public static void report(workflow, params, log) {
        if (reported) return
        reported = true
        try {
            def base = baseUrl(params)
            if (base == null) return

            def failures = failedTasks(workflow, params)
            boolean stopped = !workflow.success
            if (!failures && !stopped) return

            // Non-task failure (bad parameter, missing input, Groovy/channel error ...):
            // nothing FAILED in the trace, but the run stopped. Classify its message.
            if (!failures && stopped) {
                def msg = (workflow.errorMessage ?: workflow.errorReport ?: '').toString()
                def first = msg.readLines().find { it.trim() } ?: 'see .nextflow.log'
                def rule = classify(null, null, msg)
                if (rule[0] == A_DEBUG) {
                    def hit = SETUP_RULES.find { msg =~ (it[0] as String) }
                    if (hit) rule = [hit[1], hit[2]]
                }
                failures = [[name: 'Pipeline setup', process: 'Pipeline setup', tag: '', exit: '-', workdir: '', ignored: false,
                             error: first.trim().take(240), anchor: rule[0], title: rule[1], fatal: true]]
            }

            def colors = NfcoreTemplate.logColours(params.monochrome_logs)
            def outFiles = writeFiles(workflow, params, failures, base)
            printBlock(workflow, log, colors, failures, base, outFiles)
        } catch (Throwable e) {
            // Never let the helper mask the real outcome of the run.
            log.warn "Could not build troubleshooting links: ${e.message}"
        }
    }

    // ── helpers ──────────────────────────────────────────────────────────────────

    static String baseUrl(params) {
        def v = params.containsKey('troubleshooting_url') ? params.troubleshooting_url : DEFAULT_URL
        if (v == null || v == true) return DEFAULT_URL
        def s = v.toString().trim()
        if (v == false || s.toLowerCase() in ['false', 'off', 'none', 'no', '0', '']) return null
        return s.contains('#') ? s.substring(0, s.indexOf('#')) : s
    }

    static String link(String base, String anchor) {
        return "${base}#${anchor}"
    }

    /** First matching rule: error text, then exit code, then process name, then generic. */
    static List classify(String process, Object exit, String errText) {
        def text = errText ?: ''
        for (r in TEXT_RULES) {
            if (text =~ (r[0] as String)) return [r[1], r[2], r[0]]
        }
        def code = exit?.toString()?.isInteger() ? exit.toString().toInteger() : null
        if (code != null && EXIT_RULES.containsKey(code)) return EXIT_RULES[code]
        def simple = (process ?: '').tokenize(':') ? (process ?: '').tokenize(':')[-1] : ''
        for (r in PROCESS_RULES) {
            if (simple =~ (r[0] as String)) return [r[1], r[2]]
        }
        return [A_DEBUG, 'Inspect the failed step']
    }

    /** Path of this run's trace file (as configured), or null. */
    static Object traceFile(workflow, params) {
        try {
            def cfg = nextflow.Global.session?.config
            def tr = cfg?.navigate('trace')
            if (tr instanceof Map && tr.enabled == false) return null
            def f = (tr instanceof Map) ? tr.file : null
            if (f) return Nextflow.file(f.toString())
        } catch (Throwable ignored) { }
        if (params.outdir && params.trace_timestamp) {
            return Nextflow.file("${params.outdir}/pipeline_info/execution_trace_${params.trace_timestamp}.txt")
        }
        return null
    }

    /**
     * Tasks whose FINAL attempt failed. A task retried successfully appears first as
     * FAILED and later as COMPLETED under the same name; only the last row counts.
     */
    static List<Map> failedTasks(workflow, params) {
        def tf = traceFile(workflow, params)
        if (tf == null || !Files.exists(tf as Path)) return []
        def lines = (tf as Path).readLines()
        if (lines.size() < 2) return []
        def header = lines[0].split('\t') as List
        int iName = header.indexOf('name'), iStatus = header.indexOf('status')
        if (iName < 0 || iStatus < 0) return []
        int iExit = header.indexOf('exit'), iHash = header.indexOf('hash'), iWork = header.indexOf('workdir')

        Map<String, List> last = new LinkedHashMap<>()
        lines.drop(1).each { ln ->
            def c = ln.split('\t', -1)
            if (c.size() <= Math.max(iName, iStatus)) return
            last.remove(c[iName])            // re-insert so order follows the final attempt
            last[c[iName]] = c as List
        }
        def fatalProcess = fatalProcessName(workflow)
        def out = []
        last.each { name, c ->
            if (c[iStatus] != 'FAILED') return
            def m = (name =~ /^(.*?)(?: \((.*)\))?$/)
            def proc = m.matches() ? m.group(1) : name
            def tag  = m.matches() ? (m.group(2) ?: '') : ''
            def work = iWork >= 0 ? c[iWork] : workdirFromHash(workflow, iHash >= 0 ? c[iHash] : null)
            out << [name: name, process: proc, tag: tag,
                    exit: iExit >= 0 ? c[iExit] : '-', workdir: work ?: '',
                    fatal: !workflow.success && fatalProcess && (name == fatalProcess || name.endsWith(':' + fatalProcess)),
                    ignored: workflow.success as boolean]
        }
        // Read the error tail and classify (bounded, so a run with hundreds of ignored
        // tasks does not stall the exit; the rest still get a process-based link).
        out.eachWithIndex { t, i ->
            def tail = i < MAX_DETAILED ? errorTail(t.workdir as String) : ''
            // Remote/unreadable work dir: fall back to the 'Command error:' block of
            // Nextflow's own report for the task that stopped the run.
            if (!tail.trim() && t.fatal) tail = commandErrorBlock(workflow.errorReport) ?: (workflow.errorReport ?: '').toString()
            def rule = classify(t.process as String, t.exit, tail)
            t.anchor = rule[0]; t.title = rule[1]
            t.error = headline(tail, rule.size() > 2 ? rule[2] as String : null)
        }
        // The task that stopped the run first, then by process name.
        return out.sort { a, b -> (b.fatal <=> a.fatal) ?: (a.name <=> b.name) }
    }

    /** Lines between 'Command error:' and 'Work dir:' in Nextflow's error report. */
    static String commandErrorBlock(Object report) {
        def lines = (report ?: '').toString().readLines()
        int a = lines.findIndexOf { it.trim() == 'Command error:' }
        if (a < 0) return ''
        int b = lines.findIndexOf(a + 1) { it.trim() ==~ /(Work dir|Container|Tip):.*/ }
        return lines.subList(a + 1, b < 0 ? lines.size() : b).join('\n')
    }

    static String fatalProcessName(workflow) {
        def rep = (workflow.errorReport ?: '').toString()
        def m = rep =~ /Error executing process > '([^']+)'/
        if (!m.find()) return null
        return m.group(1)
    }

    static String workdirFromHash(workflow, String hash) {
        if (!hash || !hash.contains('/')) return ''
        try {
            def parts = hash.split('/')
            def dir = (workflow.workDir as Path).resolve(parts[0])
            def hit = Files.list(dir).withCloseable { s -> s.filter { it.fileName.toString().startsWith(parts[1]) }.findFirst() }
            return hit.isPresent() ? hit.get().toString() : ''
        } catch (Throwable ignored) { return '' }
    }

    /** Last ~64 KB of .command.err (or .command.log when stderr is empty). */
    static String errorTail(String workdir) {
        if (!workdir) return ''
        for (name in ['.command.err', '.command.log']) {
            try {
                def p = Nextflow.file("${workdir}/${name}") as Path
                if (!Files.exists(p)) continue
                long size = Files.size(p)
                if (size == 0) continue
                def bytes
                if (size <= TAIL_BYTES) {
                    bytes = Files.readAllBytes(p)
                } else {
                    def ch = Files.newByteChannel(p)
                    try {
                        ch.position(size - TAIL_BYTES)
                        def buf = java.nio.ByteBuffer.allocate(TAIL_BYTES)
                        while (buf.hasRemaining() && ch.read(buf) > 0) { }
                        bytes = buf.array()
                    } finally { ch.close() }
                }
                def text = new String(bytes, 'UTF-8')
                if (text.trim()) return text
            } catch (Throwable ignored) { }
        }
        return ''
    }

    /** The single most useful line of an error tail. */
    static String headline(String text, String matched = null) {
        def lines = (text ?: '').readLines().collect { it.replaceAll(/\u001B\[[;\d]*m/, '').trim() }.findAll { it }
        if (!lines) return ''
        def noise = /(?i)^(note:|warn(ing)?:? .*deprecat|\[mpileup\]|picked up java_tool_options|-\[|executor >)/
        def useful = lines.findAll { !(it =~ noise) }
        // The line that decided the help topic explains the failure best.
        def pick = (matched ? useful.reverse().find { it =~ matched } : null) ?: useful.reverse().find { it =~ /(?i)error|exception|fail|cannot|can't|unable|not found|denied|killed|invalid|no such|exceed|abort|terminate|space left/ } ?:
                   (useful ? useful[-1] : lines[-1])
        return pick.length() > 240 ? pick.take(237) + '...' : pick
    }

    static Map writeFiles(workflow, params, List<Map> failures, String base) {
        def res = [:]
        if (!params.outdir) return res
        def ts = params.trace_timestamp ?: new Date().format('yyyy-MM-dd_HH-mm-ss')
        try {
            def dir = Nextflow.file("${params.outdir}/pipeline_info") as Path
            if (!Files.exists(dir)) Files.createDirectories(dir)
            def html = dir.resolve("troubleshooting_${ts}.html")
            def tsv  = dir.resolve("troubleshooting_${ts}.tsv")
            html.text = renderHtml(workflow, failures, base)
            tsv.text  = renderTsv(failures, base)
            res.html = html; res.tsv = tsv
        } catch (Throwable e) {
            res.error = e.message
        }
        return res
    }

    static String renderTsv(List<Map> failures, String base) {
        def rows = ['process\ttag\texit\tstatus\terror\tworkdir\thelp_title\thelp_url']
        failures.each { t ->
            rows << [t.process, t.tag, t.exit, t.fatal ? 'stopped_run' : (t.ignored ? 'ignored' : 'failed'),
                     (t.error ?: '').replace('\t', ' '), t.workdir, t.title, link(base, t.anchor as String)].join('\t')
        }
        return rows.join('\n') + '\n'
    }

    static String esc(Object s) {
        return (s ?: '').toString().replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;').replace('"', '&quot;')
    }

    static String finishedAt(workflow) {
        try {
            if (workflow.complete) return workflow.complete.format(java.time.format.DateTimeFormatter.ofPattern('yyyy-MM-dd HH:mm:ss'))
        } catch (Throwable ignored) { }
        return new Date().format('yyyy-MM-dd HH:mm:ss')
    }

    static String workLink(String wd) {
        if (!wd) return ''
        // Local work dirs become file:// links (opens the folder); remote ones are shown as text.
        if (wd.startsWith('/')) return "<a href=\"file://${esc(wd)}\"><code>${esc(wd)}</code></a>"
        return "<code>${esc(wd)}</code>"
    }

    static String renderHtml(workflow, List<Map> failures, String base) {
        def status = workflow.success ? (failures ? 'completed, with failed optional steps' : 'completed') : 'stopped with an error'
        def byAnchor = failures.groupBy { it.anchor }
        def rows = failures.collect { t ->
            """<tr class="${t.fatal ? 'fatal' : ''}">
  <td><b>${esc(t.process.toString().tokenize(':')[-1])}</b><div class="dim">${esc(t.process)}</div></td>
  <td>${esc(t.tag)}</td>
  <td class="num">${esc(t.exit)}</td>
  <td>${t.fatal ? '<span class="pill red">stopped run</span>' : (t.ignored ? '<span class="pill">ignored</span>' : '<span class="pill">failed</span>')}</td>
  <td><code class="err">${esc(t.error)}</code><div class="dim">${workLink(t.workdir as String)}</div></td>
  <td><a href="${esc(link(base, t.anchor as String))}" target="_blank" rel="noopener">${esc(t.title)} &#8599;</a></td>
</tr>"""
        }.join('\n')
        def topics = byAnchor.collect { a, ts ->
            "<li><a href=\"${esc(link(base, a as String))}\" target=\"_blank\" rel=\"noopener\">${esc(ts[0].title)}</a> <span class=\"dim\">(${ts.size()} task${ts.size() == 1 ? '' : 's'})</span></li>"
        }.join('\n')
        return """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>TaxTriage troubleshooting</title>
<style>
:root{--bg:#fff;--fg:#1d2330;--dim:#6b7280;--line:#e5e7eb;--accent:#2563eb;--red:#b91c1c;--redbg:#fef2f2;--pill:#f3f4f6}
@media (prefers-color-scheme: dark){:root{--bg:#111318;--fg:#e5e7eb;--dim:#9ca3af;--line:#2a2f3a;--accent:#60a5fa;--red:#f87171;--redbg:#2a1515;--pill:#1f2430}}
body{margin:0;background:var(--bg);color:var(--fg);font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif}
main{max-width:1200px;margin:0 auto;padding:24px 16px}
h1{font-size:20px;margin:0 0 4px} .dim{color:var(--dim);font-size:12px} a{color:var(--accent)}
table{width:100%;border-collapse:collapse;margin-top:16px} th,td{text-align:left;vertical-align:top;padding:8px;border-bottom:1px solid var(--line)}
th{font-size:12px;text-transform:uppercase;letter-spacing:.04em;color:var(--dim)} td.num{font-variant-numeric:tabular-nums}
tr.fatal{background:var(--redbg)} code{font:12px ui-monospace,SFMono-Regular,Menlo,monospace;overflow-wrap:anywhere} code.err{white-space:pre-wrap}
.pill{display:inline-block;padding:1px 8px;border-radius:10px;background:var(--pill);font-size:12px;white-space:nowrap} .pill.red{color:var(--red)}
.wrap{overflow-x:auto} ul{margin:8px 0 0;padding-left:20px}
@media (max-width:720px){table,tbody,tr,td{display:block;width:auto}thead{display:none}tr{padding:10px 0;border-bottom:1px solid var(--line)}td{border:0;padding:2px 0}td.num:before{content:"exit ";color:var(--dim)}}
</style></head><body><main>
<h1>TaxTriage run ${esc(workflow.runName)}: ${esc(status)}</h1>
<div class="dim">${esc(failures.size())} failed task(s) &middot; finished ${esc(finishedAt(workflow))} &middot; launch dir <code>${esc(workflow.launchDir)}</code></div>
<h2 style="font-size:15px;margin-top:20px">Where to start</h2>
<ul>
${topics}
<li><a href="${esc(link(base, A_IGNORED))}" target="_blank" rel="noopener">Steps that failed but the run still finished</a></li>
<li><a href="${esc(link(base, A_DEBUG))}" target="_blank" rel="noopener">Debugging a specific step</a> &middot; <a href="${esc(link(base, A_HELP))}" target="_blank" rel="noopener">Getting help</a></li>
</ul>
<div class="wrap"><table>
<thead><tr><th>Step</th><th>Sample / tag</th><th>Exit</th><th>Status</th><th>Last error line &middot; work dir</th><th>Help</th></tr></thead>
<tbody>
${rows}
</tbody></table></div>
</main></body></html>
"""
    }

    static void printBlock(workflow, log, Map colors, List<Map> failures, String base, Map outFiles) {
        def c = colors
        def sb = new StringBuilder()
        def setupOnly = failures.size() == 1 && failures[0].process == 'Pipeline setup'
        def what = setupOnly ? 'the run stopped before any task failed' : "${failures.size()} failed task(s)"
        sb << "-${c.purple}[${workflow.manifest.name}]${c.reset} ${c.bold}Troubleshooting:${c.reset} ${what}"
        if (workflow.success) sb << " ${c.dim}(optional steps; the run still completed)${c.reset}"
        sb << '\n'
        failures.take(MAX_DETAILED).each { t ->
            def label = t.fatal ? "${c.red}[stopped the run]${c.reset}" :
                        (workflow.success ? "${c.yellow}[ignored]${c.reset}" : "${c.yellow}[failed]${c.reset}")
            sb << "  ${c.red}✗${c.reset} ${t.name}  exit ${t.exit}  ${label}\n"
            if (t.error)   sb << "      ${c.dim}error:${c.reset} ${t.error}\n"
            if (t.workdir) sb << "      ${c.dim}work: ${c.reset} ${t.workdir}\n"
            sb << "      ${c.dim}help: ${c.reset} ${c.cyan}${link(base, t.anchor as String)}${c.reset}\n"
        }
        if (failures.size() > MAX_DETAILED) {
            sb << "  ... and ${failures.size() - MAX_DETAILED} more (see the files below)\n"
        }
        if (outFiles.html) {
            def p = outFiles.html as Path
            def shown = p.toUri().toString()
            sb << "  ${c.dim}Full list with links:${c.reset} ${shown}\n"
        }
        sb << "  ${c.dim}Debugging a step:${c.reset} ${link(base, A_DEBUG)}"
        log.info sb.toString()
    }
}
