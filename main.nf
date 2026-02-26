nextflow.enable.dsl = 2

// ──────────────────────────────────────────────────────────────────────────────
// OpenCare banner (visual only; toggle with --no_logo)
// ──────────────────────────────────────────────────────────────────────────────
import groovy.transform.Field

@Field final String BGWHITE = '\u001B[47m'
@Field final String CYAN    = '\u001B[36m'
@Field final String GREEN   = '\u001B[32m'
@Field final String RESET   = '\u001B[0m'

String makeBox(String text, int pad = 4) {
    def inner = (' ' * pad) + text + (' ' * pad)
    def width = inner.size()
    def top   = "╔" + ("═" * width) + "╗"
    def mid   = "║" + inner + "║"
    def bot   = "╚" + ("═" * width) + "╝"

    int cols   = (System.getenv('COLUMNS') ?: '80') as int
    int boxLen = width + 2
    int left   = Math.max(0, ((cols - boxLen) / 2) as int)
    def M      = ' ' * left

    return [
        "${M}${BGWHITE}${CYAN}${top}${RESET}",
        "${M}${BGWHITE}${CYAN}${mid}${RESET}",
        "${M}${BGWHITE}${CYAN}${bot}${RESET}"
    ].join('\n')
}

def OPENCARE_BANNER = """
${makeBox('OpenCare Genomics Pipeline — v1.0.0')}

${GREEN}✔ OpenCare validation ready${RESET}
"""

// ──────────────────────────────────────────────────────────────────────────────
// params & defaults (only those used by main.nf; others live in sub-workflows)
// ──────────────────────────────────────────────────────────────────────────────
params.outdir         = params.outdir         ?: 'results'
params.patient_id     = params.patient_id     ?: 'PATIENT01'
params.pathway_db     = params.pathway_db     ?: null
params.gene_domains   = params.gene_domains   ?: null
params.mcode          = params.mcode          ?: null
params.arm_script        = params.arm_script        ?: "${projectDir}/scripts/call_arm_events.py"
params.bench_metrics_tsv = params.bench_metrics_tsv ?: "${params.outdir}/bench_metrics.tsv"

params.oncokb_token   = params.oncokb_token   ?: null
params.ncbi_key       = params.ncbi_key       ?: null
params.ncbi_email     = params.ncbi_email     ?: null
params.enable_online  = params.enable_online  ?: null

params.reads          = params.reads          ?: null
params.ref_fa         = params.ref_fa         ?: null
params.targets_bed    = params.targets_bed    ?: null
params.vep_cache      = params.vep_cache      ?: null
params.knowledge_tsv  = params.knowledge_tsv  ?: null

params.report_builder = params.report_builder ?: 'html'
params.build_html     = params.build_html     ?: "${baseDir}/HTML/build_html.py"

// ──────────────────────────────────────────────────────────────────────────────
// lightweight validation (avoid Nextflow file() to keep this minimal)
// ──────────────────────────────────────────────────────────────────────────────
if( !params.reads )  log.warn 'Provide --reads'
if( !params.ref_fa ) log.warn 'Provide --ref_fa'

if( !params.containsKey('no_logo') )              params.no_logo = false
if( !params.containsKey('region') )               params.region  = null
if( !params.containsKey('discard_trimmed_pass') ) params.discard_trimmed_pass = true

// Initial placeholder creation as File objects
if( !params.containsKey('empty_json') ) {
  def p = file("${projectDir}/.empty.json")
  if( !p.exists() ) p.text = ''
  params.empty_json = p
}
if( !params.containsKey('empty_tsv') ) {
  def p = file("${projectDir}/.empty.tsv")
  if( !p.exists() ) p.text = ''
  params.empty_tsv = p
}

// Optional existence note for the reference
if( params.ref_fa && !(new File(params.ref_fa as String).exists()) ) {
  log.warn "ref_fa not found: ${params.ref_fa}"
}
else if( params.ref_fa ) {
  log.info "ref_fa found"
}

// Reaffirm placeholders as path strings (kept as-is to preserve behavior)
if( !params.empty_json ) {
  def p = new File("${projectDir}/.empty.json")
  if( !p.exists() ) p.text = ''     // create a tiny placeholder
  params.empty_json = p.path        // pass the path string
}
if( !params.empty_tsv ) {
  def p = new File("${projectDir}/.empty.tsv")
  if( !p.exists() ) p.text = ''
  params.empty_tsv = p.path
}

// Optional banner
if( !params.no_logo ) println OPENCARE_BANNER

// ──────────────────────────────────────────────────────────────────────────────
// include sub-workflows BEFORE calling them
// ──────────────────────────────────────────────────────────────────────────────
include { CLINICAL_CORE as ClinicalCore } from './workflows/clinical_core.nf'

// ──────────────────────────────────────────────────────────────────────────────
// top-level workflow orchestration
// ──────────────────────────────────────────────────────────────────────────────
workflow {
  // (1) Build (sid, r1, r2)
  def reads_ch = Channel
    .fromFilePairs(params.reads, checkIfExists: true)
    .map { sid, files ->
      def list = (files instanceof List) ? files : [files]
      assert list.size() in [1,2] : "Expected 1 or 2 files for sample ${sid}, got ${list}"
      def r1 = list[0]
      def r2 = (list.size() > 1) ? list[1] : list[0]
      tuple(sid.toString(), r1, r2)
    }

  // (2) Reference as a singleton channel
  def ref_fa_ch = Channel.value(params.ref_fa as String)

  // (3) Call sub-workflow
  def core = ClinicalCore(reads_ch, ref_fa_ch)

  // (4) Re-expose outputs
  emit:
    bam_ch      = core.bam_ch
    vcf_ch      = core.vcf_ch
    panel_ch    = core.panel_ch
    report_json = core.report_json
    mcode_json  = core.mcode_json
    report_html = core.report_html
}
