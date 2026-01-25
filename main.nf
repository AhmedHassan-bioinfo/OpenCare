nextflow.enable.dsl = 2

// --- OpenCare banner ---
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
${makeBox('OpenCare Clinical Genomics Pipeline — v0.0.1')}

${GREEN}✔ OpenCare validation ready${RESET}
"""



/* -------- params & defaults -------- */
params.outdir        = params.outdir        ?: 'results'
params.patient_id    = params.patient_id    ?: 'PATIENT01'
params.pathway_db    = params.pathway_db    ?: null
params.gene_domains  = params.gene_domains  ?: null
params.mcode         = params.mcode         ?: null

params.oncokb_token  = params.oncokb_token  ?: null
params.ncbi_key      = params.ncbi_key      ?: null
params.ncbi_email    = params.ncbi_email    ?: null
params.enable_online = params.enable_online ?: null

params.reads         = params.reads         ?: null
params.ref_fa        = params.ref_fa        ?: null
params.targets_bed   = params.targets_bed   ?: null
params.vep_cache     = params.vep_cache     ?: null
params.knowledge_tsv = params.knowledge_tsv ?: null

params.report_builder = params.report_builder ?: 'html'
params.build_html     = params.build_html     ?: "${baseDir}/HTML/build_html.py"

/* -------- lightweight validation (no Nextflow file() here) -------- */
if (!params.reads)  log.warn 'Provide --reads'
if (!params.ref_fa) log.warn 'Provide --ref_fa'
if( !params.containsKey('no_logo') ) params.no_logo = false
if( !params.containsKey('region') )  params.region  = null
if( !params.containsKey('discard_trimmed_pass') ) params.discard_trimmed_pass = true
if( !params.containsKey('empty_json') ) {
  def p = file("${projectDir}/.empty.json"); if( !p.exists() ) p.text = ''
  params.empty_json = p
}
if( !params.containsKey('empty_tsv') ) {
  def p = file("${projectDir}/.empty.tsv");  if( !p.exists() ) p.text = ''
  params.empty_tsv = p
}
if (params.ref_fa && !(new File(params.ref_fa as String).exists())) {
  log.warn "ref_fa not found: ${params.ref_fa}"
} else if (params.ref_fa) {
  log.info "ref_fa found"
}
if( !params.empty_json ){
  def p = new File("${projectDir}/.empty.json")
  if( !p.exists() ) p.text = ''         // create a tiny placeholder
  params.empty_json = p.path            // pass the path string
}
if( !params.empty_tsv ){
  def p = new File("${projectDir}/.empty.tsv")
  if( !p.exists() ) p.text = ''
  params.empty_tsv = p.path
}
if( !params.no_logo ) println OPENCARE_BANNER

/* -------- include workflows BEFORE calling them -------- */
include { CLINICAL_CORE } from './workflows/clinical_core.nf'

/* -------- top-level workflow orchestration -------- */
workflow {
  // 1) Build (sid, r1, r2)
  Channel
    .fromFilePairs(params.reads, checkIfExists: true)
    .map { sid, files ->
      def list = (files instanceof List) ? files : [files]
      assert list.size() in [1,2] : "Expected 1 or 2 files for sample ${sid}, got ${list}"
      def r1 = list[0]
      def r2 = (list.size() > 1) ? list[1] : list[0]
      tuple(sid.toString(), r1, r2)
    }
    .set { reads_ch }

  // 2) Call sub-workflow
  def core = CLINICAL_CORE(reads_ch, params.ref_fa as String)

  // 3) Re-expose outputs
  emit:
    bam_ch      = core.bam_ch
    vcf_ch      = core.vcf_ch
    panel_ch    = core.panel_ch
    report_json = core.report_json
    mcode_json  = core.mcode_json
    report_html = core.report_html
}
