nextflow.enable.dsl = 2

/* ========================================================================== *
 * OpenCare Core
 * FASTQ → FastQC → fastp → BWA-MEM → sort/index → bcftools call
 * → (optional) VEP → (optional) PharmCAT → JSON/mCODE → Interactive HTML
 * ========================================================================== */
// direct include with a unique alias



//* ============================= PARAMETERS ================================= */
// ---- hard-required gnomAD resource (VCF + TBI) ----
if( !params.gnomad_vcf ) {
  error "Missing required parameter: --gnomad_vcf (bgzipped .vcf.gz path)"
}

def gnomad_vcf = file(params.gnomad_vcf)
def gnomad_tbi = file("${params.gnomad_vcf}.tbi")

if( !gnomad_vcf.exists() ) {
  error "gnomAD VCF not found: ${gnomad_vcf}"
}
if( !gnomad_tbi.exists() ) {
  error "gnomAD index (.tbi) not found: ${gnomad_tbi}\nCreate it with: tabix -p vcf ${gnomad_vcf}"
}
if( !gnomad_vcf.canRead() || !gnomad_tbi.canRead() ) {
  error "gnomAD files exist but are not readable: ${gnomad_vcf} / ${gnomad_tbi}"
}

/* ---------- basic/safe defaults (silence WARNs) ---------- */
if( !params.containsKey('outdir') )                 params.outdir            = 'results'
if( !params.containsKey('reads') )                  params.reads             = null
if( !params.containsKey('ref_fa') )                 params.ref_fa            = null
if( !params.containsKey('region') )                 params.region            = null
if( !params.containsKey('max_cpus') )               params.max_cpus          = 8
if( !params.containsKey('max_mem') )                params.max_mem           = '24 GB'
if( !params.containsKey('max_time') )               params.max_time          = '90h'
if( !params.containsKey('report_coding_only') )     params.report_coding_only = true
params.enable_benchmark = params.enable_benchmark ?: true
if( !params.containsKey('report_max_per_gene') )   params.report_max_per_gene = 50
if( !params.containsKey('report_max_rows') )       params.report_max_rows     = 100000
if( !params.containsKey('acc_metrics_tsv') )       params.acc_metrics_tsv     = null
if( !params.containsKey('trace_file') )             params.trace_file       = "${params.outdir ?: 'results'}/trace.tsv"
/* ---------- silence undefined-parameter WARNs (these are READ later) ---------- */
if( !params.containsKey('arm_script') )        params.arm_script        = null
if( !params.containsKey('bench_metrics_tsv') ) params.bench_metrics_tsv = null
if( !params.containsKey('knowledge_tsv') )     params.knowledge_tsv     = null
if( !params.containsKey('pon_vcf') )               params.pon_vcf               = null
if( !params.containsKey('somatic_intervals_bed') ) params.somatic_intervals_bed = null

/* --- Tumor/Normal params (optional) --- */
if( !params.containsKey('tumor_id') )               params.tumor_id            = null
if( !params.containsKey('normal_id') )              params.normal_id           = null


/* --- VEP defaults --- */
if( !params.containsKey('vep_cache') )              params.vep_cache         = null   // you used host path before
if( !params.containsKey('vep_cache_version') )      params.vep_cache_version = 110
if( !params.containsKey('vep_fasta') )              params.vep_fasta         = null
if( !params.containsKey('vep_assembly') )           params.vep_assembly      = 'GRCh38'

/* --- reporting defaults --- */
if( !params.containsKey('patient_id') )             params.patient_id        = 'PATIENT01'
if( !params.containsKey('pharmcat_jar') )           params.pharmcat_jar      = null
if( !params.containsKey('knowledge_tsv') )          params.knowledge_tsv     = null
if( !params.containsKey('pathway_db') )             params.pathway_db        = null
if( !params.containsKey('gene_domains') )           params.gene_domains      = null
if( !params.containsKey('civic_offline') )          params.civic_offline     = null
if( !params.containsKey('enable_civic') )           params.enable_civic      = true

/* --- suppress WARNs for unused params --- */
if( !params.containsKey('mcode') )                  params.mcode             = null
if( !params.containsKey('oncokb_token') )           params.oncokb_token       = System.getenv('ONCOKB_TOKEN')
if( !params.containsKey('ncbi_key') )               params.ncbi_key          = null
if( !params.containsKey('ncbi_email') )             params.ncbi_email        = null
if( !params.containsKey('enable_online') )          params.enable_online     = null
if( !params.containsKey('targets_bed') )            params.targets_bed       = null

/* --- panel options (optional) --- */
if( !params.containsKey('panel_name') )             params.panel_name        = null
if( !params.containsKey('assay') )                  params.assay             = (params.panel_bed ? 'panel' : 'wgs')
if( !params.containsKey('discard_trimmed_pass') )   params.discard_trimmed_pass = true

/* --- report selection --- */
params.empty_json = params.empty_json ?: "${baseDir}/.empty.json"
params.empty_tsv  = params.empty_tsv  ?: "${baseDir}/.empty.tsv"


/* ===================== LOCAL RESOURCE DEFAULTS (project-relative) ===================== */

// pick the first path that exists
def pickExisting = { List<String> paths ->
  paths.find { p -> p && file(p).exists() }
}

// set params[key] only if user did not provide it (or provided blank)
def setDefaultPath = { String key, List<String> candidates, boolean required=false ->
  boolean missing = (!params.containsKey(key) || params[key] == null || params[key].toString().trim() == '')
  if( missing ) {
    def chosen = pickExisting(candidates)
    if( chosen ) params[key] = chosen
  }
  if( required ) {
    def v = params[key]
    if( !v || !file(v.toString()).exists() ) {
      error "Required resource missing: ${key} (tried: ${candidates.join(' | ')})"
    }
  }
}

// helper: require tabix index for .vcf.gz
def requireTabix = { String key ->
  def v = params[key]
  if( v && v.toString().endsWith('.vcf.gz') ) {
    def idx = "${v}.tbi"
    if( !file(idx).exists() ) {
      error "Tabix index missing for ${key}: expected ${idx}"
    }
  }
}

/* ---- fixed resources under your project ---- */
setDefaultPath('germline_resource', [
  "${baseDir}/resources/gnomad/somatic-hg38_af-only-gnomad.hg38.vcf.gz"
], true)
requireTabix('germline_resource')

setDefaultPath('common_variants', [
  "${baseDir}/resources/gatk_somatic_resources/hg38/small_exac_common_3.hg38.vcf.gz"
], true)
requireTabix('common_variants')

setDefaultPath('public_pon', [
  "${baseDir}/resources/PON_hg38/1000g_pon.hg38.vcf.gz"
], true)
requireTabix('public_pon')

/* ---- PoN selection: default to your real PoN, else fall back to public_pon ---- */
def myPon = "${baseDir}/resources/pon/pon.vcf.gz"   // your real PoN

// If user did NOT explicitly provide --pon_vcf, pick defaults:
//   1) myPon if exists
//   2) public_pon otherwise
if( !params.pon_vcf || params.pon_vcf.toString().trim() == '' ) {
  params.pon_vcf = file(myPon).exists() ? myPon : (params.public_pon ?: null)
} else {
  // user provided --pon_vcf; keep it
  params.pon_vcf = params.pon_vcf.toString()
}

// sanity: require index if it's a bgzipped vcf
if( params.pon_vcf && params.pon_vcf.toString().endsWith('.vcf.gz') ) {
  def idx = "${params.pon_vcf}.tbi"
  if( !file(idx).exists() ) error "Tabix index missing for pon_vcf: expected ${idx}"
}


/* ---- other “always same” inputs you were passing (optional, no hard fail) ---- */
setDefaultPath('vep_cache',        [ "${baseDir}/vep_cache" ], false)
setDefaultPath('knowledge_tsv',    [ "${baseDir}/kb_global_template.tsv" ], false)
setDefaultPath('pathway_db',       [ "${baseDir}/resources/pathway_gene_map_v1.json" ], false)
setDefaultPath('msk_hotspot_keys', [ "${baseDir}/resources/hotspots/msk_cancer_hotspots.keys.tsv" ], false)
setDefaultPath('civic_offline',    [ "${baseDir}/resources/CivicDB/civic_offline.json" ], false)

/* ---- panel + intervals (optional defaults) ---- */
setDefaultPath('panel_bed', [ "${baseDir}/resources/panel/panel.hg38.coding+10bp.bed" ], false)
if( !params.somatic_intervals_bed || params.somatic_intervals_bed.toString().trim() == '' ) {
  params.somatic_intervals_bed = (params.panel_bed ?: null)
}

// hotspot behavior (keep your original intent)
if( !params.containsKey('use_hotspots') ) {
  params.use_hotspots = !(params.oncokb_token && params.oncokb_token.toString().trim())
}
/* --- report filtering defaults --- */
if( !params.containsKey('report_use_csq') )    params.report_use_csq    = true
if( !params.containsKey('report_csq_regex') )  params.report_csq_regex   =
  'missense_variant|frameshift_variant|stop_gained|stop_lost|start_lost|' +
  'splice_(acceptor|donor)_variant|inframe_(insertion|deletion)|' +
  'protein_altering_variant|synonymous_variant|' +
  'non_coding_transcript_exon_variant|3_prime_UTR_variant|5_prime_UTR_variant|' +
  'NMD_transcript_variant'

/* --- UMI/panel QC & resource knobs --- */
if( !params.containsKey('umi_regex') )         params.umi_regex         = null
if( !params.containsKey('mosdepth_bin') )      params.mosdepth_bin      = 'mosdepth'
if( !params.containsKey('panel_min_cov') )     params.panel_min_cov     = 100
if( !params.containsKey('panel_lod') )         params.panel_lod         = 0.02
if( !params.containsKey('bwa_k') )             params.bwa_k             = 20000000
if( !params.containsKey('bwa_threads') )       params.bwa_threads       = 4
if( !params.containsKey('enable_evidence') )   params.enable_evidence   = true
if( !params.containsKey('report_search_limit') ) params.report_search_limit = 500
if( !params.containsKey('report_keep_impact') )  params.report_keep_impact  = 'HIGH,MODERATE'
if( !params.containsKey('html_max_rows') )           params.html_max_rows         = 3000
if( !params.containsKey('drop_pharmcat_in_html') )   params.drop_pharmcat_in_html = true
/* --- Arm-level CN calling (defaults) --- */
if( !params.containsKey('arm_gain_log2') )  params.arm_gain_log2 =  0.25
if( !params.containsKey('arm_loss_log2') )  params.arm_loss_log2 = -0.25
if( !params.containsKey('arm_cover_min') )  params.arm_cover_min =  0.60
if( !params.containsKey('cytoband') )       params.cytoband      = null   // auto-pick below if null

/* ---------- container images ---------- */
def IMG_FASTQC   = 'quay.io/biocontainers/fastqc:0.11.9--0'
def IMG_FASTP    = 'quay.io/biocontainers/fastp:0.23.4--h5f740d0_0'
def IMG_BWA      = 'quay.io/biocontainers/bwa:0.7.17--he4a0461_11'
def IMG_SAMTOOLS = 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'
def IMG_BCFTOOLS = 'quay.io/biocontainers/bcftools:1.17--h3cc50cf_1'
def IMG_VEP      = 'ensemblorg/ensembl-vep:release_110.1'
def IMG_JAVA     = 'openjdk:17-jdk-slim'
def IMG_PY       = 'python:3.11'

/* ---------- report builder ---------- */
params.build_html           = params.build_html ?: "${baseDir}/HTML/build_html.py"
params.discard_trimmed_pass = params.discard_trimmed_pass ?: false

/* ============================== WORKFLOW ================================== */

workflow CLINICAL_CORE {

  take:
    reads_in   // (sid, r1, r2)
    ref_fa_in  // singleton (string or Path)

  main:
    /* ---- alias take: ports to locals to avoid scope collisions ---- */
    def ch_reads = reads_in
    def ch_refin = ref_fa_in
        // ---- Panel BED singleton channel (declare ONCE, before use) ----
    def panel_bed_ch = params.panel_bed \
      ? Channel.fromPath(params.panel_bed, checkIfExists: true) \
      : Channel.empty()


    /* ---------- reference singletons ---------- */
    def ref_fa_ch = ch_refin.map { file(it as String) }
    def faidx     = REF_FAIDX(ref_fa_ch)
    def ref_pair  = faidx.ref_pair
    def ref_dictp = REF_DICT(ref_fa_ch)   // Build .dict
    def bwa       = BWA_INDEX(ref_fa_ch)
    def bwa_dir   = bwa.index_dir

    /* ---------- QC + trimming ---------- */
    FASTQC_RAW(ch_reads)
    def trim    = TRIM_FASTP(ch_reads)
    FASTQC_TRIM(trim.reads)

    /* ---------- align → sort/index ---------- */
    def aln_sam = ALIGN_BWA(trim.reads, bwa_dir)
    def sorted  = SAMTOOLS_SORT_INDEX(aln_sam, ref_pair)

    /* ---------- mark dups ---------- */
    def md            = MARKDUPS(sorted, ref_pair)
    def md_align_ch   = md.align          // (sid, .cram, .crai)
    def md_metrics_ch = md.metrics

    /* ---------- ensure CRAI ---------- */
    def md_idx_in = md_align_ch.map { sid, cram, crai -> tuple(sid, cram) }
    def md_idx    = REINDEX_CRAM(md_idx_in)   // (sid, cram, cram.crai)

    /* ---------- optional Tumor/Normal pairing ---------- */
    def isTN   = (params.tumor_id && params.normal_id)
    def tn_pair = Channel.empty()
   
    if (isTN) {
      // Pick the exact tumor/normal samples by ID
      def tumor_ch  = md_align_ch.filter { sid, cram, crai -> sid.toString() == params.tumor_id.toString() }
      def normal_ch = md_align_ch.filter { sid, cram, crai -> sid.toString() == params.normal_id.toString() }

      // Use a dummy key to join the singletons → one 6-field tuple event
      def t_kv = tumor_ch .map { sid, cram, crai -> tuple(1, sid, cram, crai) }
      def n_kv = normal_ch.map { sid, cram, crai -> tuple(1, sid, cram, crai) }

      tn_pair = t_kv.join(n_kv, by: 0)
                    .map { _k, tsid, t_cram, t_crai, nsid, n_cram, n_crai ->
                      tuple(tsid, t_cram, t_crai, nsid, n_cram, n_crai)
                    }
    }


    /* ---------- Cytoband (singleton) ---------- */
    def asm   = (params.vep_assembly ?: 'GRCh38').toString().toUpperCase()
    def cytoGuess = (asm.contains('37') || asm.contains('HG19'))
      ? file("${baseDir}/resources/cytoband/cytoBand_hg19.tsv")
      : file("${baseDir}/resources/cytoband/cytoBand_hg38.tsv")
    def cytoband_f  = file(params.cytoband ?: cytoGuess.toString())
    def cytoband_ch = Channel.value( cytoband_f.exists() ? cytoband_f : file(params.empty_tsv) )

    // ARM script singleton (added to fix undefined variable)
    // ARM script singleton (try params, then scripts/ fallback)
    def arm_script_file = file(params.arm_script ?: "${baseDir}/resources/arm_events/call_arm_events.py")
    if( !arm_script_file.exists() ) {
      arm_script_file = file("${baseDir}/scripts/call_arm_events.py")
    }
    def arm_script_ch = Channel.value(arm_script_file)

    // ---------------- PoN (singleton tuple) ----------------
    def pon_vcf_f = file(params.pon_vcf ?: params.empty_tsv)
    if( !pon_vcf_f.exists() ) pon_vcf_f = file(params.empty_tsv)

    def pon_tbi_f = file(params.pon_vcf ? "${params.pon_vcf}.tbi" : params.empty_tsv)
    if( !pon_tbi_f.exists() ) pon_tbi_f = file(params.empty_tsv)

    def pon_ch = Channel.value( tuple(pon_vcf_f, pon_tbi_f) )

    // ------------- germline resource (singleton tuple) -------------
    def germ_vcf_f = file(params.germline_resource ?: params.empty_tsv)
    if( !germ_vcf_f.exists() ) germ_vcf_f = file(params.empty_tsv)

    def germ_tbi_f = file(params.germline_resource ? "${params.germline_resource}.tbi" : params.empty_tsv)
    if( !germ_tbi_f.exists() ) germ_tbi_f = file(params.empty_tsv)

    def germ_ch = Channel.value( tuple(germ_vcf_f, germ_tbi_f) )

    // ---------------- intervals (singleton path) ----------------
    def intervals_f = file(params.somatic_intervals_bed ?: params.empty_tsv)
    if( !intervals_f.exists() ) intervals_f = file(params.empty_tsv)

    def intervals_ch = Channel.value(intervals_f)
    // ---------- hard-required gnomAD (Mutect2 germline resource) ----------
    params.gnomad_vcf = params.gnomad_vcf ?: "${baseDir}/resources/gnomad/somatic-hg38_af-only-gnomad.hg38.vcf.gz"

    def gnomad_vcf_f = file(params.gnomad_vcf)
    def gnomad_tbi_f = file("${params.gnomad_vcf}.tbi")

    if( !gnomad_vcf_f.exists() || !gnomad_tbi_f.exists() ) {
      log.error "[FATAL] Missing gnomAD VCF or index:\n  VCF: ${gnomad_vcf_f}\n  TBI: ${gnomad_tbi_f}"
      System.exit(1)
    }


    // ---------------- somatic call ----------------
    def somatic = CALL_SOMATIC_MUTECT2(tn_pair, ref_pair, ref_dictp.ref_dict, pon_ch, germ_ch, intervals_ch)
    def germline = CALL_VARIANTS_CRAM(md_idx, ref_pair, panel_bed_ch)
    def called   = isTN ? somatic.vcf : germline.vcf


    /* ---------- optional VEP ---------- */
    def vep_in  = called.map { sid, vcf_gz, tbi -> tuple(sid, vcf_gz) }
    def vcf_ann = params.vep_cache ? VEP_ANNOTATE(vep_in) : called
    
    /* ---------- filter for HTML ---------- */
    def fhtml = FILTER_FOR_HTML(vcf_ann)

    // singleton: hotspot keys (ok if empty/missing)
    def hotspot_keys_ch = Channel.value( file(params.msk_hotspot_keys ?: params.empty_tsv) )

    // annotate the reportready VCF with MSK hotspots
    def annotated_vcf_ch = params.use_hotspots \
      ? ANNOTATE_MSK_HOTSPOTS( fhtml.combine(hotspot_keys_ch).map { sid, vcf, tbi, keys -> tuple(sid, vcf, tbi, keys) } ) \
      : fhtml


    // convenience map for steps that accept (sid, vcf)
    def annotated_vcf_2col = annotated_vcf_ch.map { sid, vcf_gz, tbi -> tuple(sid, vcf_gz) }

    /* ---------- PharmCAT (optional) ---------- */
    def pharm_json = params.pharmcat_jar
      ? PHARMCAT(annotated_vcf_2col).pharmcat_json
      : annotated_vcf_2col.map { sid, _ -> tuple(sid, file(params.empty_json)) }


    /* ---------- knowledge (singleton) ---------- */
    def kb_ch       = Channel.value(file(params.knowledge_tsv ?: params.empty_tsv))
    def hotspots_ch = Channel.value(file(params.msk_hotspot_keys ?: params.empty_tsv))

    /* ---------- pathway + gene domains (singletons) ---------- */
    def pdb_f  = file(params.pathway_db   ?: "${baseDir}/resources/pathway_gene_map_v1.json")
    def gd_f   = file(params.gene_domains ?: "${baseDir}/resources/gene_domains.json")
    def pdb_ch = Channel.value(pdb_f.exists() ? pdb_f : file(params.empty_json))
    def gd_ch  = Channel.value(gd_f.exists()  ? gd_f  : file(params.empty_json))


    /* ---------- Panel QC (optional) ---------- */
    def mos_regions_ch = Channel.empty()   // <— add this line

    def panel_qc_ch
    if (params.panel_bed) {
      def mos = PANEL_QC_MOSDEPTH(md_align_ch, file(params.panel_bed), ref_pair)

      // expose regions for downstream CN arms:
      mos_regions_ch = mos.regions         // <— add this line

      // one joined tuple per sample...
      def qc_in = mos.summary
        .join(md_metrics_ch, by: 0)
        .join(trim.json,    by: 0)
        .join(md_align_ch,  by: 0)
        .join(mos.thresholds, by: 0)
        .join(mos.regions,    by: 0)
        .map { sid, summary, metrics, fastp, cram, crai, thr, reg ->
          tuple(sid, summary, metrics, fastp, cram, crai, thr, reg)
        }
      panel_qc_ch = QC_AGGREGATE(qc_in, file(params.panel_bed), ref_pair)
    } else {
      panel_qc_ch = md_align_ch.map { sid, cram, crai -> tuple(sid, file(params.empty_json)) }
    }

    /* ---------- JSON + mCODE ---------- */
    def merged  = annotated_vcf_2col.join(pharm_json, by: 0)
                                    .map { sid, vcf_gz, pharmcat_json -> tuple(sid, vcf_gz, pharmcat_json) }
    def merged2 = merged.join(panel_qc_ch, by: 0)
                        .map { sid, vcf_gz, pharmcat_json, panel_qc_json ->
                          tuple(sid, vcf_gz, pharmcat_json, panel_qc_json)
                        }

    def rep = MAKE_JSON_SUMMARY(merged2, kb_ch, hotspot_keys_ch, pdb_ch)


    //  civic singleton
    def civic_ch = Channel.value( file(params.civic_offline) )

    //  enrich the report.json → <sid>.report.civic.json
    def civic_summary = CIVIC_OFFLINE_ENRICH( rep.report_json, civic_ch ).civic_report


    //  use the enriched JSON for HTML (the process has a `when:` guard)
    def json_for_html = params.enable_civic ? civic_summary : rep.report_json


    //  build slim/full from the chosen JSON
    def slim = MAKE_SLIM_JSON(json_for_html)
    /* ==================== ARM EVENTS INTEGRATION ==================== */
    def arm_calls_json_ch = Channel.empty()

    if (isTN && params.panel_bed) {
      // pick tumor/normal mosdepth regions
      def t_reg = mos_regions_ch
                    .filter { sid, _ -> sid.toString() == params.tumor_id.toString() }
                    .map    { sid, reg -> tuple(1, sid, reg) }
      def n_reg = mos_regions_ch
                    .filter { sid, _ -> sid.toString() == params.normal_id.toString() }
                    .map    { sid, reg -> tuple(1, sid, reg) }

      // pair tumor/normal regions → (tsid, t_regions, nsid, n_regions)
      def tn_regs = t_reg.join(n_reg, by: 0)
                        .map { _k, tsid, t_regions, nsid, n_regions ->
                          tuple(tsid, t_regions, nsid, n_regions)
                        }

      // segments from mosdepth tumor vs normal
      def tn_segments = MAKE_TN_SEGMENTS_FROM_MOSDEPTH(tn_regs)

      // tn_segments emits: (tsid, nsid, segments_tsv)
      // make a staged file object from params
      def arm_script_f  = file(params.arm_script)

    

      // (tsid, nsid, segments_tsv) + (arm_script) -> (tsid, segments_tsv, arm_script)
      def arm_in = tn_segments
        .combine(arm_script_ch)
        .map { tsid, nsid, seg_tsv, arm_py -> tuple(tsid, seg_tsv, arm_py) }

      def arm_calls = CALL_ARM_EVENTS(arm_in, cytoband_ch)

      // keep JSON for merge
      arm_calls_json_ch = arm_calls.arms.map { sid, arm_json, _tags -> tuple(sid, arm_json) }

    }
    /* ================== END ARM EVENTS INTEGRATION ================== */

    /* (3) Augment tumor slim/full JSON with ARM summary (when available) */
    def slim_for_merge_tumor
    if (isTN && params.panel_bed) {
      // grab tumor-only slim/full
      def slim_t = slim.filter { sid, _s, _f -> sid.toString() == params.tumor_id.toString() }
      // join with arm JSON; if arm JSON is missing, fall back to pass-through
      def t_join  = slim_t.join(arm_calls_json_ch, by: 0)
                          .map { sid, slim_json, full_json, arm_json ->
                            tuple(sid, slim_json, full_json, arm_json)
                          }
      def t_aug = AUGMENT_JSON_ARM_EVENTS(t_join).json
      // -> (sid, sid.report.slim.arm.json, sid.report.full.arm.json)
      slim_for_merge_tumor = t_aug
    } else {
      // no TN or no panel → just reuse the original tumor slim/full
      slim_for_merge_tumor = slim.filter { sid, _s, _f -> sid.toString() == params.tumor_id?.toString() }
    }
    /* ================== END ARM EVENTS INTEGRATION ================== */

    /* ---------- choose HTML input (TN combined vs per-sample) ---------- */
    def html_in
    if (isTN) {
      // tumor QC and normal QC are the same as before
      def qc_t = panel_qc_ch.filter { sid, _ -> sid.toString() == params.tumor_id.toString() }
      def qc_n = panel_qc_ch.filter { sid, _ -> sid.toString() == params.normal_id.toString() }

      // prefer augmented tumor slim/full (has arm_cn in summary); fallback handled above
      def pair_in = slim_for_merge_tumor
        .join(qc_t, by: 0)                                  // (sid, slim, full, t_qc)
        .map { sid, s, f, tqc -> tuple(1, sid, s, f, tqc) }
        .join(qc_n.map { sid, nqc -> tuple(1, sid, nqc) }, by: 0)
        .map { _k, tsid, s, f, tqc, _nsid, nqc ->
          tuple(tsid, s, f, tqc, params.normal_id.toString(), nqc)
        }

      def mergedTN = MERGE_TN_FOR_HTML(pair_in)             // -> (tumor_id, slim_json, full_json)
      html_in = mergedTN.map { tsid, slim_json, full_json ->
        tuple("${tsid}_vs_${params.normal_id}", slim_json, full_json)
      }
    } else {
      html_in = slim.map { sid, slim_json, full_json -> tuple(sid, slim_json, full_json) }
    }

    /* ---------- HTML report ---------- */
    def builder = Channel.fromPath(params.build_html ?: "${baseDir}/HTML/build_html.py", checkIfExists:true)

    // ADD THESE TWO LINES HERE
    def civ_f = file(params.civic_offline ?: "${baseDir}/resources/CivicDB/civic_offline.json")
    def civic_for_build = Channel.value(civ_f.exists() ? civ_f : file(params.empty_json))

    // and REPLACE the old call with this one (note the extra arg at the end)
    def built         = BUILD_HTML(html_in, builder, pdb_ch, gd_ch, civic_for_build)
    def to_inject     = built.html.map { f -> tuple(f.getName(), f) }
    def builtInjected = INJECT_PAIRED_QC_HTML(to_inject)
    // ---------------- BENCHMARK REPORT (runs once, after final HTML exists) ----------------
    // ---------------- BENCHMARK REPORT (runs once, after final HTML exists) ----------------
    def html_ch = builtInjected.html

    if( params.enable_benchmark ) {
      // gate: waits until final HTML exists (i.e., the pipeline is effectively finished)
      def gate_ch  = html_ch.take(1).map { 1 }

      // also gate the trace path emission so downstream processes don't start too early
      def trace_ch = html_ch.take(1).map { file(params.trace_file) }

      // build bench metrics from the trace (single output channel assumed)
      def bench_metrics_ch = MAKE_BENCH_METRICS(trace_ch)

      // now BENCHMARK_REPORT gets exactly 4 inputs
      // in workflow scope:
      def acc_ch = Channel.value( file(params.acc_metrics_tsv ?: params.empty_tsv) )

      BENCHMARK_REPORT(
        gate_ch,
        trace_ch,
        bench_metrics_ch,
        acc_ch,
        Channel.value(params.outdir ?: 'results')
      )

    }

  emit:
    aligned          = sorted
    bam_ch           = md_align_ch
    vcf_ch           = annotated_vcf_ch
    panel_ch         = panel_qc_ch
    report_json      = rep.report_json
    mcode_json       = rep.mcode_json
    report_html      = builtInjected.html
    report_full_json = built.full_json
}


/* ============================== PROCESSES =================================
 * (ordered roughly by pipeline flow)
 * ========================================================================== */

/* ---------- Reference prep ---------- */
process REF_FAIDX {
  tag { ref_fa.baseName }
  container "${IMG_SAMTOOLS}"
  cpus 5
  memory '12 GB'
  time '30m'
  stageInMode 'copy'

  input:
    path ref_fa
  output:
    tuple path(ref_fa), path("${ref_fa}.fai"), emit: ref_pair
    path("${ref_fa}.gzi"), optional: true,     emit: gzi
  shell:
  """
  set -euo pipefail
  samtools faidx "!{ref_fa}"
  """
}

process REF_DICT {
  tag { ref_fa.baseName }
  container "${IMG_SAMTOOLS}"
  cpus 1
  memory '1 GB'
  time '10m'
  stageInMode 'copy'

  input:
    path ref_fa

  output:
    path("*.dict"), emit: ref_dict

  shell:
  '''
  set -eo pipefail
  out="hg38.dict"
  if [ ! -s "$out" ]; then
    samtools dict -o "$out" "hg38.fa"
  fi
  '''
}

process BWA_INDEX {
  tag { reference_fasta.baseName }
  container "${IMG_BWA}"
  cpus 4
  memory '8 GB'
  stageInMode 'copy'

  input:
    path reference_fasta
  output:
    path 'bwa_index', emit: index_dir
  shell:
  """
  set -euo pipefail
  mkdir -p bwa_index
  bwa index -p bwa_index/ref "!{reference_fasta}"
  """
}

/* ---------- QC & trimming ---------- */
process FASTQC_RAW {
  container "${IMG_FASTQC}"
  stageInMode 'symlink'
  input:
    tuple val(sid), path(r1), path(r2)
  output:
    path "${sid}_R1_fastqc.html"
    path "${sid}_R1_fastqc.zip"
    path "${sid}_R2_fastqc.html"
    path "${sid}_R2_fastqc.zip"
  shell:
  """
  set -eo pipefail
  if [ -s !{r2} ] && [ "\$(readlink -f !{r1})" != "\$(readlink -f !{r2})" ]; then
    fastqc -q -o ./ !{r1} !{r2}
  else
    fastqc -q -o ./ !{r1}
    : > "!{sid}_R2_fastqc.html"
    : > "!{sid}_R2_fastqc.zip"
  fi
  r1base=\$(basename "!{r1}"); s=\${r1base%.gz}; s=\${s%.fastq}; s=\${s%.fq}
  mv "\${s}_fastqc.html" "!{sid}_R1_fastqc.html" || true
  mv "\${s}_fastqc.zip"  "!{sid}_R1_fastqc.zip"  || true
  if ls *_fastqc.html >/dev/null 2>&1; then
    r2base=\$(basename "!{r2}"); s=\${r2base%.gz}; s=\${s%.fastq}; s=\${s%.fq}
    mv "\${s}_fastqc.html" "!{sid}_R2_fastqc.html" || true
    mv "\${s}_fastqc.zip"  "!{sid}_R2_fastqc.zip"  || true
  fi
  """
}

process TRIM_FASTP {
  tag { sid }
  container "${IMG_FASTP}"
  cpus 6
  memory '24 GB'
  time '6h'
  input:
    tuple val(sid), path(r1), path(r2)
  output:
    tuple val(sid), path("${sid}.R1.trim.fastq.gz"), path("${sid}.R2.trim.fastq.gz"), emit: reads
    tuple val(sid), path("${sid}.fastp.json"), emit: json
    tuple val(sid), path("${sid}.fastp.html"), emit: html
    tuple val(sid), path("${sid}.fastp.log") , emit: log
  shell:
  """
  set -euo pipefail
  if [ -s "!{r2}" ] && [ "\$(readlink -f "!{r1}")" != "\$(readlink -f "!{r2}")" ]; then
    fastp \
      -i "!{r1}" -I "!{r2}" \
      -o "!{sid}.R1.trim.fastq.gz" -O "!{sid}.R2.trim.fastq.gz" \
      -w !{task.cpus} -q 20 -l 30 \
      -j "!{sid}.fastp.json" -h "!{sid}.fastp.html" \
      2> >(tee "!{sid}.fastp.log" >&2)
  else
    fastp \
      -i "!{r1}" \
      -o "!{sid}.R1.trim.fastq.gz" \
      -w !{task.cpus} -q 20 -l 30 \
      -j "!{sid}.fastp.json" -h "!{sid}.fastp.html" \
      2> >(tee "!{sid}.fastp.log" >&2)
    : > "!{sid}.R2.trim.fastq.gz"
  fi
  """
}

process FASTQC_TRIM {
  container "${IMG_FASTQC}"
  stageInMode 'symlink'
  input:
    tuple val(sid), path(r1), path(r2)
  output:
    path "${sid}_R1_fastqc.html"
    path "${sid}_R1_fastqc.zip"
    path "${sid}_R2_fastqc.html"
    path "${sid}_R2_fastqc.zip"
  shell:
  """
  set -euo pipefail
  if [ -s "!{r2}" ] && [ "\$(readlink -f "!{r1}")" != "\$(readlink -f "!{r2}")" ]; then
    fastqc -q -o ./ "!{r1}" "!{r2}"
    r1base=\$(basename "!{r1}"); r1stem=\${r1base%.gz}; r1stem=\${r1stem%.fastq}; r1stem=\${r1stem%.fq}
    r2base=\$(basename "!{r2}"); r2stem=\${r2base%.gz}; r2stem=\${r2stem%.fastq}; r2stem=\${r2stem%.fq}
    mv "\${r1stem}_fastqc.html" "!{sid}_R1_fastqc.html"
    mv "\${r1stem}_fastqc.zip"  "!{sid}_R1_fastqc.zip"
    mv "\${r2stem}_fastqc.html" "!{sid}_R2_fastqc.html"
    mv "\${r2stem}_fastqc.zip"  "!{sid}_R2_fastqc.zip"
  else
    fastqc -q -o ./ "!{r1}"
    r1base=\$(basename "!{r1}"); r1stem=\${r1base%.gz}; r1stem=\${r1stem%.fastq}; r1stem=\${r1stem%.fq}
    mv "\${r1stem}_fastqc.html" "!{sid}_R1_fastqc.html"
    mv "\${r1stem}_fastqc.zip"  "!{sid}_R1_fastqc.zip"
    : > "!{sid}_R2_fastqc.html"
    : > "!{sid}_R2_fastqc.zip"
  fi
  """
}

/* ---------- Alignment & post-processing ---------- */
process ALIGN_BWA {
  tag { sid }
  container "${IMG_BWA}"
  cpus { params.bwa_threads }
  memory { params.max_mem }
  time { params.max_time }
  stageInMode 'symlink'
  input:
    tuple val(sid), path(r1), path(r2)
    path bwa_dir
  output:
    tuple val(sid), path("${sid}.aln.sam.gz")
  shell:
  """
  set -euo pipefail
  RG_STR='@RG\\tID:!{sid}\\tSM:!{sid}\\tPL:ILLUMINA'
  if [ -s "!{r2}" ] && [ "\$(readlink -f "!{r1}")" != "\$(readlink -f "!{r2}")" ]; then
    bwa mem -K !{params.bwa_k} -Y -t !{task.cpus} -R "\$RG_STR" \
      "!{bwa_dir}/ref" "!{r1}" "!{r2}" | gzip -1 > "!{sid}.aln.sam.gz"
  else
    bwa mem -K !{params.bwa_k} -Y -t !{task.cpus} -R "\$RG_STR" \
      "!{bwa_dir}/ref" "!{r1}" | gzip -1 > "!{sid}.aln.sam.gz"
  fi
  """
}

process SAMTOOLS_SORT_INDEX {
  tag { sid }
  container "${IMG_SAMTOOLS}"
  cpus { params.max_cpus }
  memory { params.max_mem }
  time { params.max_time }
  stageInMode 'copy'
  publishDir "${params.outdir}/align", mode: 'copy', overwrite: true
  input:
    tuple val(sid), path(sam_gz)
    tuple path(ref_fa), path(ref_fai)
  output:
    tuple val(sid), path("${sid}.sorted.cram"), path("${sid}.sorted.cram.crai")
  shell:
  """
  set -euo pipefail
  mkdir -p tmp
  gzip -cd "!{sam_gz}" \
  | samtools sort -@ !{task.cpus} -m 768M \
      -O CRAM --reference "!{ref_fa}" \
      -T "tmp/!{sid}" -o "!{sid}.sorted.cram" -
  samtools index -@ !{task.cpus} "!{sid}.sorted.cram"
  """
}

process MARKDUPS {
  tag { sid }
  container 'broadinstitute/gatk:4.5.0.0'
  cpus { params.max_cpus }
  memory '24 GB'
  time '8h'
  stageInMode 'copy'
  input:
    tuple val(sid), path(cram), path(crai)
    tuple path(ref_fa), path(ref_fai)
  output:
    tuple val(sid), path("${sid}.markdup.cram"), path("${sid}.markdup.cram.crai"), emit: align
    tuple val(sid), path("${sid}.markdup.metrics.txt"), emit: metrics
  shell:
  """
  set -euo pipefail
  TMPDIR="\${TMPDIR:-\$PWD/tmp}"; mkdir -p "\$TMPDIR"
  gatk --java-options "-Djava.io.tmpdir=\$TMPDIR" MarkDuplicates \
    -I "!{cram}" -O "!{sid}.markdup.cram" -M "!{sid}.markdup.metrics.txt" \
    --REFERENCE_SEQUENCE "!{ref_fa}" --CREATE_INDEX true --TMP_DIR "\$TMPDIR"
  if [ ! -f "!{sid}.markdup.cram.crai" ]; then
    gatk BuildBamIndex -I "!{sid}.markdup.cram" || true
  fi
  if [ ! -f "!{sid}.markdup.cram.crai" ] && [ -f "!{sid}.markdup.cram.bai" ]; then
    ln -sf "!{sid}.markdup.cram.bai" "!{sid}.markdup.cram.crai"
  fi
  ls -l "!{sid}.markdup.cram" "!{sid}.markdup.cram.crai"
  """
}

process REINDEX_CRAM {
  tag { sid }
  container "${IMG_SAMTOOLS}"
  cpus 3
  time '1h'
  stageInMode 'link'
  input:
    tuple val(sid), path(cram)
  output:
    tuple val(sid), path(cram), path("${cram}.crai")
  shell:
  """
  set -euo pipefail
  samtools index -@ !{task.cpus} -c "!{cram}"
  """
}

/* ---------- Variant calling ---------- */
process CALL_VARIANTS_CRAM {
  tag { sid }
  container "${IMG_BCFTOOLS}"
  cpus 3
  memory '6 GB'
  time '8h'
  stageInMode 'copy'
  publishDir "${params.outdir}/vcf", mode:'copy', overwrite:true
  input:
    tuple val(sid), path(cram), path(crai)
    tuple path(ref_fa), path(ref_fai)
    path panel_bed

  output:
    tuple val(sid), path("${sid}.vcf.gz"), path("${sid}.vcf.gz.tbi"), emit: vcf
  shell:
  '''
  set -euo pipefail
  set -x

  BED_IN="!{panel_bed}"
  BED_ARG=""

  if [ -n "$BED_IN" ] && [ -s "$BED_IN" ]; then
    # get contig names from CRAM header
    cut -f1 "!{ref_fai}" | sort -u > __ref.contigs

    # first contig in BED (ignore comments)
    BED_CONTIG=$(awk 'NF>=1 && $1!~/^#/ {print $1; exit}' "$BED_IN")

    # if BED contig style doesn't match CRAM contigs, convert BED contigs
    if ! grep -qxF "$BED_CONTIG" __ref.contigs; then
      if grep -q '^chr' __ref.contigs; then
        # CRAM uses chr*, BED doesn't -> add chr
        awk 'BEGIN{OFS="\t"} $1!~/^#/{ $1 = ($1 ~ /^chr/ ? $1 : "chr"$1) } {print}' "$BED_IN" > __bed.tmp
      else
        # CRAM uses no chr, BED uses chr* -> remove chr
        awk 'BEGIN{OFS="\t"} $1!~/^#/{ sub(/^chr/,"",$1) } {print}' "$BED_IN" > __bed.tmp
      fi
      BED_IN="__bed.tmp"
    fi

    # sort + dedup using reference contig order (portable sort; no GNU -V)
    # rank contigs by appearance in hg38.fa.fai
    awk -v OFS='\t' '
      NR==FNR {rank[$1]=NR; next}
      $0 ~ /^#/ {next}
      {
        r = rank[$1]
        if (r=="") r=999999
        print r, $0
      }
    ' "hg38.fa.fai" "$BED_IN" > __bed.ranked

    # now sort by: contig-rank, start, end; then remove rank column
    sort -k1,1n -k3,3n -k4,4n __bed.ranked \
      | cut -f2- \
      | awk '!seen[$0]++' > __bed.sorted

    BED_ARG="-R __bed.sorted"
  fi

  REG_ARG=""
  if [ -n "!{ params.region ?: '' }" ]; then
    REG_ARG="-r !{params.region}"
  fi

  bcftools mpileup -Ob -f "!{ref_fa}" \$REG_ARG \$BED_ARG \
      -q 20 -Q 20 -a AD,DP,SP "!{cram}" -o mpileup.bcf 2> mpileup.log || {
    echo "--- mpileup.log (head) ---" >&2
    sed -n '1,200p' mpileup.log >&2 || true
    exit 1
  }

  [ -s mpileup.bcf ] || { echo "ERROR: mpileup.bcf not created" >&2; sed -n '1,200p' mpileup.log >&2 || true; exit 1; }

  # call -> sort -> bgzip VCF  (FIX: use sid-based name)
  bcftools call -mv -Ou mpileup.bcf 2> call.log \
    | bcftools sort -Oz -o "!{sid}.vcf.gz" -T . 2>> call.log || {
      echo "--- call.log (head) ---" >&2
      sed -n '1,200p' call.log >&2 || true
      exit 1
    }

  # now index works
  bcftools index -t "!{sid}.vcf.gz" 2> index.log || {
    echo "--- index.log (head) ---" >&2
    sed -n '1,200p' index.log >&2 || true
    exit 1
  }
  '''
}

process CALL_SOMATIC_MUTECT2 {
  tag { "${tsid}_vs_${nsid}" }
  container 'broadinstitute/gatk:4.5.0.0'
  cpus 6
  memory '16 GB'
  time '12h'
  stageInMode 'copy'
  publishDir "${params.outdir}/vcf_somatic", mode: 'copy', overwrite: true

  input:
    tuple val(tsid), path(t_cram), path(t_crai), val(nsid), path(n_cram), path(n_crai)
    tuple path(ref_fa), path(ref_fai)
    path ref_dict

    // stageAs prevents ".empty.tsv" collisions when these are placeholders
    tuple path(pon_vcf,  stageAs: 'pon.vcf.gz'),
          path(pon_tbi,  stageAs: 'pon.vcf.gz.tbi')

    tuple path(germ_vcf, stageAs: 'germline_resource.vcf.gz'),
          path(germ_tbi, stageAs: 'germline_resource.vcf.gz.tbi')

    path(intervals_bed,  stageAs: 'somatic_intervals.bed')


  output:
    tuple val(tsid), path("${tsid}.somatic.vcf.gz"), path("${tsid}.somatic.vcf.gz.tbi"), emit: vcf

  shell:
  '''
  set -euo pipefail

  [ -s "!{t_cram}.crai" ] || ln -sf "!{t_crai}" "!{t_cram}.crai"
  [ -s "!{n_cram}.crai" ] || ln -sf "!{n_crai}" "!{n_cram}.crai"
  [ -s "!{ref_fa}.fai" ]  || ln -sf "!{ref_fai}" "!{ref_fa}.fai"
  [ -s "!{ref_fa}.dict" ] || ln -sf "!{ref_dict}" "!{ref_fa}.dict"

  PON_ARG=()
  if [ -s "!{pon_vcf}" ] && [[ "!{pon_vcf}" == *.vcf.gz ]]; then
    [ -s "!{pon_tbi}" ] && ln -sf "!{pon_tbi}" "!{pon_vcf}.tbi" || true
    PON_ARG=(--panel-of-normals "!{pon_vcf}")
  fi

  # --- HARD REQUIRED germline resource (gnomAD) ---
  [ -s "!{germ_vcf}" ] || { echo "[ERROR] Missing germline VCF: !{germ_vcf}" >&2; exit 1; }
  [ -s "!{germ_tbi}" ] || { echo "[ERROR] Missing germline TBI: !{germ_tbi}" >&2; exit 1; }

  # ensure expected index naming exists
  ln -sf "!{germ_tbi}" "!{germ_vcf}.tbi" || true

  GERM_ARG=(--germline-resource "!{germ_vcf}")


  INT_ARG=()
  if [ -s "!{intervals_bed}" ]; then
    INT_ARG=(-L "!{intervals_bed}")
  fi

  gatk Mutect2 \
    -R "!{ref_fa}" \
    -I "!{t_cram}" -tumor "!{tsid}" \
    -I "!{n_cram}" -normal "!{nsid}" \
    "${PON_ARG[@]}" "${GERM_ARG[@]}" "${INT_ARG[@]}" \
    -O "!{tsid}.unfiltered.vcf.gz" 2> mutect2.log \
  || { echo '--- Mutect2 log ---' >&2; sed -n '1,200p' mutect2.log >&2; exit 1; }


  gatk FilterMutectCalls \
    -R "!{ref_fa}" \
    -V "!{tsid}.unfiltered.vcf.gz" \
    -O "!{tsid}.somatic.vcf.gz" 2> filter.log \
  || { echo '--- FilterMutectCalls log ---' >&2; sed -n '1,200p' filter.log >&2; exit 1; }

  gatk IndexFeatureFile -I "!{tsid}.somatic.vcf.gz"
  '''
}



/* ---------- Annotation / filtering ---------- */
process VEP_ANNOTATE {
  tag { "vep:${sid}" }
  container "${IMG_VEP}"
  cpus 6
  memory '24 GB'
  time '6h'
  publishDir "${params.outdir}/vep", mode: 'copy', overwrite: true
  input:
    tuple val(sid), path(vcf_gz)
  output:
    tuple val(sid), path("${sid}.vep.vcf.gz"), path("${sid}.vep.vcf.gz.tbi")
  when:
    params.vep_cache != null
  shell:
  """
  set -euo pipefail

  CACHE_DIR="!{params.vep_cache}"
  SPECIES="homo_sapiens"
  ASSEMBLY="!{params.vep_assembly}"
  VER="!{params.vep_cache_version}"


  exp="\${CACHE_DIR}/\${SPECIES}/\${VER}_\${ASSEMBLY}"
  if [ ! -d "\$exp" ]; then
    echo "ERROR: VEP cache dir not found at: \$exp" >&2
    ls -l "\${CACHE_DIR}/\${SPECIES}" || true
    exit 2
  fi

  EXTRA="--cache_version \$VER"

  FASTA="!{ params.vep_fasta ?: '' }"
  if [ -n "\$FASTA" ] && [ "\$FASTA" != "null" ] && [ -e "\$FASTA" ]; then
    EXTRA="\$EXTRA --fasta \$FASTA"
  fi

  vep \
    --offline --cache --dir_cache "\$CACHE_DIR" \
    --assembly "\$ASSEMBLY" --species "\$SPECIES" \
    --format vcf --vcf --compress_output bgzip \
    --everything --fork !{task.cpus} \
    \$EXTRA \
    --input_file "!{vcf_gz}" \
    --output_file "!{sid}.vep.vcf.gz"

  tabix -p vcf "!{sid}.vep.vcf.gz"
  """
}

process FILTER_FOR_HTML {
  tag { sid }
  container "${IMG_BCFTOOLS}"
  cpus 4
  memory '8 GB'
  time '2h'
  publishDir "${params.outdir}/vep", mode: 'copy', overwrite: true
  input:
    tuple val(sid), path(vcf_gz), path(tbi)
  output:
    tuple val(sid), path("${sid}.reportready.vcf.gz"), path("${sid}.reportready.vcf.gz.tbi")
  shell:
  '''
  set -euo pipefail

  IN="!{vcf_gz}"
  SID="!{sid}"
  OUT="${SID}.reportready.vcf.gz"

  USE_CSQ="!{ params.report_use_csq ? '1' : '0' }"
  KEEP_IMPACT="!{params.report_keep_impact ?: 'HIGH,MODERATE'}"
  CSQ_REGEX="!{params.report_csq_regex ?: 'missense_variant|frameshift_variant|stop_gained|start_lost|stop_lost|splice_(acceptor|donor)_variant|inframe_(insertion|deletion)'}"
  HOT_KEYS="!{ params.msk_hotspot_keys ?: '' }"
  USE_HOTSPOTS="!{ params.use_hotspots ? '1' : '0' }"

  bcftools view -s "$SID" "$IN" -Ou \
  | bcftools +fill-tags -Ou -- -t AC,AN,AF \
  | bcftools view -f PASS -i 'AC>0 && (FMT/DP>=8 || INFO/DP>=8 || (FMT/AD[0:0] + FMT/AD[0:1])>=8)' \
  | bcftools view -Oz -o pre.vcf.gz
  tabix -f -p vcf pre.vcf.gz

  if [ "$(bcftools query -l pre.vcf.gz | wc -l)" -ne 1 ]; then
    echo "$SID" > samples.txt
    bcftools reheader -s samples.txt pre.vcf.gz -Oz -o pre.tmp.vcf.gz
    mv pre.tmp.vcf.gz pre.vcf.gz
    tabix -f -p vcf pre.vcf.gz
  fi

  if [ "$USE_CSQ" = "1" ]; then
    bcftools +split-vep -d -f '%CHROM\t%POS\t%REF\t%ALT\t%IMPACT\t%Consequence\n' pre.vcf.gz > __ann.tsv
    # KEEP_IMPACT を本当に使う AWK（小修正）
    awk -F'\t' -v OFS='\t' -v keep_imp="$KEEP_IMPACT" '
      BEGIN{
        n=split(keep_imp, a, /,/);
        for(i=1;i<=n;i++){ gsub(/^ *| *$/,"",a[i]); IMP[a[i]]=1 }
      }
      { key=$1"\t"$2"\t"$3"\t"$4; if ($5 in IMP) SEL[key]=1 }
      END{ for(k in SEL) print k }
    ' __ann.tsv | sort -u > __keep.keys

    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' pre.vcf.gz | sort -u > __all.keys
    grep -Fxf __keep.keys __all.keys > __sel.keys || true

    if [ -s __sel.keys ]; then
      bcftools view -T __sel.keys -Oz -o "$OUT" pre.vcf.gz
      tabix -f -p vcf "$OUT"
    else
      cp -f pre.vcf.gz "$OUT"
      tabix -f -p vcf "$OUT"
    fi
  else
    cp -f pre.vcf.gz "$OUT"
    tabix -f -p vcf "$OUT"
  fi

  if [ "$USE_HOTSPOTS" = "1" ] && [ -n "$HOT_KEYS" ] && [ -s "$HOT_KEYS" ]; then
    TMP=$(mktemp -d __hot.XXXXXX)
    bcftools +split-vep -d -f '%CHROM|%POS|%REF|%ALT|%SYMBOL|%HGVSp\n' "$OUT" > "$TMP/ann.tsv"
    python3 - "$HOT_KEYS" "$TMP/ann.tsv" > "$TMP/msk_hits.coords" <<'PY'
import sys, csv
keys=set()
with open(sys.argv[1], newline='') as f:
    r=csv.reader(f, delimiter='\t')
    for row in r:
        if len(row) >= 2:
            keys.add((row[0], row[1]))
hits=set()
with open(sys.argv[2], newline='') as f:
    r=csv.reader(f, delimiter='|')
    for row in r:
        if len(row) < 6: continue
        chrom, pos, ref, alt, gene, prot = row[:6]
        if ':p.' in prot: prot = 'p.' + prot.split(':p.', 1)[1]
        if not prot.startswith('p.'): continue
        if (gene, prot) in keys:
            hits.add((chrom, int(pos), ref, alt, f"{gene}:{prot}"))
for chrom, pos, ref, alt, key in sorted(hits, key=lambda x: (x[0], x[1], x[2], x[3], x[4])):
    print(f"{chrom}\t{pos}\t{ref}\t{alt}\t{key}")
PY
    bgzip -c "$TMP/msk_hits.coords" > "$TMP/msk_hits.coords.gz"
    tabix -s1 -b2 -e2 "$TMP/msk_hits.coords.gz"
    printf '##INFO=<ID=MSK_HOTSPOT_KEY,Number=.,Type=String,Description="Matched MSK Cancer Hotspots (GENE:p.*)">\n' > "$TMP/msk.header"
    bcftools annotate \
      -a "$TMP/msk_hits.coords.gz" \
      -c CHROM,POS,REF,ALT,INFO/MSK_HOTSPOT_KEY \
      -h "$TMP/msk.header" \
      -Oz -o "${SID}.reportready.tmp.vcf.gz" "$OUT"
    tabix -f -p vcf "${SID}.reportready.tmp.vcf.gz"
    mv -f "${SID}.reportready.tmp.vcf.gz" "$OUT"
    mv -f "${SID}.reportready.tmp.vcf.gz.tbi" "$OUT.tbi"
  fi

  echo "reportready_rows=$(bcftools view -H "$OUT" | wc -l || true)" >&2
  echo "hotspots=$(bcftools view -H -i 'INFO/MSK_HOTSPOT_KEY!=""' "$OUT" | wc -l || true)" >&2
  '''
}

process ANNOTATE_MSK_HOTSPOTS {
  tag { sid }
  container "${IMG_BCFTOOLS}"
  cpus 2
  memory '2 GB'
  time '1h'
  publishDir "${params.outdir}/vep", mode: 'copy', overwrite: true

  input:
    tuple val(sid), path(vcf_gz), path(tbi), path(keys)

  output:
    tuple val(sid), path("${sid}.hotspot.vcf.gz"), path("${sid}.hotspot.vcf.gz.tbi")

  shell:
  """
  set -euo pipefail

  IN="!{vcf_gz}"
  OUT="!{sid}.hotspot.vcf.gz"
  KEYS="!{keys}"

  [[ -s "\$KEYS" ]] || { echo "[hotspot] Missing/empty keys: \$KEYS" >&2; exit 2; }

  # ---- Normalize keys -> __keys.tsv (GENE \\t p.X...) ----
  tr -d '\\015' < "\$KEYS" \
  | awk -F'\\t| +' 'BEGIN{OFS="\\t"}
      NF>=2 {
        g=toupper(${'$'}1); p=${'$'}2
        sub(/^.*:p\\./,"p.",p)
        gsub(/Ala/,"A",p); gsub(/Arg/,"R",p); gsub(/Asn/,"N",p); gsub(/Asp/,"D",p)
        gsub(/Cys/,"C",p); gsub(/Gln/,"Q",p); gsub(/Glu/,"E",p); gsub(/Gly/,"G",p)
        gsub(/His/,"H",p); gsub(/Ile/,"I",p); gsub(/Leu/,"L",p); gsub(/Lys/,"K",p)
        gsub(/Met/,"M",p); gsub(/Phe/,"F",p); gsub(/Pro/,"P",p); gsub(/Ser/,"S",p)
        gsub(/Thr/,"T",p); gsub(/Trp/,"W",p); gsub(/Tyr/,"Y",p); gsub(/Val/,"V",p)
        gsub(/Stop|Ter/,"*",p); gsub(/Sec/,"U",p)
        if (g!="" && p ~ /^p\\./) print g,p
      }' \
  | LC_ALL=C sort -u > __keys.tsv

  # ---- Split VEP CSQ fields ----
  bcftools +split-vep -d -f "%CHROM\\t%POS\\t%REF\\t%ALT\\t%SYMBOL\\t%HGVSp\\n" "\$IN" > __ann.raw.tsv

  # ---- Aggregate by allele (CHROM:POS:REF:ALT) ----
  awk -F'\\t' -v OFS='\\t' '
    NR==FNR { K[${'$'}1 FS ${'$'}2]=1; next }
    function normprot(s){
      sub(/^.*:p\\./,"p.",s)
      gsub(/Ala/,"A",s); gsub(/Arg/,"R",s); gsub(/Asn/,"N",s); gsub(/Asp/,"D",s)
      gsub(/Cys/,"C",s); gsub(/Gln/,"Q",s); gsub(/Glu/,"E",s); gsub(/Gly/,"G",s)
      gsub(/His/,"H",s); gsub(/Ile/,"I",s); gsub(/Leu/,"L",s); gsub(/Lys/,"K",s)
      gsub(/Met/,"M",s); gsub(/Phe/,"F",s); gsub(/Pro/,"P",s); gsub(/Ser/,"S",s)
      gsub(/Thr/,"T",s); gsub(/Trp/,"W",s); gsub(/Tyr/,"Y",s); gsub(/Val/,"V",s)
      gsub(/Stop|Ter/,"*",s); gsub(/Sec/,"U",s)
      return s
    }
    {
      chr=${'$'}1; pos=${'$'}2; ref=${'$'}3; alt=${'$'}4; gene=toupper(${ '$'}5); prot=normprot(${ '$'}6)
      if (gene!="" && prot ~ /^p\\./) {
        if ((gene FS prot) in K) {
          key = gene ":" prot
          loc = chr FS pos FS ref FS alt
          if (A[loc] == "") A[loc] = key
          else if (index("," A[loc] ",", "," key ",")==0) A[loc] = A[loc] "," key
        }
      }
    }
    END { for (loc in A) { split(loc,f,FS); print f[1],f[2],f[3],f[4],A[loc] } }
  ' __keys.tsv __ann.raw.tsv > __anno.tsv

  # ---- Header describing the INFO tag ----
  printf '##INFO=<ID=MSK_HOTSPOT_KEY,Number=.,Type=String,Description="Matched MSK Cancer Hotspots (GENE:p.*)">\\n' > __hdr.txt

  # ---- Annotate ----
  if [ -s __anno.tsv ]; then
    LC_ALL=C sort -k1,1 -k2,2n -k3,3 -k4,4 __anno.tsv > __anno.sorted.tsv
    bgzip -f -c __anno.sorted.tsv > __anno.tsv.gz
    tabix -f -s1 -b2 -e2 __anno.tsv.gz
    bcftools annotate -a __anno.tsv.gz -c CHROM,POS,REF,ALT,INFO/MSK_HOTSPOT_KEY -h __hdr.txt "\$IN" -Oz -o "\$OUT"
  else
    bcftools annotate -h __hdr.txt "\$IN" -Oz -o "\$OUT"
  fi

  tabix -f -p vcf "\$OUT"

  echo "[hotspot] n_hits=\$(bcftools view -H -i 'INFO/MSK_HOTSPOT_KEY!=""' \"\$OUT\" | wc -l || true)" >&2
  echo "[hotspot] unique_keys=\$(bcftools query -i 'INFO/MSK_HOTSPOT_KEY!=""' -f '%INFO/MSK_HOTSPOT_KEY\n' \"\$OUT\" | tr ',' '\n' | LC_ALL=C sort -u | wc -l || true)" >&2
  """
}










/* ---------- Panel QC & aggregation ---------- */
process PANEL_QC_MOSDEPTH {
  tag { "panelqc_${sid}" }
  stageInMode 'link'
  cpus 6
  memory '24 GB'
  time '6h'
  publishDir "${params.outdir}/qc_panel", mode:'copy', overwrite:true
  input:
    tuple val(sid), path(aln), path(aln_index)
    path  panel_bed
    tuple path(ref_fa), path(ref_fai)
  output:
    tuple val(sid), path("panel_qc.mosdepth.summary.txt"), emit: summary
    tuple val(sid), path("panel_qc.regions.bed.gz")      , emit: regions
    tuple val(sid), path("panel_qc.thresholds.bed.gz")   , emit: thresholds
    tuple val(sid), path("${sid}.panel_qc.json")         , emit: panel_json
  shell:
  """
  set -euo pipefail
  MOS="!{ params.mosdepth_bin ?: 'mosdepth' }"
  SAM="!{ params.samtools_bin ?: 'samtools' }"

  FA_ARG=""
  case "!{aln}" in
    *.cram) FA_ARG="-f !{ref_fa}" ;;
  esac

  [ -s "!{ref_fa}.fai" ] || "\$SAM" faidx "!{ref_fa}"

  case "!{aln}" in
    *.cram) "\$SAM" index -@ !{task.cpus} -c "!{aln}" ;;
    *.bam)  "\$SAM" index -@ !{task.cpus}     "!{aln}" ;;
  esac

  "\$MOS" -t !{task.cpus} \${FA_ARG} --by "!{panel_bed}" --thresholds 20,50,100 panel_qc "!{aln}"
  [ -f panel_qc.mosdepth.summary.txt ] || mv panel_qc*.summary.txt panel_qc.mosdepth.summary.txt

  python3 - << 'PY'
import json, glob
thr = int(!{params.panel_min_cov})
out = {"threshold": thr, "mosdepth": {}}
for s in glob.glob("panel_qc*.summary.txt"):
    for line in open(s):
        p = line.rstrip().split("\t")
        if len(p) >= 4 and p[0] == "total_regions": out["mosdepth"]["regions"] = int(p[1])
        if len(p) >= 7 and p[0] == "total_bases":   out["mosdepth"]["bases"]   = int(p[1])
open("panel_qc.json","w").write(json.dumps(out, indent=2))
PY

  mv -f panel_qc.json "!{sid}.panel_qc.json"
  """
}

process QC_AGGREGATE {
  tag { sid }
  cpus 4
  memory '4 GB'
  time '2h'
  stageInMode 'copy'
  publishDir "${params.outdir}/qc_panel", mode:'copy', overwrite:true

  input:
    tuple val(sid), path(summary_txt), path(markdup_metrics), path(fastp_json), path(cram), path(crai), path(thresholds_bedgz), path(regions_bedgz)
    path panel_bed
    tuple path(ref_fa), path(ref_fai)

  output:
    tuple val(sid), path("${sid}.panel_qc.json")

  shell:
  """
  set -euo pipefail

  export SAM="!{ params.samtools_bin ?: 'samtools' }"
  export CPUS="!{ task.cpus }"

  # Avoid hard-link pitfalls on WSL2
  if [ ! -e "!{cram}.crai" ] && [ -e "!{crai}" ]; then cp -f "!{crai}" "!{cram}.crai"; fi
  if [ ! -e "!{ref_fa}.fai" ] && [ -e "!{ref_fai}" ]; then cp -f "!{ref_fai}" "!{ref_fa}.fai"; fi

  python3 - <<'PY'
import gzip, json, statistics as stats, os, re, subprocess

summary = "!{summary_txt}"
markdup = "!{markdup_metrics}"
fastp   = "!{fastp_json}"
thresh  = "!{thresholds_bedgz}"
regions = "!{regions_bedgz}"
sam     = os.environ.get("SAM","samtools")
cram    = "!{cram}"
panel   = "!{panel_bed}"
cpus    = os.environ.get("CPUS","1")  # read from env (string OK for subprocess)

def _to_int(x):
    try: return int(x)
    except: return None

def _to_float(x):
    try: return float(x)
    except: return None

total_mapped = 0
try:
    out = subprocess.check_output([sam, "idxstats", cram], text=True, stderr=subprocess.DEVNULL)
    for line in out.splitlines():
        if line.strip():
            p = line.split('\\t')
            if len(p) >= 3:
                total_mapped += (_to_int(p[2]) or 0)
except Exception:
    pass

ont = 0
try:
    ont = int(subprocess.check_output([sam, "view","-@",cpus,"-c","-F","260","-L",panel,cram],
                                      text=True, stderr=subprocess.DEVNULL).strip() or "0")
except Exception:
    pass

def percent(x, total): return round((x/total)*100.0, 1) if total and x is not None else None
on_target_rate = percent(ont, total_mapped)

total_len = 0
sum_cov   = 0.0
region_means = []
with gzip.open(regions, "rt") as f:
    for line in f:
        if not line.strip() or line.startswith("#"): continue
        p = re.split(r"\\s+", line.strip())
        if len(p) < 3: continue
        start = _to_int(p[1]); end = _to_int(p[2])
        if start is None or end is None: continue
        L = max(0, end - start)
        nums = [_to_float(x) for x in p[3:] if re.fullmatch(r"[-+]?\\d*\\.?\\d+(?:[eE][-+]?\\d+)?", x)]
        mean = nums[-1] if nums else None
        if L > 0 and mean is not None:
            total_len += L
            sum_cov   += mean * L
            region_means.append(mean)

mean_cov   = round(sum_cov/total_len, 1) if total_len else None
median_cov = round(stats.median(region_means), 1) if region_means else None

bases_total = 0
ge = {20:0, 50:0, 100:0}
try:
    with gzip.open(thresh, "rt") as f:
        for line in f:
            if not line.strip() or line.startswith("#"): continue
            p = re.split(r"\\s+", line.strip())
            if len(p) < 3: continue
            s = _to_int(p[1]); e = _to_int(p[2]); L = max(0, e - s)
            vals = [_to_int(x) for x in p[3:] if re.fullmatch(r"-?\\d+", x)]
            if L > 0:
                bases_total += L
                for i, thr in enumerate((20, 50, 100)):
                    if i < len(vals) and vals[i] is not None:
                        ge[thr] += max(0, vals[i])
except Exception:
    pass

def pct(x): return round(100.0 * x / bases_total, 1) if bases_total and x is not None else None
ge20, ge50, ge100 = pct(ge[20]), pct(ge[50]), pct(ge[100])

dup_pct = None
try:
    def split_metrics(s): return s.strip().split('\\t') if '\\t' in s else re.split(r"\\s+", s.strip())
    with open(markdup) as fh:
        for line in fh:
            if line.startswith("LIBRARY"):
                header = split_metrics(line)
                idx = header.index("PERCENT_DUPLICATION")
                row = split_metrics(next(fh))
                dup_pct = round(float(row[idx]) * 100.0, 1)
                break
except Exception:
    pass

q30 = None
try:
    j = json.load(open(fastp))
    q30 = round(float(j["summary"]["after_filtering"]["q30_rate"]) * 100.0, 1)
except Exception:
    pass

out = {
  "mean_coverage": mean_cov,
  "median_coverage": median_cov,
  "ge20_pct": ge20,
  "ge50_pct": ge50,
  "ge100_pct": ge100,
  "on_target_rate_pct": on_target_rate,
  "contamination_pct": None,
  "duplication_pct": dup_pct,
  "insert_size_median_bp": None,
  "q30_pct": q30
}
open("panel_qc.json","w").write(json.dumps(out, indent=2))
PY

  mv -f panel_qc.json "!{sid}.panel_qc.json"
  """
}

/* ---------- Report JSON shaping ---------- */
process MAKE_SLIM_JSON {
  tag { sid }
  container "${IMG_PY}"
  cpus 1
  memory '24 GB'
  time '30m'
  input:
    tuple val(sid), path(full_json)
  output:
    tuple val(sid), path("${sid}.report.slim.json"), path("${sid}.report.full.json")
  shell:
  """
  set -euo pipefail
  python3 - << 'PY'
import json, shutil
sid = "!{sid}"
N = int("!{params.html_max_rows}")
drop_pharm = True if "!{params.drop_pharmcat_in_html}" == "true" else False
data = json.load(open("!{full_json}"))
rows = data.get("rows", [])
data["rows"] = rows[:N]
summ = data.setdefault("summary", {})
tops = summ.get("top_variants", [])
summ["top_variants"] = tops[:min(len(data["rows"]), len(tops))]
if drop_pharm:
    data["pharmcat"] = None
open(f"{sid}.report.slim.json","w").write(json.dumps(data, separators=(',',':')))
shutil.copyfile("!{full_json}", f"{sid}.report.full.json")
PY
  """
}

process MERGE_TN_FOR_HTML {
  tag { tsid }
  container 'python:3.11'
  cpus 1
  memory '1 GB'
  time '5m'
  input:
    tuple val(tsid), path(tumor_slim), path(tumor_full),
         path(tumor_qc,  stageAs: "tumor_qc.json"),
         val(nsid), path(normal_qc, stageAs: "normal_qc.json")
  output:
    tuple val(tsid), path("${tsid}.report.slim.json"), path("${tsid}.report.full.json")
  shell:
  """
  set -euo pipefail
  python3 - <<'PY'
import json, os
tsid = "!{tsid}"; nsid = "!{nsid}"
data = json.load(open("!{tumor_slim}"))
full = json.load(open("!{tumor_full}"))
tqc  = json.load(open("tumor_qc.json"))  if os.path.getsize("tumor_qc.json")  else {}
nqc  = json.load(open("normal_qc.json")) if os.path.getsize("normal_qc.json") else {}
data.setdefault("qc", data.get("qc", {}))
data["meta"]       = dict(data.get("meta", {}), pair={"tumor_id": tsid, "normal_id": nsid})
data["qc_pair"]    = { tsid: tqc, nsid: nqc }
data["qc_compare"] = { tsid: tqc, nsid: nqc }
data["qc_samples"] = [{"id": tsid, **tqc}, {"id": nsid, **nqc}]
full.setdefault("meta", {})
full["meta"]       = dict(full["meta"], pair={"tumor_id": tsid, "normal_id": nsid})
full["qc_pair"]    = { tsid: tqc, nsid: nqc }
full["qc_compare"] = { tsid: tqc, nsid: nqc }
full["qc_samples"] = [{"id": tsid, **tqc}, {"id": nsid, **nqc}]
json.dump(data, open(f"{tsid}.report.slim.json","w"), separators=(',',':'))
json.dump(full, open(f"{tsid}.report.full.json","w"))
PY
  """
}

/* ---------- Report HTML build & injection ---------- */
process BUILD_HTML {
  tag { sample_id }
  container 'python:3.11'
  cpus 6
  memory '20 GB'
  time '30m'

  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { fn ->
    def s = fn.toString()
    if (s.startsWith('OpenCare_') && s.endsWith('_report.html')) {
      def sid = s.substring(9, s.length() - 12)
      return sid + "/" + s
    }
    if (s.startsWith('out_')) {
      int slash = s.indexOf('/')
      if (slash > 4) {
        def sid  = s.substring(4, slash)
        def rest = s.substring(slash + 1)
        return sid + "/" + rest
      }
    }
    return s
  }

  input:
    tuple val(sample_id), path(report_json_slim), path(report_json_full)
    path  build_script
    path  pathway_db
    path  gene_domains
    path  civic_offline_json

  output:
    path "OpenCare_*_report.html", emit: html
    path "out_*/OpenCare_*_report.full.json", emit: full_json
    path "OpenCARE2_report.full.json", emit: full_json_copy


  shell:
  '''
  set -euo pipefail
  safe_id="!{ sample_id.toString().replaceAll(/[^A-Za-z0-9._-]/,'_') }"
  mkdir -p "modules_${safe_id}" "out_${safe_id}"

  # ---- Build PD (pathway DB flag) only if JSON is a non-empty list/dict
  PD=()   # <-- ARRAY (not a string)
  if [ -s "!{pathway_db}" ]; then
    if python3 - "!{pathway_db}" <<'PY'
import json, sys
p = sys.argv[1]
try:
    d = json.load(open(p, 'r', encoding='utf-8'))
except Exception as e:
    print("[diag] pathway_db read error:", e); sys.exit(2)
t = type(d).__name__
n = len(d) if hasattr(d, '__len__') else 'NA'
sample = (d[:2] if isinstance(d, list) else list(d.items())[:2] if isinstance(d, dict) else '')
print(f"[diag] pathway_db type={t} len={n} sample={sample}")
ok = isinstance(d, (list, dict)) and bool(d)
sys.exit(0 if ok else 1)
PY
    then
      PD+=(--pathway-db "!{pathway_db}")
    else
      echo "[diag] pathway_db is empty/invalid → using builder's internal DB" >&2
    fi
  else
    echo "[diag] pathway_db not found or zero bytes: !{pathway_db}" >&2
  fi

  # ---- Gene domains flag (presence check only)
  GD=()   # <-- ARRAY (not a string)
  [ -s "!{gene_domains}" ] && GD+=(--gene-domains "!{gene_domains}")

  # ---- Environment knobs
  export ENABLE_EVIDENCE='!{ params.enable_evidence ? "1" : "0" }'
  if [ -n "!{ params.oncokb_token ?: '' }" ]; then
    export ONCOKB_TOKEN="!{ params.oncokb_token }"
  fi
  if [ -s "!{civic_offline_json}" ]; then
    export CIVIC_OFFLINE="!{civic_offline_json}"
  fi

  # (debug) show what we will pass
  echo "[debug] builder args: ${PD[*]} ${GD[*]}"

  # ---- Build HTML (script writes into out_${safe_id})
  python3 "!{build_script}" \
    --report "!{report_json_slim}" \
    --modules "modules_${safe_id}" \
    --outdir  "out_${safe_id}" \
    "${PD[@]}" "${GD[@]}"

  # ---- Copy full JSON to the canonical name the HTML should link to
  cp -f "!{report_json_full}" "out_${safe_id}/OpenCare_${safe_id}_report.full.json"
  cp -f "!{report_json_full}" "OpenCARE2_report.full.json"

  # ---- Find produced HTML and normalize to OpenCare_${safe_id}_report.html
  html_src=""
  if [ -s "out_${safe_id}/OpenCare_${safe_id}_report.html" ]; then
    html_src="out_${safe_id}/OpenCare_${safe_id}_report.html"
  else
    html_src="$(ls -1t out_${safe_id}/OpenCare*_report*.html 2>/dev/null | head -n1 || true)"
    [ -z "$html_src" ] && html_src="$(ls -1t out_${safe_id}/*.html 2>/dev/null | head -n1 || true)"
  fi

  if [ -z "${html_src}" ] || [ ! -s "${html_src}" ]; then
    echo "ERROR: No HTML produced in out_${safe_id}" >&2
    echo "--- out_${safe_id} listing ---" >&2
    ls -l out_${safe_id} >&2 || true
    exit 2
  fi

  # ---- Ensure the JSON download button points at the canonical full JSON
  sed -E -i "s|const JSON_DL = ['\"][^'\"]+['\"];|const JSON_DL = 'OpenCare_${safe_id}_report.full.json';|g" "${html_src}" || true

  # ---- Write canonical copies
  cp -f "${html_src}" "out_${safe_id}/OpenCare_${safe_id}_report.html"
  cp -f "out_${safe_id}/OpenCare_${safe_id}_report.html" "OpenCare_${safe_id}_report.html"

  echo "HTML source: ${html_src}"
  echo "Wrote (canonical): out_${safe_id}/OpenCare_${safe_id}_report.html"
  echo "Wrote (root copy):  OpenCare_${safe_id}_report.html"
  '''
}



process INJECT_PAIRED_QC_HTML {
  tag { out_name }
  container 'python:3.11'
  cpus 1
  memory '1 GB'
  time '5m'
  stageInMode 'copy'

  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { fn ->
    def s = fn.toString()
    if (s.startsWith('OpenCare_') && s.endsWith('_report.html')) {
      def sid = s.substring(9, s.length() - 12)
      return sid + "/" + s
    }
    return s
  }

  input:
    tuple val(out_name), path(html_in, stageAs: 'orig.html')

  output:
    path "${out_name}", emit: html

  shell:
  $/
  set -euo pipefail
  cp -f orig.html in.html

  python3 - <<'PY'
import os

p = 'in.html'
html = open(p, 'r', encoding='utf-8', errors='ignore').read()

# NOTE: keep ONLY ONE <script> opening tag
patch = r"""<script>(function(){
  function ready(fn){ if(document.readyState!=='loading') fn(); else document.addEventListener('DOMContentLoaded', fn); }
  function ensureJSON(cb){
    if (window.reportFull || window.report){ cb(window.reportFull || window.report); return; }
    const htmlFile = (location.pathname.split('/').pop() || '');
    const fullJson = htmlFile.replace(/_report\.html$/,'_report.full.json');
    fetch(fullJson).then(r => r.ok ? r.json() : null).then(j => {
      if (j){ window.reportFull = j; cb(j); }
    }).catch(()=>{});
  }
  function fmt(x){ return (x==null || x==='') ? '—' : x; }
  function gradeLabel(cls){ return cls==='good' ? 'Good' : cls==='acc' ? 'Acceptable' : 'Concerning'; }

  function evalOnT(x){
    if (x==null) return {cls:'bad', label:'Concerning', reason:'Not reported'};
    if (x>=80)   return {cls:'good', label:'Good',       reason:'Efficient capture'};
    if (x>=65)   return {cls:'acc',  label:'Acceptable', reason:'Acceptable capture; some off-target reads'};
    return        {cls:'bad', label:'Concerning', reason:'Poor capture; coverage uneven'};
  }
  function evalDup(x){
    if (x==null) return {cls:'bad', label:'Concerning', reason:'Not reported'};
    if (x<20)    return {cls:'good', label:'Good',       reason:'Low duplication'};
    if (x<40)    return {cls:'acc',  label:'Acceptable', reason:'Moderate duplication'};
    return        {cls:'bad', label:'Concerning', reason:'High duplication; reduces unique depth'};
  }
  function evalQ30(x){
    if (x==null) return {cls:'bad', label:'Concerning', reason:'Not reported'};
    if (x>=90)   return {cls:'good', label:'Good',       reason:'Base quality high'};
    if (x>=85)   return {cls:'acc',  label:'Acceptable', reason:'Moderate base quality'};
    return        {cls:'bad', label:'Concerning', reason:'Sub-optimal base quality; filters impacted'};
  }
  function evalCovMean(x){
    if (x==null) return {cls:'bad', label:'Concerning', reason:'Not reported'};
    if (x>=50)   return {cls:'good', label:'Good',       reason:'Ample depth for somatic calling'};
    if (x>=30)   return {cls:'acc',  label:'Acceptable', reason:'Adequate; borderline for low-VAF'};
    return        {cls:'bad', label:'Concerning',        reason:'Insufficient depth; low-VAF sensitivity reduced'};
  }
  function evalCovMedian(x){
    if (x==null) return {cls:'bad', label:'Concerning', reason:'Not reported'};
    if (x>=50)   return {cls:'good', label:'Good',       reason:'Median depth is strong'};
    if (x>=30)   return {cls:'acc',  label:'Acceptable', reason:'Median depth borderline'};
    return        {cls:'bad', label:'Concerning',        reason:'Median depth low'};
  }

  const evalGE20  = x => (x==null) ? {cls:'bad',label:'Concerning',reason:'Not reported'}
    : x>=80 ? {cls:'good',label:'Good',reason:'Most targets ≥20×'}
    : x>=60 ? {cls:'acc', label:'Acceptable',reason:'Moderate breadth ≥20×'}
            : {cls:'bad', label:'Concerning',reason:'Large under-covered fraction <20×'};

  const evalGE50  = x => (x==null) ? {cls:'bad',label:'Concerning',reason:'Not reported'}
    : x>=60 ? {cls:'good',label:'Good',reason:'Good high-depth fraction (≥50×)'}
    : x>=40 ? {cls:'acc', label:'Acceptable',reason:'Moderate ≥50× fraction'}
            : {cls:'bad', label:'Concerning',reason:'Limited ≥50× coverage'};

  const evalGE100 = x => (x==null) ? {cls:'bad',label:'Concerning',reason:'Not reported'}
    : x>=30 ? {cls:'good',label:'Good',reason:'Many regions reach ≥100×'}
    : x>=15 ? {cls:'acc', label:'Acceptable',reason:'Some regions reach ≥100×'}
            : {cls:'bad', label:'Concerning',reason:'Deep (≥100×) coverage limited'};

  function capsuleTile(title, value, ev, unit){
    var d=document.createElement('div');
    d.className='qcbox ' + (ev.cls || '');
    var v = (value==null? '—' : (unit==='%'? Math.round(value) + '%' : value));
    d.innerHTML =
      '<div class="grade">' + (ev.label || '') + '</div>' +
      '<div><strong>' + title + ':</strong> ' + v + '</div>' +
      '<div><small class="muted">' + (ev.reason || '') + '</small></div>';
    return d;
  }

  function render(r){
    if(!r || !r.meta || !r.meta.pair || !r.qc_pair) return;

    var qs=document.getElementById('qcSummary'); if(qs) qs.style.display='none';
    var chips=document.getElementById('qcChips'); if(chips) chips.style.display='none';

    var tId=r.meta.pair.tumor_id, nId=r.meta.pair.normal_id;
    var t=r.qc_pair[tId]||{}, n=r.qc_pair[nId]||{};
    var panel=document.getElementById('qcPanel');

    if(panel){
      panel.innerHTML='';
      var wrap=document.createElement('div'); wrap.className='qc-2col';

      function block(name,q){
        var box=document.createElement('div'); box.className='qc-sample';
        var h=document.createElement('h4'); h.textContent=name; box.appendChild(h);
        var grid=document.createElement('div'); grid.className='qcboxes';

        grid.appendChild(capsuleTile('Mean coverage',     q.mean_coverage,        evalCovMean(q.mean_coverage),           '×'));
        grid.appendChild(capsuleTile('Median coverage',   q.median_coverage,      evalCovMedian(q.median_coverage),      '×'));
        grid.appendChild(capsuleTile('≥20× (%)',          q.ge20_pct,             evalGE20(q.ge20_pct),                  '%'));
        grid.appendChild(capsuleTile('≥50× (%)',          q.ge50_pct,             evalGE50(q.ge50_pct),                  '%'));
        grid.appendChild(capsuleTile('≥100× (%)',         q.ge100_pct,            evalGE100(q.ge100_pct),                '%'));
        grid.appendChild(capsuleTile('On-target (%)',     q.on_target_rate_pct,   evalOnT(q.on_target_rate_pct),         '%'));
        grid.appendChild(capsuleTile('Duplication (%)',   q.duplication_pct,      evalDup(q.duplication_pct),            '%'));
        grid.appendChild(capsuleTile('Q30 (%)',           q.q30_pct,              evalQ30(q.q30_pct),                    '%'));

        box.appendChild(grid);
        return box;
      }

      wrap.appendChild(block(tId, t));
      wrap.appendChild(block(nId, n));
      panel.appendChild(wrap);

      var qcBlock=document.getElementById('qcBlock');
      if(qcBlock) qcBlock.style.display='';
    }

    if(!document.getElementById('qc2colStyle')){
      var style=document.createElement('style'); style.id='qc2colStyle';
      style.textContent=[
        '.qc-2col{display:grid;grid-template-columns:1fr 1fr;gap:12px}',
        '.qc-sample h4{margin:0 0 6px}'
      ].join('\\n');
      document.head.appendChild(style);
    }
  }

  ready(function(){
    if (window.reportFull || window.report) render(window.reportFull || window.report);
    else ensureJSON(render);
  });
})();</script>
"""

if '</body>' in html:
  html = html.replace('</body>', patch + '\n</body>')
else:
  html = html + '\n' + patch

open(p, 'w', encoding='utf-8').write(html)
PY

  mv -f in.html "!{out_name}"
  /$
}




/* ---------- PharmCAT (optional) ---------- */
process PHARMCAT {
  tag { "pharmcat:${sid}" }
  container "pgkb/pharmcat:3.0.1"
  cpus 4
  memory '8 GB'
  time '3h'
  publishDir "${params.outdir}/pharmcat", mode: 'copy', overwrite: true
  input:
    tuple val(sid), path(vcf)
  output:
    tuple val(sid), path("${sid}.pharmcat.report.html"), emit: pharmcat_html
    tuple val(sid), path("${sid}.pharmcat.report.json"), emit: pharmcat_json
  shell:
  """
  set -euo pipefail
  IN="input.vcf"
  case "!{vcf}" in
    *.gz|*.bgz)  gunzip -c "!{vcf}" > "\$IN" ;;
    *)           cp "!{vcf}" "\$IN" ;;
  esac
  pharmcat_pipeline "\$IN" -bf "!{sid}" -reporterHtml -reporterJson
  mv "!{sid}.report.html" "!{sid}.pharmcat.report.html"
  mv "!{sid}.report.json" "!{sid}.pharmcat.report.json" || \
    cp "!{sid}".*.json "!{sid}.pharmcat.report.json" 2>/dev/null || true
  """
}



process CIVIC_OFFLINE_ENRICH {
  tag { sid }
  container "${IMG_PY}"
  publishDir "${params.outdir}/json", mode: 'copy', overwrite: true

  input:
    tuple val(sid), path(report_json)
    path civic_json

  output:
    tuple val(sid), path("${sid}.report.civic.json"), emit: civic_report

  when:
    params.enable_civic

  shell:
  $/
  set -euo pipefail

  python3 - <<'PY'
import json, re, sys
from pathlib import Path

sid         = "!{sid}"
report_json = Path("!{report_json}")
civic_json  = Path("!{civic_json}")
out_json    = Path(f"{sid}.report.civic.json")

R       = json.loads(report_json.read_text(encoding="utf-8"))
IDX_RAW = json.loads(civic_json.read_text(encoding="utf-8"))

AA3 = {
    "Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G",
    "His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S",
    "Thr":"T","Trp":"W","Tyr":"Y","Val":"V","Ter":"*","Stop":"*","Sec":"U","Pyl":"O"
}
AA3_RE = re.compile(r'\\b(?:' + '|'.join(map(re.escape, AA3.keys())) + r')\\b', re.I)

def _tri_to_one_token(tok: str) -> str:
    return AA3_RE.sub(lambda m: AA3.get(m.group(0).capitalize(), m.group(0)), tok)

def _aa_key(hgvsp: str) -> str:
    # Return compact AA key like 'R175H' from p.Arg175His / p.R175H / ENSP:...
    s = str(hgvsp or "").strip().split(":")[-1]
    if not s:
        return ""
    s = _tri_to_one_token(s).replace(" ", "")
    m = re.search(r'(?i)\\bp\\.([A-Z\\*])(\\d+)([A-Z\\*])\\b', s)
    if not m:
        return ""
    f, pos, t = m.groups()
    return f"{f.upper()}{pos}{t.upper()}"

def _candidate_keys(hgvsp: str):
    out = set()
    s = str(hgvsp or "").strip()
    if not s:
        return []
    s = s.replace(" ", "")
    last = s.split(':')[-1]
    out.add(last)

    m3 = re.search(r'(?i)\\bp\\.([A-Za-z]{3})(\\d+)([A-Za-z]{3}|\\*)\\b', last)
    if m3:
        f3,pos,t3 = m3.groups()
        f1 = AA3.get(f3.capitalize(), f3[0]).upper()
        t1 = '*' if t3 == '*' else AA3.get(t3.capitalize(), t3[0]).upper()
        out.add(f"p.{f1}{pos}{t1}")

    m1 = re.search(r'(?i)\\bp\\.([A-Za-z\\*])(\\d+)([A-Za-z\\*])\\b', last)
    if m1:
        f1,pos,t1 = m1.groups()
        out.add(f"p.{f1.upper()}{pos}{t1.upper()}")

    for t in list(out):
        if t.lower().startswith('p.'):
            out.add(t[2:])  # prefixless (R175H)

    if ':' in s and last:
        out.add(s)  # keep ENSP:...:p.Arg175His too

    casey = set()
    for k in out:
        casey.update({k, k.upper(), k.lower(), k.title()})
    return list(casey)

def _norm_ev(s: str) -> str:
    if not s: return ""
    return re.sub(r'(?i)\\bLEVEL[\\s_]*', '', s).strip().upper()

def _eid(url: str) -> str:
    # Normalize CIViC evidence id from URL (or fall back to canonical URL).
    if not url: return ""
    u = url.strip().rstrip('/')
    m = re.search(r'/evidence_items/(\\d+)$', u)
    return m.group(1) if m else u

LABEL_RE = re.compile(r'(FDA|EMA|\\blabel\\b|NCCN|ESMO|guideline)', re.I)
CIVIC_T1 = re.compile(r'\\b(A|LEVEL[_ ]?A|B|LEVEL[_ ]?B)\\b', re.I)
CIVIC_T2 = re.compile(r'\\b(C|LEVEL[_ ]?C|D|LEVEL[_ ]?D)\\b', re.I)

def promote_tier(row: dict):
    ts = row.get('treatments') or []
    joined = ' '.join(f"{t.get('evidence','')} {t.get('source_url','')} {t.get('source','')}" for t in ts)
    if LABEL_RE.search(joined): row['tier'] = 'T1'; return
    if any(CIVIC_T1.search(_norm_ev(t.get('evidence',''))) for t in ts): row['tier'] = 'T1'; return
    if any(CIVIC_T2.search(_norm_ev(t.get('evidence',''))) for t in ts): row['tier'] = 'T2'; return

# --------- NEW: make civic_offline format-agnostic (flat "GENE|VAR" OR nested) ----------
def _norm_var_keys(v: str):
    s = str(v or "").strip()
    if not s:
        return []
    s = s.replace(" ", "")
    s1 = _tri_to_one_token(s)
    keys = set([s, s1])

    for x in list(keys):
        keys.add(x.upper()); keys.add(x.lower())
        keys.add(x.replace("P.","p.")); keys.add(x.replace("p.","P."))
        keys.add(x.replace("P.","").replace("p.",""))
        if x.lower().startswith("p."):
            keys.add(x[2:])  # prefixless

    return [k for k in keys if k]

def _as_list(x):
    if x is None: return []
    return x if isinstance(x, list) else [x]

def build_canonical_index(raw):
    """
    Returns: IDX[GENE][VARKEY] -> list[hit_dict]
    Supports:
      - flat dict:  { "TP53|P.R175H": [...], ... }
      - nested dict:{ "TP53": { "p.R175H": [...], ... }, ... }
    """
    idx = {}
    if not isinstance(raw, dict):
        print("[civic] civic_offline is not a dict; cannot index", file=sys.stderr)
        return idx

    sample_keys = list(raw.keys())[:200]
    is_flat = any(isinstance(k, str) and "|" in k for k in sample_keys)

    if is_flat:
        for k, payload in raw.items():
            if not isinstance(k, str) or "|" not in k:
                continue
            gene, var = k.split("|", 1)
            g = gene.strip().upper()
            v = var.strip()
            if not g or not v:
                continue
            idx.setdefault(g, {})
            hits = [h for h in _as_list(payload) if isinstance(h, dict)]
            if not hits:
                continue
            for vk in _norm_var_keys(v):
                idx[g].setdefault(vk, []).extend(hits)
        print(f"[civic] index_mode=flat genes={len(idx)}", file=sys.stderr)
    else:
        for gene, sub in raw.items():
            if not isinstance(sub, dict):
                continue
            g = str(gene).strip().upper()
            if not g:
                continue
            idx.setdefault(g, {})
            for var, payload in sub.items():
                v = str(var).strip()
                hits = [h for h in _as_list(payload) if isinstance(h, dict)]
                if not v or not hits:
                    continue
                for vk in _norm_var_keys(v):
                    idx[g].setdefault(vk, []).extend(hits)
        print(f"[civic] index_mode=nested genes={len(idx)}", file=sys.stderr)

    return idx

IDX = build_canonical_index(IDX_RAW)

def civic_hits(gene, hgvsp):
    if not gene or not hgvsp: return []
    G = str(gene).strip().upper()
    sub = IDX.get(G)
    if not isinstance(sub, dict):
        return []

    keys = _candidate_keys(hgvsp)

    out = []
    for k in keys:
        for kk in _norm_var_keys(k):
            v = sub.get(kk)
            if v:
                out.extend(v)

    # Dedup hits by evidence signature
    uniq, seen = [], set()
    for h in out:
        if not isinstance(h, dict):
            continue
        url = (h.get('source_url') or h.get('url') or h.get('link') or '').strip()
        eid = _eid(url)
        evl = _norm_ev(h.get('evidence') or h.get('evidence_level') or h.get('evidence_label') or h.get('level') or '')
        th  = (h.get('therapy') or h.get('drug') or h.get('agent') or h.get('treatment') or '').strip().lower()
        pm  = (h.get('pmid') or h.get('PMID') or '').strip()
        sig = (eid, evl, th, pm)
        if sig not in seen:
            seen.add(sig)
            uniq.append(h)
    return uniq
# --------------------------------------------------------------------------------------

CT_IGNORE = re.compile(r'^(any|pan-?cancer|all tumor|n/?a)\\b', re.I)

def _merge_comment(existing: str, addition: str) -> str:
    ex = (existing or '').strip()
    if not ex:
        return addition
    return ex if addition.lower() in ex.lower() else (ex + " | " + addition)

def add_ct_to_comment(row: dict):
    # If tier is T1, add Cancer Type(s) from treatments into clinical_comment/comment.
    tier = (row.get('tier') or '').upper()
    if not tier.startswith('T1'):
        return
    ts = row.get('treatments') or []
    cts = []
    for t in ts:
        ct = (t.get('cancer_type') or '').strip()
        if ct and not CT_IGNORE.match(ct):
            ct = re.sub(r'\\s+', ' ', ct)
            if ct not in cts:
                cts.append(ct)
    if not cts:
        return
    addition = "Cancer Type: " + ", ".join(cts)
    new_comment = _merge_comment(row.get('clinical_comment') or row.get('comment'), addition)
    row['clinical_comment'] = new_comment
    row['comment'] = new_comment

promoted = 0
for row in R.get('rows', []):
    g = row.get('gene') or row.get('GENE')
    v = row.get('hgvsp') or row.get('HGVSp') or row.get('variant')

    hits = civic_hits(g, v)

    want = _aa_key(v)
    def hit_key(h):
        return _aa_key(h.get('variant') or h.get('variant_name') or h.get('name') or h.get('hgvsp'))
    exact = [h for h in hits if want and hit_key(h) == want]
    if exact:
        hits = exact

    if hits:
        tx_seen, tx = set(), []
        for h in hits:
            if not isinstance(h, dict):
                continue
            t = {
                "therapy":     (h.get("therapy") or h.get("drug") or h.get("agent") or h.get("treatment") or "").strip(),
                "evidence":    _norm_ev(h.get("evidence") or h.get("evidence_level") or h.get("evidence_label") or h.get("level") or ""),
                "cancer_type": (h.get("disease") or h.get("cancer") or h.get("tumor_type") or h.get("cancer_type") or "").strip(),
                "source_url":  (h.get("source_url") or h.get("url") or h.get("link") or "").strip(),
                "source":      (h.get("source") or h.get("database") or h.get("db") or "CIViC").strip(),
                "pmid":        (h.get("pmid") or h.get("PMID") or "").strip(),
                "direction":   (h.get("direction") or h.get("response") or h.get("response_type") or "").strip(),
                "match":       ("Exact" if want and hit_key(h) == want else "Gene-level"),
            }
            sig = (_eid(t["source_url"]), t["evidence"], t["therapy"].lower(), t["pmid"])
            if (t["therapy"] or t["evidence"] or t["source_url"] or t["pmid"]) and sig not in tx_seen:
                tx_seen.add(sig)
                tx.append(t)

        if tx:
            cur = row.get("treatments")
            if not isinstance(cur, list):
                cur = []
            cur.extend(tx)
            row["treatments"] = cur

        if not row.get("association"):
            assoc = next((t.get("direction") for t in tx if t.get("direction")), None)
            assoc = assoc or next((t.get("evidence") for t in tx if t.get("evidence")), None)
            if assoc:
                row["association"] = assoc

        if not row.get("cancer_type"):
            merged_tx = row.get("treatments") or []
            ct = next((t.get("cancer_type") for t in merged_tx if t.get("cancer_type")), None)
            if ct:
                row["cancer_type"] = ct

    before = (row.get("tier") or "")
    promote_tier(row)
    if before.upper().startswith("T1") and not row.get("treatments"):
        row["tier"] = before
    if row.get('tier') != before:
        promoted += 1

for row in R.get('rows', []):
    add_ct_to_comment(row)

out_json.write_text(json.dumps(R, ensure_ascii=False, indent=2), encoding='utf-8')
print(f"[civic] promoted_to_T1={promoted}", file=sys.stderr)
PY
  /$
}




process MAKE_JSON_SUMMARY {
  tag "summary_json"
  container "${IMG_PY}"
  cpus 6
  memory '24 GB'
  time '6h'
  publishDir "${params.outdir}/json", mode: 'copy', overwrite: true

  // env values as plain strings to keep Groovy happy
  env 'REPORT_USE_CSQ', (params.report_use_csq ? '1' : '0')
  env 'REPORT_CODING_ONLY',   (params.report_coding_only ? '1' : '0')
  env 'REPORT_KEEP_IMPACT',   "${params.report_keep_impact ?: 'HIGH,MODERATE'}"
  env 'REPORT_MAX_PER_GENE',  "${params.report_max_per_gene ?: 50}"
  env 'REPORT_MAX_ROWS',      "${params.report_max_rows     ?: 100000}"
  env 'REPORT_SEARCH_LIMIT',  "${params.report_search_limit ?: 500}"
  env 'ENABLE_EVIDENCE',      (params.enable_evidence ? '1' : '0')
  env 'ONCOKB_TOKEN',         "${params.oncokb_token ?: ''}"

  input:
    tuple val(sid), path(vep_vcf), path(pharmcat_json), path(panel_qc_json)
    path kb_tsv
    path hotspots_keys
    path pathway_db
  output:
    tuple val(sid), path("${sid}.report.json"),       emit: report_json
    tuple val(sid), path("${sid}.mcode_bundle.json"), emit: mcode_json

  shell:
  $/
  set -euo pipefail

  # Normalize inputs to predictable names
  cp -f "!{vep_vcf}"       vcf_in
  cp -f "!{pharmcat_json}" pharmcat.json  2>/dev/null || :
  cp -f "!{kb_tsv}"        knowledge.tsv  2>/dev/null || :
  cp -f "!{panel_qc_json}" panel_qc.json  2>/dev/null || :
  cp -f "!{hotspots_keys}" hotspots.keys.tsv 2>/dev/null || :   # optional
  cp -f "!{pathway_db}" pathway_db.json 2>/dev/null || :
  export SID='!{sid}'
  export PATIENT_ID='!{params.patient_id}'
  export VEP_ASSEMBLY='!{params.vep_assembly}'
  export ASSAY='!{params.assay}'
  export PANEL_NAME='!{ params.panel_name ?: "" }'
  
  python3 - <<'PY'
import os, json, gzip, re, statistics as stats
from pathlib import Path

# ---------------- helpers ----------------
def open_vcf(p: Path):
    with p.open('rb') as fh:
        magic = fh.read(2)

    # IMPORTANT: single backslashes here
    if magic == b'\x1f\x8b':
        return gzip.open(p, 'rt', encoding='utf-8', errors='replace')
    return p.open('rt', encoding='utf-8', errors='replace')

def parse_info(s):
    d={}
    for it in (s or '').split(';'):
        if not it: continue
        if '=' in it:
            k,v=it.split('=',1); d[k]=v
        else:
            d[it]=True
    return d
import re, json
from pathlib import Path

def load_pathway_map(path: Path):
    mp = {}
    if not path.exists() or path.stat().st_size == 0:
        return mp
    try:
        raw = json.loads(path.read_text(errors='ignore') or 'null')
    except Exception:
        return mp

    def add(g, pathways):
        g = (g or '').strip().upper()
        if not g:
            return
        parts = []
        if isinstance(pathways, str):
            parts = [re.sub(r'\s+', ' ', x.strip()) for x in re.split(r'[;,]', pathways) if x.strip()]
        elif isinstance(pathways, list):
            for x in pathways:
                if isinstance(x, str) and x.strip():
                    parts.append(re.sub(r'\s+', ' ', x.strip()))
        if parts:
            mp[g] = parts

    if isinstance(raw, list):
        items = raw
    elif isinstance(raw, dict):
        items = raw.get('genes') or raw.get('records') or raw.get('data')
        if not isinstance(items, list):
            # treat dict as {gene: [..]} or {gene: {"pathways": ..}}
            for k, v in raw.items():
                if isinstance(v, dict):
                    add(k, v.get('pathways') or v.get('pathway') or v.get('PATHWAYS'))
                else:
                    add(k, v)
            return mp
    else:
        return mp

    for rec in (items or []):
        if not isinstance(rec, dict):
            continue
        g  = rec.get('gene_symbol') or rec.get('gene') or rec.get('GENE') or rec.get('symbol')
        ps = rec.get('pathways') or rec.get('pathway') or rec.get('PATHWAYS') or rec.get('Pathway') or rec.get('pathways_list')
        add(g, ps)
    return mp


PATHWAY_MAP = load_pathway_map(Path("pathway_db.json"))

def parse_fmt(fmt, sample):
    m={}
    if not fmt or not sample: return m
    ks=fmt.split(':'); vs=sample.split(':')
    for i,k in enumerate(ks):
        if i<len(vs): m[k]=vs[i]
    return m

def first_ab(ad):
    if not ad: return None
    try:
        xs=[int(x) for x in ad.split(',') if x!='.']
        if len(xs)>=2:
            ref,alt=xs[0],xs[1]; tot=ref+alt
            return (alt/tot) if tot>0 else None
    except: pass
    return None

def terms_from(c):
    for t in (c or '').replace(',', '&').split('&'):
        t=t.strip()
        if t: yield t

AA3_TO1 = {
    'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Gln':'Q','Glu':'E','Gly':'G',
    'His':'H','Ile':'I','Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P','Ser':'S',
    'Thr':'T','Trp':'W','Tyr':'Y','Val':'V','Ter':'*','Stop':'*'
}

def normalize_hgvsp_short(hgvsp: str):
    if not hgvsp:
        return None
    s = hgvsp.split(':', 1)[-1].strip()
    if not s.startswith('p.'):
        return None
    for k, v in AA3_TO1.items():
        s = re.sub(rf'\b{k}\b', v, s)
    return s


def load_hotspot_keys(path: Path):
    """Read a TSV of (Gene, HGVSp_short) and return a set of (GENE, p.X123Y)."""
    pairs = set()
    if not path.exists() or path.stat().st_size == 0:
        return pairs
    lines = path.read_text(errors='ignore').splitlines()
    if not lines:
        return pairs

    hdr = [c.strip().lower() for c in lines[0].split('\t')]
    header_like = any(c in {'hugo_symbol','gene','symbol','hgvsp_short','protein_change','hgvsp','aa_change','protein'} for c in hdr)

    if header_like:
        start = 1
        col_gene  = next((i for i,c in enumerate(hdr) if c in {'hugo_symbol','gene','symbol'}), None)
        col_hgvsp = next((i for i,c in enumerate(hdr) if c in {'hgvsp_short','protein_change','hgvsp','aa_change','protein'}), None)
        if col_gene is None or col_hgvsp is None:
            return pairs
    else:
        start = 0
        col_gene, col_hgvsp = 0, 1

    for ln in lines[start:]:
        if not ln.strip():
            continue
        parts = ln.split('\t')
        if len(parts) <= max(col_gene, col_hgvsp):
            continue
        g = parts[col_gene].strip().upper()
        aa_norm = normalize_hgvsp_short(parts[col_hgvsp].strip())
        if g and aa_norm:
            pairs.add((g, aa_norm))
    return pairs

hotspot_pairs = load_hotspot_keys(Path("hotspots.keys.tsv"))



def safe_get(d, *keys, default=None):
    for k in keys:
        if isinstance(d, dict) and k in d: return d[k]
    return default

def qc_bucket(value, good_min=None, warn_min=None, good_max=None, warn_max=None, fmt=None, label=None):
    if value is None:
        return {"value": None, "status": "na", "comment": "No value", "thresholds": {}}
    val_disp = f"{value:.1f}" if fmt == "1dp" else f"{value:.0f}"
    th = {}
    if good_min is not None: th["good_min"] = good_min
    if warn_min is not None: th["warn_min"] = warn_min
    if good_max is not None: th["good_max"] = good_max
    if warn_max is not None: th["warn_max"] = warn_max

    status = "good"
    if warn_min is not None and value < warn_min: status = "poor"
    elif good_min is not None and value < good_min: status = "warn"
    if warn_max is not None and value > warn_max: status = "poor"
    elif good_max is not None and value > good_max: status = "warn"

    direction = "below" if (good_min is not None and value < good_min) else ("above" if (good_max is not None and value > good_max) else "near")
    comment = f"{label or 'Metric'} {val_disp} meets target" if status == "good" else f"{label or 'Metric'} {val_disp} {direction} target"
    return {"value": value, "status": status, "thresholds": th, "comment": comment}

# --------------- params/env ----------------
sid        = os.environ.get('SID','')
patient_id = os.environ.get('PATIENT_ID','')
assembly   = os.environ.get('VEP_ASSEMBLY','')
assay      = os.environ.get('ASSAY','')
panel_name = os.environ.get('PANEL_NAME') or None

MIN_QUAL     = float(os.getenv('REPORT_MIN_QUAL', '30'))
MIN_DP       = int(os.getenv('REPORT_MIN_DP',   '10'))
REQUIRE_PASS = os.getenv('REPORT_REQUIRE_PASS', '1') == '1'
APPLY_AB     = os.getenv('REPORT_APPLY_AB',     '0') == '1'
AB_MIN       = float(os.getenv('REPORT_AB_MIN', '0.25'))
AB_MAX       = float(os.getenv('REPORT_AB_MAX', '0.75'))
CODING_ONLY  = os.getenv('REPORT_CODING_ONLY',  '0') == '1'
KEEP_IMPACT  = {s.strip().upper() for s in os.getenv('REPORT_KEEP_IMPACT','').split(',') if s.strip()}
MAX_PER_GENE = int(os.getenv('REPORT_MAX_PER_GENE','0') or '0')
MAX_ROWS     = int(os.getenv('REPORT_MAX_ROWS','0') or '0')
ENABLE_EVID  = os.getenv('ENABLE_EVIDENCE','0') in {'1','true','TRUE','yes','on'}

vcf_path = Path("vcf_in")
phc_path = Path("pharmcat.json")
panel_qc = Path("panel_qc.json")
hotkeys  = Path("hotspots.keys.tsv")
kb_path  = Path("knowledge.tsv")



# ------------- load TSG list (optional) -------------
tsg_genes=set()
if kb_path.exists() and kb_path.stat().st_size>0:
    try:
        header=None
        for ln in kb_path.read_text(errors='ignore').splitlines():
            if header is None:
                header=[c.strip().lower() for c in ln.split('\\t')]
                col_gene = next((i for i,c in enumerate(header) if c in {'gene','symbol','hugo_symbol'}), None)
                # role-like columns
                col_role = next((i for i,c in enumerate(header) if 'tsg' in c or 'role' in c or 'class' in c), None)
                continue
            if not ln.strip(): continue
            parts=ln.split('\\t')
            try:
                g = parts[col_gene].strip().upper() if col_gene is not None else ''
                role = (parts[col_role].strip().lower() if col_role is not None else '')
            except IndexError:
                continue
            if g and ('tsg' in role or 'tumor' in role or 'suppress' in role):
                tsg_genes.add(g)
    except Exception:
        pass

# --------------- VEP typing helpers ---------------
vep_priority = [
    "transcript_ablation","splice_acceptor_variant","splice_donor_variant",
    "stop_gained","frameshift_variant","stop_lost","start_lost",
    "inframe_insertion","inframe_deletion","missense_variant",
    "protein_altering_variant","synonymous_variant"
]
type_map = {
    "missense_variant":"missense",
    "frameshift_variant":"frameshift",
    "stop_gained":"nonsense",
    "stop_lost":"stop_lost",
    "start_lost":"start_lost",
    "inframe_insertion":"inframe_ins",
    "inframe_deletion":"inframe_del",
    "protein_altering_variant":"protein_altering",
    "synonymous_variant":"synonymous",
    "splice_acceptor_variant":"splice_acceptor",
    "splice_donor_variant":"splice_donor"
}
codingish = set(type_map.keys()) | {"protein_altering_variant"}
order = {t:i for i,t in enumerate(vep_priority)}

def best_vep_term(anns):
    best=None; best_i=10**9
    for a in (anns or []):
        for t in terms_from(a.get('Consequence')):
            i=order.get(t, 10**6)
            if i<best_i:
                best_i=i; best=t
    return best

def assign_tier(r):
    # T1 if guideline/label evidence exists (only if enrichment ever adds treatments)
    ev_texts = [ (t.get("evidence") or "") + " " + (t.get("source") or "")
                 for t in (r.get("treatments") or []) ]
    if any(re.search(r'(?:FDA|EMA|NCCN|ESMO|Guideline|Level\\s*[AB]\\b)', e, re.I) for e in ev_texts):
        return "T1"
    # T2 for hotspot/known oncogenic/TSG truncation
    if r.get("known_label"):
        return "T2"
    return "T3+"

# --------------- parse VCF ----------------
rows=[]; csq_fields=[]; gene_counts={}
csq_re = re.compile(r'##INFO=<ID=CSQ[^>]*Description="[^"]*Format:\s*([^"]+)"')

with open_vcf(vcf_path) as fh:
    for raw in fh:
        line = raw.lstrip(" \t\r\ufeff")

        if not line.strip():
            continue

        # --- parse VEP CSQ header to learn field order ---
        if line.startswith("##INFO=<ID=CSQ"):
            m = csq_re.search(line)
            if m:
                csq_fields = [x.strip() for x in m.group(1).split("|")]
            continue

        # pass other meta/header lines
        if line.startswith("#"):
            continue

        parts = line.rstrip('\n').split('\t')
        if len(parts) < 8:
            continue


        chrom, pos, _id, ref, alt, qstr, filt, info = parts[:8]
        fmt  = parts[8] if len(parts) > 8 else ''
        samp = parts[9] if len(parts) > 9 else ''
        info_d = parse_info(info); fmt_d = parse_fmt(fmt, samp)

        try: q = None if qstr=='.' else float(qstr)
        except: q = None

        # depth
        dp = None
        if fmt_d.get('DP') and fmt_d['DP']!='.':
            try: dp=int(fmt_d['DP'])
            except: dp=None
        if dp is None and info_d.get('DP') and info_d['DP']!='.':
            try: dp=int(info_d['DP'])
            except: dp=None
        if dp is None and fmt_d.get('AD'):
            try:
                ds=[int(x) for x in fmt_d['AD'].split(',') if x!='.']
                dp = sum(ds) if ds else None
            except: dp=None

        gt = fmt_d.get('GT') or ''
        ab = first_ab(fmt_d.get('AD'))

        # CSQ ann
        anns=[]
        csq = None
        for it in info.split(';'):
            if it.startswith('CSQ='):
                csq = it[4:]; break
        if csq and csq_fields:
            for rec in csq.split(','):
                vs = rec.split('|')
                anns.append({k: (vs[i] if i<len(vs) else '') for i,k in enumerate(csq_fields)})

        gene = ''
        cons = ''
        hgvsp = ''
        impact = ''
        if anns:
            a0 = anns[0]
            gene  = a0.get('SYMBOL','') or ''
            cons  = a0.get('Consequence','') or ''
            hgvsp = a0.get('HGVSp','') or a0.get('HGVSp_VEP','') or ''
            impact = (a0.get("IMPACT") or "")

        row = {
            "chrom": chrom, "pos": int(pos), "ref": ref, "alt": alt,
            "gene": gene, "consequence": cons, "hgvsp": hgvsp,
            "qual": q, "filter": filt, "depth": dp,
            "type": None, "coding": None,
            "impact": impact
        }
        h_short = normalize_hgvsp_short(hgvsp)  # e.g. p.R175H (if possible)
        row["hgvsp_short"] = h_short

        # Observed Variant display name for the UI
        if gene and h_short:
            row["variant"] = f"{gene} {h_short}"
        elif gene and hgvsp:
            # keep full if short fails
            row["variant"] = f"{gene} {hgvsp.split(':')[-1]}"
        else:
            # fallback to locus-style
            row["variant"] = f"{chrom}:{pos} {ref}>{alt}"

        # some UIs use different keys
        row["observed_variant"] = row["variant"]

        if ab is not None:
            row["ab"] = round(ab, 4)

        # VEP typing
        best = best_vep_term(anns)
        row["type"]   = type_map.get(best, best or "other")
        row["coding"] = any(t in codingish for a in anns for t in terms_from(a.get('Consequence')))

        # --- Hotspot detection + pathways mapping ---
        g_up    = (gene or '').strip().upper()
        h_short = normalize_hgvsp_short(hgvsp)  # e.g. ENSP:… → p.V600E
        keys_str = (info_d.get('MSK_HOTSPOT_KEY') or '').strip()

        # Pathways (per-gene)
        if g_up in PATHWAY_MAP:
            row["pathways"] = PATHWAY_MAP[g_up]

        # Hotspot (MSK)
        if keys_str:
            row["known_label"] = "Hotspot"
            row.setdefault("tags", []).append("MSK_hotspot")
            row["hotspot_keys"] = [s for s in keys_str.split(",") if s]
        elif g_up and h_short and (g_up, h_short) in hotspot_pairs:
            row["known_label"] = "Hotspot"
            row.setdefault("tags", []).append("MSK_hotspot")

        # TSG truncation rule (if knowledge base provides TSG info)
        trunc_terms = {"stop_gained","frameshift_variant","splice_acceptor_variant","splice_donor_variant","start_lost"}
        if g_up and g_up in tsg_genes and any(t in trunc_terms for a in anns for t in terms_from(a.get('Consequence'))):
            row["known_label"] = row.get("known_label") or "TSG_trunc"

        # Metrics bundle for UI
        row["metrics"] = {"DP": dp, "QUAL": q, "VAF": ab}
        ad_raw = fmt_d.get('AD')
        if ad_raw and ad_raw != '.':
            try:
                row["metrics"]["AD"] = [int(x) for x in ad_raw.split(',') if x != '.'][:2]
            except:
                pass

        # Assign tier AFTER known_label is set
        row["tier"] = assign_tier(row)
        

        # gating filters (skip append if fails)
        if REQUIRE_PASS and (filt or '').upper() not in {'PASS','.'}:
            continue
        if (q is not None) and (q < MIN_QUAL):
            continue
        if dp is None or dp < MIN_DP:
            continue
        if APPLY_AB and gt in {'0/1','1/0'} and ab is not None and not (AB_MIN <= ab <= AB_MAX):
            continue
        if CODING_ONLY and csq_fields:
            coding_terms = {"missense_variant","frameshift_variant","stop_gained","stop_lost","start_lost",
                            "inframe_insertion","inframe_deletion","protein_altering_variant","synonymous_variant"}
            all_cons = ' '.join(a.get('Consequence','') for a in anns)
            if not any(t in all_cons for t in coding_terms):
                continue
        if KEEP_IMPACT:
            ims = {(a.get('IMPACT') or '').upper() for a in anns if a.get('IMPACT')}
            if not ims.intersection(KEEP_IMPACT):
                continue

        if MAX_PER_GENE and g_up:
            n = gene_counts.get(g_up, 0)
            if n >= MAX_PER_GENE:
                continue
            gene_counts[g_up] = n + 1
        row["is_hotspot"]   = (row.get("known_label") == "Hotspot")
        tier = str(row.get("tier", "")).strip().upper()   # e.g. "T1", "T3+"
        m = re.match(r"^T(\d+)", tier)                    # captures 1,2,3...
        row["tier_numeric"] = int(m.group(1)) if m else 9
        rows.append(row)
        if MAX_ROWS and len(rows) >= MAX_ROWS:
            break


# -------- Summary --------
snv = sum(1 for r in rows if len(r["ref"])==1 and len(r["alt"])==1)
ind = len(rows) - snv
def is_ts(r,a): return (r,a) in {("A","G"),("G","A"),("C","T"),("T","C")}
titv_num = sum(1 for r in rows if len(r["ref"])==1 and len(r["alt"])==1 and is_ts(r["ref"],r["alt"]))
titv_den = max(1, snv - titv_num)
titv = (titv_num / titv_den) if snv else None
try:
    med_dp = stats.median([r["depth"] for r in rows if r.get("depth") is not None]) if rows else None
except:
    med_dp=None

# QC evaluation (with comments)
qc_panel = None
try:
    if panel_qc.exists() and panel_qc.stat().st_size > 0:
        qc_panel = json.loads(panel_qc.read_text())
except Exception:
    qc_panel = {"_error": "Failed to read panel_qc.json"}

qc_eval = {"status": "na", "metrics": []}
if isinstance(qc_panel, dict):
    qc_metrics = [
        {"name": "Mean coverage",   **qc_bucket(safe_get(qc_panel,"mean_coverage","mean_cov", default=None),   good_min=50, warn_min=30, fmt="1dp", label="Mean coverage (×)"), "unit":"x"},
        {"name": "Median coverage", **qc_bucket(safe_get(qc_panel,"median_coverage","median_cov", default=None),good_min=50, warn_min=30, fmt="0dp", label="Median coverage (×)"), "unit":"x"},
        {"name": "≥20× (%)",        **qc_bucket(safe_get(qc_panel,"ge20_pct","pct_20x","ge20x_pct", default=None),        good_min=80, warn_min=60, fmt="0dp", label="≥20× (%)"), "unit":"%"},
        {"name": "≥50× (%)",        **qc_bucket(safe_get(qc_panel,"ge50_pct","pct_50x","ge50x_pct", default=None),        good_min=60, warn_min=40, fmt="0dp", label="≥50× (%)"), "unit":"%"},
        {"name": "≥100× (%)",       **qc_bucket(safe_get(qc_panel,"ge100_pct","pct_100x","ge100x_pct", default=None),     good_min=30, warn_min=15, fmt="0dp", label="≥100× (%)"), "unit":"%"},
        {"name": "On-target (%)",   **qc_bucket(safe_get(qc_panel,"on_target_rate_pct","on_target_pct","ontarget_pct", default=None), good_min=50, warn_min=40, fmt="1dp", label="On-target (%)"), "unit":"%"},
        {"name": "Duplication (%)", **qc_bucket(safe_get(qc_panel,"duplication_pct","dup_pct", default=None),             good_max=20, warn_max=30, fmt="1dp", label="Duplication (%)"), "unit":"%"},
        {"name": "Q30 (%)",         **qc_bucket(safe_get(qc_panel,"q30_pct","pct_q30_bases", default=None),               good_min=85, warn_min=80, fmt="1dp", label="Q30 (%)"), "unit":"%"},
    ]

    worst = "good"
    for m in qc_metrics:
        if m["status"] == "poor":
            worst = "poor"
            break
        if m["status"] == "warn":
            worst = "warn"
    qc_eval = {"status": worst, "metrics": qc_metrics}

# Top variants summary (independent of qc_panel)
summary_top = []
top = [r for r in rows if (r.get("tier_numeric", 9) <= 2)]
top.sort(key=lambda r: (r.get("tier_numeric", 9),
                        -(r.get("metrics", {}).get("DP") or 0),
                        -(r.get("metrics", {}).get("VAF") or 0)))
summary_top = top[:50]


# -------- Report JSON --------
report = {
  "meta": {
    "report_version": "1.0.0",
    "reference": assembly,
    "pipeline": "OpenCare  (panel-enabled)",
    "pipeline_version": "v1.0.0",
    "panel_name": panel_name,
    "assay": assay,
    "has_vep": True if os.path.exists(str(vcf_path)) else False,
    "filters": {"min_qual": MIN_QUAL, "min_dp": MIN_DP, "require_pass": REQUIRE_PASS,
                "apply_ab": APPLY_AB, "ab_min": AB_MIN, "ab_max": AB_MAX, "coding_only": CODING_ONLY},
    "disclaimer": "Not for clinical decision-making without local validation.",
    "evidence_enabled": ENABLE_EVID
  },
  "patient": {"id": patient_id, "report_date": None},
  "alt_id": sid,
  "qc": {"panel": qc_panel, "panel_eval": qc_eval},
  "summary": {
      "tumor_content": None,
      "microbial_species": None,
      "mutation_signatures": [],
      "sv_burden": 0,
      "counts": {"snv": snv, "indel": ind, "titv": titv, "median_dp": med_dp},
      "top_variants": summary_top
  },
  "rows": rows,
  "pharmcat": (json.loads(Path("pharmcat.json").read_text()) if (Path("pharmcat.json").exists() and Path("pharmcat.json").stat().st_size>0) else None)
}

Path(f"{sid}.report.json").write_text(json.dumps(report, indent=2))

# -------- mCODE bundle --------
bundle={"resourceType":"Bundle","type":"collection","entry":[{"resource":{"resourceType":"Patient","id":patient_id or sid}}]}
for i,r in enumerate(rows,1):
    bundle["entry"].append({
      "resource":{
        "resourceType":"Observation","id":f"var-{i}","status":"final",
        "code":{"text":"Genetic variant"},
        "subject":{"reference":f"Patient/{patient_id or sid}"},
        "component":[
          {"code":{"text":"Gene"},"valueString":r.get("gene","")},
          {"code":{"text":"HGVSp"},"valueString":r.get("hgvsp","")},
          {"code":{"text":"Consequence"},"valueString":r.get("consequence","")},
          {"code":{"text":"Location"},"valueString":f"{r['chrom']}:{r['pos']} {r['ref']}>{r['alt']}"}
        ]
      }
    })
Path(f"{sid}.mcode_bundle.json").write_text(json.dumps(bundle, indent=2))
PY

  ls -lah
  /$
}


process CALL_ARM_EVENTS {
  tag { sid }
  container 'python:3.11'
  cpus 1
  memory '1 GB'
  time '15m'
  stageInMode 'copy'
  publishDir "${params.outdir}/cn", mode: 'copy', overwrite: true

  input:
    tuple val(sid), path(segments_tsv), path(arm_py)   // <— script staged here
    path  cytoband

  output:
    tuple val(sid), path("${sid}.arm_calls.json"), path("${sid}.arm_calls.txt"), emit: arms

  shell:
  """
  set -euo pipefail

  # Sanitize cytoband: drop blank lines and rows with empty band 'name' (4th column).
  awk -F '\\t' 'NF>=4 && \$4!="" {print}' "!{cytoband}" | tr -d '\\r' > cytoband.clean.tsv

  python3 "!{arm_py}" \\
    --segments "!{segments_tsv}" \\
    --cytoband "cytoband.clean.tsv" \\
    --log2-gain "!{params.arm_gain_log2}" \\
    --log2-loss "!{params.arm_loss_log2}" \\
    --min-cover-frac "!{params.arm_cover_min}" \\
    --out-json "!{sid}.arm_calls.json" \\
    --out-tags "!{sid}.arm_calls.txt"
  """
}



process AUGMENT_JSON_ARM_EVENTS {
  tag { sid }
  container 'python:3.11'
  cpus 1
  memory '1 GB'
  time '5m'
  stageInMode 'copy'
  publishDir "${params.outdir}/cn", mode: 'copy', overwrite: true
  input:
    tuple val(sid), path(slim_json), path(full_json), path(arm_json)

  output:
    tuple val(sid),
          path("${sid}.report.slim.arm.json"),
          path("${sid}.report.full.arm.json"), emit: json

  shell:
  """
  set -euo pipefail
  python3 - <<'PY' "!{slim_json}" "!{full_json}" "!{arm_json}" "!{sid}"
import json, sys

slim_p, full_p, arms_p, sid = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

with open(slim_p, encoding="utf-8") as f: slim = json.load(f)
with open(full_p, encoding="utf-8") as f: full = json.load(f)
try:
    with open(arms_p, encoding="utf-8") as f: arms = json.load(f) or {}
except Exception:
    arms = {}

for rep in (slim, full):
    rep.setdefault("summary", {})
    rep["summary"]["arm_cn"] = arms.get("arms", {})
    rep["summary"]["arm_cn_summary"] = arms.get("summary", "")

open(f"{sid}.report.slim.arm.json","w",encoding="utf-8").write(json.dumps(slim, ensure_ascii=False, indent=2))
open(f"{sid}.report.full.arm.json","w",encoding="utf-8").write(json.dumps(full, ensure_ascii=False, indent=2))
PY
  """
}

process MAKE_TN_SEGMENTS_FROM_MOSDEPTH {
  tag { "${tsid}_vs_${nsid}" }
  container 'python:3.11'
  stageInMode 'copy'
  publishDir "${params.outdir}/cn", mode: 'copy', overwrite: true

  input:
    tuple val(tsid),
          path(t_regions, stageAs: "tumor.regions.bed.gz"),
          val(nsid),
          path(n_regions, stageAs: "normal.regions.bed.gz")

  output:
    tuple val(tsid), val(nsid), path("${tsid}.segments.tsv")

  shell:
  '''
  set -euo pipefail

  python3 - "tumor.regions.bed.gz" "normal.regions.bed.gz" "!{tsid}.segments.tsv" <<'PY'
import gzip, math, sys, re

t_path, n_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]

def read_mosdepth_regions(p):
    d = {}
    with gzip.open(p, "rt") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()   # <-- handles tabs OR spaces
            if len(parts) < 4:
                continue
            chrom = parts[0]
            try:
                start = int(parts[1]); end = int(parts[2])
            except ValueError:
                continue

            mean = None
            for x in reversed(parts):
                try:
                    mean = float(x)
                    break
                except ValueError:
                    continue
            if mean is None:
                continue

            d[(chrom, start, end)] = mean
    return d

def chrom_key(ch):
    c = ch[3:] if ch.startswith("chr") else ch
    if c.isdigit():
        return (0, int(c))
    cU = c.upper()
    if cU == "X": return (1, 23)
    if cU == "Y": return (1, 24)
    if cU in ("M","MT"): return (1, 25)
    return (2, cU)

T = read_mosdepth_regions(t_path)
N = read_mosdepth_regions(n_path)

keys = sorted(set(T.keys()) & set(N.keys()), key=lambda k: (chrom_key(k[0]), k[1], k[2]))

PC = 1.0
with open(out_path, "w", encoding="utf-8") as out:
    for chrom, start, end in keys:
        t = T[(chrom, start, end)]
        n = N[(chrom, start, end)]
        log2 = math.log2((t + PC) / (n + PC))
        print(chrom, start, end, f"{log2:.6f}", sep="\t", file=out)

print(f"[tn_segments] wrote={out_path} n={len(keys)}", file=sys.stderr)
PY
  '''
}

process MAKE_BENCH_METRICS {
  tag "bench_metrics"
  container 'python:3.11'
  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
    path trace_tsv

  output:
    path "bench_metrics.tsv"

  script:
  """
  python3 -m pip install --no-cache-dir pandas
  python3 - <<'PY'
  import pandas as pd, re

  df = pd.read_csv("${trace_tsv}", sep="\\t")

  def to_seconds(x):
    if pd.isna(x): return 0.0
    s = str(x)
    try:
      return pd.to_timedelta(s).total_seconds()
    except Exception:
      pass
    # fallback for strings like "1h 2m 3s"
    mult = {'s':1,'m':60,'h':3600,'d':86400}
    parts = re.findall(r'(\\d+(?:\\.\\d+)?)([smhd])', s)
    return sum(float(v)*mult[u] for v,u in parts) if parts else 0.0

  if 'realtime' in df.columns:
    df['realtime_s'] = df['realtime'].apply(to_seconds)
  else:
    df['realtime_s'] = 0.0

  cpus = df['cpus'].fillna(1) if 'cpus' in df.columns else 1
  cpu_hours  = (cpus * df['realtime_s']).sum() / 3600.0
  wall_hours = df['realtime_s'].sum() / 3600.0

  completed = int((df['status'] == 'COMPLETED').sum()) if 'status' in df.columns else len(df)

  rows = [
    ("tasks_total", len(df)),
    ("tasks_completed", completed),
    ("wall_hours_total", round(wall_hours, 4)),
    ("cpu_hours_total", round(cpu_hours, 4)),
  ]

  with open("bench_metrics.tsv","w") as f:
    f.write("metric\\tvalue\\n")
    for k,v in rows:
      f.write(f"{k}\\t{v}\\n")
  PY
  """
}




process BENCHMARK_REPORT {
  tag "benchmark_report"
  container 'python:3.11'
  cpus 1
  memory '1 GB'
  time '10m'
  stageInMode 'copy'
  publishDir "${params.outdir ?: 'results'}/benchmark", mode: 'copy', overwrite: true

  input:
    val  gate
    path trace_tsv
    path bench_metrics_tsv
    path acc_metrics_tsv
    val  outdir_path

  output:
    path "workstation_specs.txt"
    path "table2_benchmark_hcc1395.tsv"
    path "benchmark_report.tsv"
    path "benchmark_trace_by_process.tsv"

  when:
    params.enable_benchmark

  shell:
  $/
  set -euo pipefail

  export TRACE="!{trace_tsv}"
  export RUNTIME="!{bench_metrics_tsv}"
  export ACC="!{acc_metrics_tsv}"
  export OUTDIR="!{outdir_path}"

  python3 - <<'PY'
import csv, os, re, subprocess
from datetime import datetime
from pathlib import Path
from collections import defaultdict

TRACE   = Path(os.environ["TRACE"])
RUNTIME = Path(os.environ["RUNTIME"])
ACC     = Path(os.environ.get("ACC",""))
OUTDIR  = os.environ.get("OUTDIR",".")

def sh(cmd):
    try:
        return subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT, text=True).strip()
    except Exception:
        return ""

def truthy(s):
    return str(s or "").strip().lower() in ("true","1","yes","y")

def parse_dt(s):
    if not s or s == "-":
        return None
    s = str(s).strip()
    fmts = [
        "%Y-%m-%d %H:%M:%S.%f",
        "%Y-%m-%d %H:%M:%S",
        "%Y-%m-%dT%H:%M:%S.%f",
        "%Y-%m-%dT%H:%M:%S",
    ]
    for f in fmts:
        try:
            return datetime.strptime(s, f)
        except Exception:
            pass
    return None

def parse_duration_to_seconds(s):
    if not s or s == "-":
        return None
    s = str(s).strip()
    m = re.match(r'^(\d+):(\d+):(\d+(?:\.\d+)?)$', s)
    if m:
        h, mnt, sec = m.groups()
        return int(h)*3600 + int(mnt)*60 + float(sec)

    total = 0.0
    found = False
    for num, unit in re.findall(r"([0-9]*\.?[0-9]+)\s*(ms|s|m|h|d)", s):
        v = float(num); found = True
        if unit == "ms": total += v/1000.0
        elif unit == "s": total += v
        elif unit == "m": total += v*60.0
        elif unit == "h": total += v*3600.0
        elif unit == "d": total += v*86400.0
    if found:
        return total

    try:
        return float(s)
    except Exception:
        return None

def parse_mem_to_gb(s):
    if s is None:
        return None
    s = str(s).strip()
    if not s or s == "-":
        return None

    m = re.match(r'^([0-9.]+)\s*([KMGTP]?B)$', s, re.I)
    if m:
        v = float(m.group(1)); u = m.group(2).upper()
        if u == "KB": return v / (1024**2)
        if u == "MB": return v / 1024
        if u == "GB": return v
        if u == "TB": return v * 1024
        if u == "PB": return v * 1024 * 1024
        return None

    # fall back: bytes
    try:
        return float(s) / (1024**3)
    except Exception:
        return None

def parse_io_to_gb(s):
    if s is None: return None
    s = str(s).strip()
    if not s or s == "-": return None
    try:
        return float(s) / (1024**3)
    except Exception:
        return None

def proc_name(full):
    full = str(full or "")
    return full.split(" (", 1)[0].strip() if " (" in full else full.strip()

def df_stats(path):
    out = sh(f"df -Pk --output=size,avail '{path}' 2>/dev/null | tail -n 1")
    parts = out.split()
    if len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit():
        total_kb = float(parts[0]); avail_kb = float(parts[1])
        return (f"{total_kb/1024/1024:.2f}", f"{avail_kb/1024/1024:.2f}")

    out = sh(f"df -Pk '{path}' 2>/dev/null | tail -n 1")
    toks = out.split()
    nums = [t for t in toks if re.fullmatch(r"\d+", t)]
    if len(nums) >= 3:
        total_kb = float(nums[0])
        avail_kb = float(nums[2])
        return (f"{total_kb/1024/1024:.2f}", f"{avail_kb/1024/1024:.2f}")

    return ("","")

# ---------------- trace summary ----------------
tasks = []
start_dt = None
end_dt = None

tasks_total = tasks_ok = tasks_fail = tasks_cached = 0
sum_realtime_s = 0.0
sum_cpu_s = 0.0
max_realtime_s = 0.0
max_peak_rss_gb = 0.0
max_peak_vmem_gb = 0.0
sum_r_gb = 0.0
sum_w_gb = 0.0

by_proc = defaultdict(lambda: {"n":0,"cached":0,"rt_s":0.0,"cpu_s":0.0,"max_rss":0.0,"max_vmem":0.0})

with TRACE.open(newline='', encoding='utf-8', errors='ignore') as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        tasks.append(row)

for row in tasks:
    tasks_total += 1
    status = str(row.get("status","") or "").strip()
    status_u = status.upper()

    cached = (status_u == "CACHED") or truthy(row.get("cached",""))
    if cached: tasks_cached += 1

    if status_u in ("COMPLETED","CACHED"):
        tasks_ok += 1
    elif status_u in ("FAILED","ABORTED"):
        tasks_fail += 1

    sdt = parse_dt(row.get("start","")) or parse_dt(row.get("submit",""))
    edt = parse_dt(row.get("complete",""))
    if sdt and (start_dt is None or sdt < start_dt): start_dt = sdt
    if edt and (end_dt is None or edt > end_dt): end_dt = edt

    rt_s = parse_duration_to_seconds(row.get("realtime","")) or parse_duration_to_seconds(row.get("duration",""))
    if rt_s is not None:
        sum_realtime_s += rt_s
        max_realtime_s = max(max_realtime_s, rt_s)

    cpu_s = parse_duration_to_seconds(row.get("cputime","")) or parse_duration_to_seconds(row.get("cpu",""))
    if cpu_s is None and rt_s is not None:
        try:
            cpus = int(float(row.get("cpus","") or ""))
        except Exception:
            cpus = None
        if cpus is not None:
            cpu_s = rt_s * cpus
    if cpu_s is not None:
        sum_cpu_s += cpu_s

    prss_gb = parse_mem_to_gb(row.get("peak_rss","") or row.get("rss",""))
    pvm_gb  = parse_mem_to_gb(row.get("peak_vmem","") or row.get("vmem",""))
    if prss_gb is not None: max_peak_rss_gb = max(max_peak_rss_gb, prss_gb)
    if pvm_gb  is not None: max_peak_vmem_gb = max(max_peak_vmem_gb, pvm_gb)

    rchar = parse_io_to_gb(row.get("rchar",""))
    wchar = parse_io_to_gb(row.get("wchar",""))
    if rchar is not None: sum_r_gb += rchar
    if wchar is not None: sum_w_gb += wchar

    p = proc_name(row.get("name",""))
    by_proc[p]["n"] += 1
    if cached: by_proc[p]["cached"] += 1
    if rt_s is not None: by_proc[p]["rt_s"] += rt_s
    if cpu_s is not None: by_proc[p]["cpu_s"] += cpu_s
    if prss_gb is not None: by_proc[p]["max_rss"] = max(by_proc[p]["max_rss"], prss_gb)
    if pvm_gb  is not None: by_proc[p]["max_vmem"] = max(by_proc[p]["max_vmem"], pvm_gb)

wall_s = (end_dt - start_dt).total_seconds() if (start_dt and end_dt) else None

pipeline = {
    "nf_tasks": tasks_total,
    "nf_cached_tasks": tasks_cached,
    "nf_sum_realtime_h": round(sum_realtime_s/3600, 3),
    "nf_sum_cpu_h": round(sum_cpu_s/3600, 3),
    "nf_max_peak_rss_gb": round(max_peak_rss_gb, 3),
    "nf_max_peak_vmem_gb": round(max_peak_vmem_gb, 3),
    "nf_sum_rchar_gb": round(sum_r_gb, 3),
    "nf_sum_wchar_gb": round(sum_w_gb, 3),
}

# per-process
outp = Path("benchmark_trace_by_process.tsv")
hdrp = ["process","tasks","cached","realtime_h_sum","cpu_h_sum","peak_rss_gb_max","peak_vmem_gb_max"]
outp.write_text("\t".join(hdrp) + "\n", encoding="utf-8")
with outp.open("a", encoding="utf-8") as w:
    for p in sorted(by_proc.keys()):
        d = by_proc[p]
        w.write("\t".join([
            p, str(d["n"]), str(d["cached"]),
            f"{d['rt_s']/3600:.3f}",
            f"{d['cpu_s']/3600:.3f}",
            f"{d['max_rss']:.3f}",
            f"{d['max_vmem']:.3f}",
        ]) + "\n")

# workstation specs
wsl = bool(sh("grep -qi microsoft /proc/version && echo yes || true"))
cpu_model   = sh("lscpu | sed -n 's/^Model name:\\s*//p' | head -n 1")
cpu_sockets = sh("lscpu | sed -n 's/^Socket(s):\\s*//p' | head -n 1")
cpu_cores   = sh("lscpu | sed -n 's/^Core(s) per socket:\\s*//p' | head -n 1")
cpu_threads = sh("lscpu | sed -n 's/^CPU(s):\\s*//p' | head -n 1")

mem_total_gb = ""
mem_kb = sh("awk '/^MemTotal:/ {print $2}' /proc/meminfo 2>/dev/null")
if mem_kb.isdigit():
    mem_total_gb = f"{int(mem_kb)/1024/1024:.2f}"

work_total, work_avail = df_stats(".")
out_total, out_avail   = df_stats(OUTDIR)

gpu_name = sh("nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -n 1")
gpu_drv  = sh("nvidia-smi --query-gpu=driver_version --format=csv,noheader 2>/dev/null | head -n 1")
cuda_ver = sh("nvidia-smi 2>/dev/null | sed -n 's/.*CUDA Version: \\([0-9.]*\\).*/\\1/p' | head -n 1")
java_ver = sh("java -version 2>&1 | head -n 1")

os_pretty = ""
try:
    txt = Path("/etc/os-release").read_text(encoding="utf-8", errors="ignore")
    m = re.search(r'^PRETTY_NAME="?([^"\n]*)"?', txt, flags=re.M)
    if m: os_pretty = m.group(1).strip()
except Exception:
    pass

with open("workstation_specs.txt","w", encoding="utf-8") as o:
    def w(k,v): o.write(f"{k}: {v}\n")
    w("pipeline_start_time", start_dt.isoformat(sep=" ") if start_dt else "")
    w("pipeline_end_time", end_dt.isoformat(sep=" ") if end_dt else "")
    w("pipeline_wall_clock_seconds", f"{wall_s:.3f}" if wall_s is not None else "")
    w("tasks_total", tasks_total)
    w("tasks_succeeded", tasks_ok)
    w("tasks_failed", tasks_fail)
    w("tasks_cached", tasks_cached)
    w("total_cpu_time_seconds", f"{sum_cpu_s:.3f}")
    w("sum_realtime_seconds", f"{sum_realtime_s:.3f}")
    w("max_realtime_seconds", f"{max_realtime_s:.3f}")
    w("max_peak_rss_gb", f"{max_peak_rss_gb:.3f}")
    w("max_peak_vmem_gb", f"{max_peak_vmem_gb:.3f}")
    w("sum_rchar_gb", f"{sum_r_gb:.3f}")
    w("sum_wchar_gb", f"{sum_w_gb:.3f}")
    w("hostname", sh("hostname"))
    w("user", sh("whoami"))
    w("os_name", os_pretty)
    w("kernel_version", sh("uname -r"))
    w("wsl_detected", str(wsl).lower())
    w("cpu_model", cpu_model)
    w("cpu_sockets", cpu_sockets)
    w("cpu_cores_physical", cpu_cores)
    w("cpu_threads_logical", cpu_threads)
    w("mem_total_gb", mem_total_gb)
    w("disk_workdir_total_gb", work_total)
    w("disk_workdir_avail_gb", work_avail)
    w("disk_outdir_total_gb", out_total)
    w("disk_outdir_avail_gb", out_avail)
    w("gpu_name", gpu_name)
    w("gpu_driver_version", gpu_drv)
    w("cuda_version", cuda_ver)
    w("java_version", java_ver)

# benchmark_report.tsv merge
bench_rows = []
bench_hdr = []
if RUNTIME.exists() and RUNTIME.stat().st_size > 0:
    with RUNTIME.open(newline='', encoding='utf-8', errors='ignore') as f:
        r = csv.DictReader(f, delimiter='\t')
        bench_hdr = r.fieldnames or []
        for row in r:
            row.update(pipeline)
            bench_rows.append(row)

out = Path("benchmark_report.tsv")
final_hdr = (bench_hdr or []) + [k for k in pipeline.keys() if k not in (bench_hdr or [])]
out.write_text("\t".join(final_hdr) + "\n", encoding="utf-8")
with out.open("a", encoding="utf-8") as w:
    for row in bench_rows:
        w.write("\t".join(str(row.get(k,"")) for k in final_hdr) + "\n")

# table2
ROWS = [
    ("High-confidence ∩ callable (DP≥10 in tumor & normal)", "SNV",  "Overall"),
    ("High-confidence ∩ callable (DP≥10 in tumor & normal)", "Indel","Overall"),
    ("Protein-coding exons ±10 bp (restricted)",             "SNV",  "Overall"),
    ("Protein-coding exons ±10 bp (restricted)",             "Indel","Overall"),
    ("High-confidence ∩ callable (DP≥10 in tumor & normal)", "SNV",  "<5%"),
    ("High-confidence ∩ callable (DP≥10 in tumor & normal)", "SNV",  "5–10%"),
    ("High-confidence ∩ callable (DP≥10 in tumor & normal)", "SNV",  "10–20%"),
    ("High-confidence ∩ callable (DP≥10 in tumor & normal)", "SNV",  "≥20%"),
    ("High-confidence ∩ callable (DP≥10 in tumor & normal)", "Indel","<5%"),
    ("High-confidence ∩ callable (DP≥10 in tumor & normal)", "Indel","5–10%"),
    ("High-confidence ∩ callable (DP≥10 in tumor & normal)", "Indel","10–20%"),
    ("High-confidence ∩ callable (DP≥10 in tumor & normal)", "Indel","≥20%"),
]

metrics_map = {}
if ACC and str(ACC) and ACC.exists() and ACC.stat().st_size > 0:
    with ACC.open(newline='', encoding='utf-8', errors='ignore') as f:
        r = csv.DictReader(f, delimiter='\t')
        fields = r.fieldnames or []
        f_lut = {x.lower(): x for x in fields}
        need = {"region_restriction","variant_class","vaf_stratum","precision","recall","f1"}
        if need.issubset(set(f_lut.keys())):
            RR = f_lut["region_restriction"]
            VC = f_lut["variant_class"]
            VS = f_lut["vaf_stratum"]
            PR = f_lut["precision"]
            RE = f_lut["recall"]
            F1 = f_lut["f1"]
            for row in r:
                key = (row.get(RR,""), row.get(VC,""), row.get(VS,""))
                metrics_map[key] = (row.get(PR,""), row.get(RE,""), row.get(F1,""))

with open("table2_benchmark_hcc1395.tsv","w", newline="", encoding="utf-8") as o:
    w = csv.writer(o, delimiter='\t')
    w.writerow(["Region restriction","Variant class","VAF stratum","Precision","Recall","F1"])
    for k in ROWS:
        p, r, f1 = metrics_map.get(k, ("","",""))
        w.writerow([k[0], k[1], k[2], p, r, f1])
PY
  /$
}
