nextflow.enable.dsl = 2

/*
  OpenCare Clinical Core 
  FASTQ -> FastQC -> fastp -> BWA-MEM -> sort/index -> bcftools call
  -> (optional) VEP -> (optional) PharmCAT -> JSON/mCODE -> Interactive HTML

  Notes / run-time hygiene:
   - Every script block starts with `set -euo pipefail`
   - FASTQC handles missing/SE R2 safely
   - BWA index step re-uses prebuilt indices if present
*/

/* ---------- params (safe defaults, silence WARNs) ---------- */
if( !params.containsKey('outdir') )         params.outdir         = 'results'
if( !params.containsKey('reads') )          params.reads          = null
if( !params.containsKey('ref_fa') )         params.ref_fa         = null
if( !params.containsKey('region') )         params.region         = null
if( !params.containsKey('max_cpus') )       params.max_cpus       = 8
if( !params.containsKey('max_mem') )        params.max_mem        = '24 GB'
if( !params.containsKey('max_time') )       params.max_time       = '90h'
if( !params.containsKey('report_coding_only') ) params.report_coding_only = true
if( params.panel_bed && !file(params.panel_bed).exists() )
  error "panel_bed not found: ${params.panel_bed}"

if( !params.containsKey('vep_cache') )      params.vep_cache      = null
if( !params.containsKey('vep_assembly') )   params.vep_assembly   = 'GRCh38'
if( !params.containsKey('patient_id') )     params.patient_id     = 'PATIENT01'
if( !params.containsKey('pharmcat_jar') )   params.pharmcat_jar   = null
if( !params.containsKey('knowledge_tsv') )  params.knowledge_tsv  = null
if( !params.containsKey('pathway_db') )   params.pathway_db   = null
if( !params.containsKey('gene_domains') ) params.gene_domains = null

/* these aren’t used now, but defining them kills WARNs */
if( !params.containsKey('mcode') )          params.mcode          = null
if( !params.containsKey('oncokb_token') )   params.oncokb_token   = null
if( !params.containsKey('ncbi_key') )       params.ncbi_key       = null
if( !params.containsKey('ncbi_email') )     params.ncbi_email     = null
if( !params.containsKey('enable_online') )  params.enable_online  = null
if( !params.containsKey('targets_bed') )    params.targets_bed    = null

/* --- panel options (optional) --- */
if( !params.containsKey('panel_bed') )      params.panel_bed      = null
if( !params.containsKey('panel_name') )     params.panel_name     = null
if( !params.containsKey('assay') )          params.assay          = (params.panel_bed ? 'panel' : 'wgs')
if( params.targets_bed && !params.panel_bed ) params.panel_bed = params.targets_bed
if( !params.containsKey('discard_trimmed_pass') )  params.discard_trimmed_pass = true
    
// define once
if( !params.empty_json ){
  def p = file("${projectDir}/.empty.json"); if( !p.exists() ) p.text = ''
  params.empty_json = p
}
if( !params.empty_tsv ){
  def p = file("${projectDir}/.empty.tsv");  if( !p.exists() ) p.text = ''
  params.empty_tsv = p
}

if( !params.containsKey('report_use_csq') )  params.report_use_csq = true
if( !params.containsKey('report_csq_regex') ) params.report_csq_regex =
  'missense_variant|frameshift_variant|stop_gained|stop_lost|start_lost|' +
  'splice_(acceptor|donor)_variant|inframe_(insertion|deletion)|' +
  'protein_altering_variant|synonymous_variant|' +
  'non_coding_transcript_exon_variant|3_prime_UTR_variant|5_prime_UTR_variant|' +
  'NMD_transcript_variant'

/* Optional UMI/panel QC */
if( !params.containsKey('umi_regex') )      params.umi_regex      = null
if( !params.containsKey('mosdepth_bin') )   params.mosdepth_bin   = 'mosdepth'
if( !params.containsKey('panel_min_cov') )  params.panel_min_cov  = 100
if( !params.containsKey('panel_lod') )      params.panel_lod      = 0.02
if( !params.containsKey('bwa_k') )        params.bwa_k = 20000000   // 20M is sane
if( !params.containsKey('bwa_threads') )  params.bwa_threads = 4    // 2–4 keeps RAM tame
if( !params.containsKey('enable_evidence') ) params.enable_evidence = true
// let users set the token via env or CLI using  --oncokb_token "$YOUR_TOKEN"
if( !params.containsKey('oncokb_token') ) params.oncokb_token = System.getenv('ONCOKB_TOKEN')

/* ---------- container images ---------- */
def IMG_FASTQC   = 'quay.io/biocontainers/fastqc:0.11.9--0'
def IMG_FASTP    = 'quay.io/biocontainers/fastp:0.23.4--h5f740d0_0'
def IMG_BWA      = 'quay.io/biocontainers/bwa:0.7.17--he4a0461_11'
def IMG_SAMTOOLS = 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'
def IMG_BCFTOOLS = 'quay.io/biocontainers/bcftools:1.17--h3cc50cf_1'
def IMG_VEP      = 'ensemblorg/ensembl-vep:release_110.1'
def IMG_JAVA     = 'openjdk:17-jdk-slim'
def IMG_PY       = 'python:3.11'

params.build_html = params.build_html ?: "${baseDir}/HTML/build_html.py"
params.discard_trimmed_pass = params.discard_trimmed_pass ?: false




/* ---------- helpers ---------- */
/* Build (sid, r1, r2) tuples from a file pattern like .../*_{1,2}.fastq.gz */
def readsFromPattern(String pattern) {
  // Build (sid, r1, r2) exactly once
Channel
    .fromFilePairs(pattern, checkIfExists: true)   // no flat:true
    .map { sid, files ->
      // files is usually a List<Path>, but be robust:
      def list = (files instanceof List) ? files : [ files ]
      assert list.size() in [1,2] : "Expected 1 or 2 files for sample ${sid}, got ${list}"
      def r1 = list[0]
      def r2 = (list.size() > 1) ? list[1] : list[0]   // single-end -> reuse R1 as R2
      tuple( sid.toString(), r1, r2 )
    }
}


/* ---------------- workflow ---------------- */
workflow CLINICAL_CORE {
  take:
    reads_ch        // tuples: (sid, r1, r2)
    ref_fa_path     // string or path to FASTA

  main:
    // sanity
    if (!params.ref_fa)     error 'Provide --ref_fa'
    if (!params.reads)      error 'Provide --reads'
    if (!params.patient_id) error 'Provide --patient_id'
    if (params.vep_cache && !file(params.vep_cache).exists())
      error "VEP cache path not found: ${params.vep_cache}"

    // --- refs (each built once, then fanned out) ---
    ref_fa_ch   = Channel.value( file(ref_fa_path as String) )
    def faidx   = REF_FAIDX(ref_fa_ch)
    ref_pair_ch = faidx.ref_pair                  // (ref_fa, ref_fai)

    def idx     = BWA_INDEX(ref_fa_ch)
    bwa_dir_ch  = idx.index_dir                   // e.g. path("bwa_index")

    // QC(raw)
    FASTQC_RAW(reads_ch)

    // trim
    reads_for_trim = params.umi_regex ? UMI_EXTRACT(reads_ch) : reads_ch
    TRIM = TRIM_FASTP(reads_for_trim)       // multi-out
    FASTQC_TRIM(TRIM.reads)                 // use only the reads output

    // align
    ALN_SAM = ALIGN_BWA(TRIM.reads, bwa_dir_ch)
    ALN     = SAMTOOLS_SORT_INDEX(ALN_SAM, ref_pair_ch)

    // mark dups
    def md         = MARKDUPS(ALN, ref_pair_ch)
    def MD_DEFAULT = md.align
    def MD_METRICS = md.metrics

    // reindex (ensure *.crai)
    def MD_FOR_INDEX = MD_DEFAULT.map { sid, cram, crai -> tuple(sid, cram) }
    MD_IDX = REINDEX_CRAM(MD_FOR_INDEX)

    // call variants
    def CALLV       = CALL_VARIANTS_CRAM(MD_IDX, ref_pair_ch)
    def CH_VEP_IN   = CALLV.map { sid, vcf_gz, tbi -> tuple(sid, vcf_gz) }
    def CH_VCF_ANN  = params.vep_cache ? VEP_ANNOTATE(CH_VEP_IN) : CALLV

    // filter for HTML
    def FHTML          = FILTER_FOR_HTML(CH_VCF_ANN)
    def CH_SUMMARY_VCF = FHTML.map { sid, vcf, _tbi -> tuple(sid, vcf) }

    // PharmCAT (optional) -- define ONCE
    PHARM_JSON = params.pharmcat_jar
      ? PHARMCAT(CH_SUMMARY_VCF).pharmcat_json
      : CH_SUMMARY_VCF.map { sid, _ -> tuple(sid, file(params.empty_json)) }

    // knowledge file (singleton)
    KB_FILE = Channel.value( (params.knowledge_tsv ?: params.empty_tsv) as String )

    // ---- Panel QC (optional, brace-free) ----
    def mos           = params.panel_bed ? PANEL_QC_MOSDEPTH(MD_DEFAULT, file(params.panel_bed), ref_pair_ch) : null
def PANEL_JSON_CH = params.panel_bed
    ? QC_AGGREGATE(
        mos.summary,          // (sid, panel_qc.mosdepth.summary.txt)
        MD_METRICS,           // (sid, *.markdup.metrics.txt)
        TRIM.json,            // (sid, *.fastp.json)
        MD_DEFAULT,           // (sid, *.markdup.cram, *.markdup.cram.crai)
        file(params.panel_bed),
        ref_pair_ch,          // (ref_fa, ref_fai)
        mos.thresholds,       // (sid, panel_qc.thresholds.bed.gz)
        mos.regions           // (sid, panel_qc.regions.bed.gz)
      )
    : MD_DEFAULT.map { sid, cram, crai -> tuple(sid, file(params.empty_json)) }


    // ---- Build JSON + mCODE ----
    def SUMMARY_WR = MAKE_JSON_SUMMARY(
      CH_SUMMARY_VCF,
      PHARM_JSON,
      KB_FILE,
      PANEL_JSON_CH
    )

    // ---- HTML report ----
    def builder_ch = Channel.fromPath(params.build_html ?: "${baseDir}/HTML/build_html.py", checkIfExists:true)
    def annotated  = SUMMARY_WR.report_json.map { sid, rep -> tuple(sid, rep) }
    def pdb_file = file( params.pathway_db ?: "${baseDir}/resources/pathway_gene_map_v1.json" )
    def gd_file  = file( params.gene_domains ?: "${baseDir}/resources/gene_domains.json" )

    def pathway_db_ch   = Channel.value( pdb_file.exists() ? pdb_file : file(params.empty_json) )
    def gene_domains_ch = Channel.value( gd_file.exists()  ? gd_file  : file(params.empty_json) )

    def BUILT_HTML = BUILD_HTML(annotated, builder_ch, pathway_db_ch, gene_domains_ch)

  emit:
    aligned     = ALN
    bam_ch      = MD_DEFAULT
    vcf_ch      = CH_VCF_ANN
    panel_ch    = Channel.empty()   // placeholder until wired
    report_json = SUMMARY_WR.report_json
    mcode_json  = SUMMARY_WR.mcode_json
    report_html = BUILT_HTML
}


process VEP_ANNOTATE {
  tag { "vep:${sid}" }
  container "${IMG_VEP}"
  cpus 6
  memory '24 GB'
  time '6h'
  publishDir "${params.outdir}/vep", mode: 'copy', overwrite: true

  // Expect (sid, vcf_gz); if you currently have (sid, vcf_gz, tbi),
  // map it upstream:  CALL.map { sid, v, t -> tuple(sid, v) }
  input:
    tuple val(sid), path(vcf_gz)

  // Return a single tuple carrying sid + both outputs
  output:
    tuple val(sid), path("${sid}.vep.vcf.gz"), path("${sid}.vep.vcf.gz.tbi")

  // Only run if a cache path was provided
  when:
    params.vep_cache != null

  shell:
  '''
  set -euo pipefail

  # Run VEP offline with cache; compress with bgzip; then index with tabix
  vep \
    --offline --cache --dir_cache !{params.vep_cache} \
    --assembly !{params.vep_assembly} --species homo_sapiens \
    --format vcf --vcf --compress_output bgzip \
    --everything --fork !{task.cpus} \
    --input_file !{vcf_gz} \
    --output_file "!{sid}.vep.vcf.gz"

  tabix -p vcf "!{sid}.vep.vcf.gz"
  '''
}

process PANEL_QC_MOSDEPTH {
  tag { "panelqc_${sid}" }
  // run on host so we can use params.mosdepth_bin
  stageInMode 'link'                          // keep original filenames
  cpus 6; memory '24 GB'; time '6h'
  publishDir "${params.outdir}/qc_panel", mode:'copy', overwrite:true

  input:
    tuple val(sid), path(aln), path(aln_index)
    path  panel_bed
    tuple path(ref_fa), path(ref_fai)

  output:
    tuple val(sid), path("panel_qc.mosdepth.summary.txt"), emit: summary
    tuple val(sid), path("panel_qc.regions.bed.gz")      , emit: regions
    tuple val(sid), path("panel_qc.thresholds.bed.gz")   , emit: thresholds
    tuple val(sid), path("panel_qc.json")                , emit: panel_json

  shell:
'''
set -euo pipefail

MOS="!{ params.mosdepth_bin ?: 'mosdepth' }"
SAM="!{ params.samtools_bin ?: 'samtools' }"

# If CRAM, give mosdepth the reference
FA_ARG=""
case "!{aln}" in
  *.cram) FA_ARG="-f !{ref_fa}" ;;
esac

# Ensure FASTA index exists in this dir
[ -s "!{ref_fa}.fai" ] || "$SAM" faidx "!{ref_fa}"

# *** Rebuild the alignment index locally so htslib finds it ***
case "!{aln}" in
  *.cram) "$SAM" index -@ !{task.cpus} -c "!{aln}" ;;   # writes <aln>.crai here
  *.bam)  "$SAM" index -@ !{task.cpus}     "!{aln}" ;;  # writes <aln>.bai
esac

# Quick debug
ls -lh "!{aln}" !{aln}.crai* !{aln}.bai* "!{ref_fa}" "!{ref_fa}.fai" >&2 || true

"$MOS" -t !{task.cpus} ${FA_ARG} \
  --by "!{panel_bed}" \
  --thresholds 20,50,100 \
  panel_qc "!{aln}"

# Normalize summary filename if needed
[ -f panel_qc.mosdepth.summary.txt ] || mv panel_qc*.summary.txt panel_qc.mosdepth.summary.txt

# Tiny JSON digest (use Nextflow interpolation, not ${...})
python3 - << 'PY'
import json, glob
thr = int(!{params.panel_min_cov})
out = {"threshold": thr, "mosdepth": {}}
for s in glob.glob("panel_qc*.summary.txt"):
    for line in open(s):
        p = line.rstrip().split("\t")
        if len(p) >= 4 and p[0] == "total_regions":
            out["mosdepth"]["regions"] = int(p[1])
        if len(p) >= 7 and p[0] == "total_bases":
            out["mosdepth"]["bases"] = int(p[1])
open("panel_qc.json","w").write(json.dumps(out, indent=2))
PY
'''

}

process REF_FAIDX {
  tag { ref_fa.baseName }
  container 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'
  cpus 5
  memory '12 GB'
  time '30m'
  stageInMode 'copy'

  
  input:
    path ref_fa

  output:
    tuple path(ref_fa), path("${ref_fa}.fai"), emit: ref_pair
    path("${ref_fa}.gzi"), optional: true,     emit: gzi

 
  script:
  """
  set -euo pipefail
  echo "PWD: \$(pwd)"
  echo "Listing:"
  ls -l
  echo "Ref path: ${ref_fa}"
  echo "Ref basename: ${ref_fa.baseName}"
  samtools --version

  samtools faidx "${ref_fa}"

  echo 'Result:'
  ls -l "${ref_fa}.fai"
  """
}

process BWA_INDEX {
  tag { reference_fasta.baseName }
  container 'quay.io/biocontainers/bwa:0.7.17--he4a0461_11'
  cpus 4
  memory '8 GB'
  stageInMode 'copy'

  input:
    path reference_fasta

  output:
    path 'bwa_index', emit: index_dir   // contains files with prefix bwa_index/ref

  shell:
  '''
  set -euo pipefail
  mkdir -p bwa_index
  bwa index -p bwa_index/ref "!{reference_fasta}"
  '''
}


process FASTQC_RAW {
  stageInMode 'symlink'
  input: tuple val(sid), path(r1), path(r2)
  output:
    path "${sid}_R1_fastqc.html"
    path "${sid}_R1_fastqc.zip"
    path "${sid}_R2_fastqc.html"
    path "${sid}_R2_fastqc.zip"
  shell:
  '''
  set -euo pipefail
  if [ -s !{r2} ] && [ "$(readlink -f !{r1})" != "$(readlink -f !{r2})" ]; then
    fastqc -q -o ./ !{r1} !{r2}
  else
    fastqc -q -o ./ !{r1}
    : > "!{sid}_R2_fastqc.html"
    : > "!{sid}_R2_fastqc.zip"
  fi
  # rename only the produced reports (tiny), not the inputs
  r1base=$(basename "!{r1}"); s=${r1base%.gz}; s=${s%.fastq}; s=${s%.fq}
  mv "${s}_fastqc.html" "!{sid}_R1_fastqc.html" || true
  mv "${s}_fastqc.zip"  "!{sid}_R1_fastqc.zip"  || true
  if ls *_fastqc.html >/dev/null 2>&1; then
    r2base=$(basename "!{r2}"); s=${r2base%.gz}; s=${s%.fastq}; s=${s%.fq}
    mv "${s}_fastqc.html" "!{sid}_R2_fastqc.html" || true
    mv "${s}_fastqc.zip"  "!{sid}_R2_fastqc.zip"  || true
  fi
  '''
}

process FASTQC_TRIM {
  stageInMode 'symlink'
  input:
    tuple val(sid), path(r1), path(r2)

  output:
    path "${sid}_R1_fastqc.html"
    path "${sid}_R1_fastqc.zip"
    path "${sid}_R2_fastqc.html"
    path "${sid}_R2_fastqc.zip"

  shell:
  '''
  set -euo pipefail

  if [ -s "!{r2}" ] && [ "$(readlink -f "!{r1}")" != "$(readlink -f "!{r2}")" ]; then
    # paired-end
    fastqc -q -o ./ "!{r1}" "!{r2}"

    r1base=$(basename "!{r1}")
    r1stem=${r1base%.gz}; r1stem=${r1stem%.fastq}; r1stem=${r1stem%.fq}

    r2base=$(basename "!{r2}")
    r2stem=${r2base%.gz}; r2stem=${r2stem%.fastq}; r2stem=${r2stem%.fq}

    mv "${r1stem}_fastqc.html" "!{sid}_R1_fastqc.html"
    mv "${r1stem}_fastqc.zip"  "!{sid}_R1_fastqc.zip"
    mv "${r2stem}_fastqc.html" "!{sid}_R2_fastqc.html"
    mv "${r2stem}_fastqc.zip"  "!{sid}_R2_fastqc.zip"
  else
    # single-end
    fastqc -q -o ./ "!{r1}"

    r1base=$(basename "!{r1}")
    r1stem=${r1base%.gz}; r1stem=${r1stem%.fastq}; r1stem=${r1stem%.fq}

    mv "${r1stem}_fastqc.html" "!{sid}_R1_fastqc.html"
    mv "${r1stem}_fastqc.zip"  "!{sid}_R1_fastqc.zip"

    : > "!{sid}_R2_fastqc.html"
    : > "!{sid}_R2_fastqc.zip"
  fi
  '''
}



process TRIM_FASTP {
  tag { sid }
  container IMG_FASTP
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
  '''
  set -euo pipefail

  if [ -s !{r2} ] && [ "$(readlink -f !{r1})" != "$(readlink -f !{r2})" ]; then
    # Paired-end: make two FIFOs using BusyBox-friendly mktemp
    fifo1="$(mktemp -u -p "${TMPDIR:-/tmp}" fastp1.XXXXXX)"; fifo1="${fifo1}.fq.gz"
    fifo2="$(mktemp -u -p "${TMPDIR:-/tmp}" fastp2.XXXXXX)"; fifo2="${fifo2}.fq.gz"
    mkfifo "$fifo1" "$fifo2"

    # Drain FIFOs to /dev/null
    cat "$fifo1" > /dev/null & C1=$!
    cat "$fifo2" > /dev/null & C2=$!

    fastp \
      -i !{r1} -I !{r2} \
      -o "$fifo1" -O "$fifo2" \
      -w !{task.cpus} -q 20 -l 30 \
      -j "!{sid}.fastp.json" \
      -h "!{sid}.fastp.html" \
      2> >(tee "!{sid}.fastp.log" >&2)

    wait $C1 $C2 || true
    rm -f "$fifo1" "$fifo2"

    ln -sf "$(basename "!{r1}")" "!{sid}.R1.trim.fastq.gz"
    ln -sf "$(basename "!{r2}")" "!{sid}.R2.trim.fastq.gz"

  else
    # Single-end: one FIFO
    fifo="$(mktemp -u -p "${TMPDIR:-/tmp}" fastp.XXXXXX)"; fifo="${fifo}.fq.gz"
    mkfifo "$fifo"
    cat "$fifo" > /dev/null & C=$!

    fastp \
      -i !{r1} \
      -o "$fifo" \
      -w !{task.cpus} -q 20 -l 30 \
      -j "!{sid}.fastp.json" \
      -h "!{sid}.fastp.html" \
      2> >(tee "!{sid}.fastp.log" >&2)

    wait $C || true
    rm -f "$fifo"

    ln -sf "$(basename "!{r1}")" "!{sid}.R1.trim.fastq.gz"
    : > "!{sid}.R2.trim.fastq.gz"   # sentinel for SE
  fi
  '''
}




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
  stageInMode 'link'
  publishDir "${params.outdir}/align", mode: 'copy', overwrite: true

  input:
    tuple val(sid), path(sam_gz)
    tuple path(ref_fa), path(ref_fai)

  output:
    tuple val(sid), path("${sid}.sorted.cram"), path("${sid}.sorted.cram.crai")

  shell:
  '''
  set -euo pipefail
  mkdir -p tmp

  # keep per-thread memory modest; adjust if you have RAM headroom
  gzip -cd "!{sam_gz}" \
  | samtools sort -@ !{task.cpus} -m 768M \
      -O CRAM --reference "!{ref_fa}" \
      -T "tmp/!{sid}" -o "!{sid}.sorted.cram" -

  samtools index -@ !{task.cpus} "!{sid}.sorted.cram"
  '''
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
  '''
  set -euo pipefail
  TMPDIR="${TMPDIR:-$PWD/tmp}"; mkdir -p "$TMPDIR"

  gatk --java-options "-Djava.io.tmpdir=$TMPDIR" MarkDuplicates \
    -I "!{cram}" \
    -O "!{sid}.markdup.cram" \
    -M "!{sid}.markdup.metrics.txt" \
    --REFERENCE_SEQUENCE "!{ref_fa}" \
    --CREATE_INDEX true \
    --TMP_DIR "$TMPDIR"

  # Ensure the expected CRAI exists (some environments don’t leave it)
  if [ ! -f "!{sid}.markdup.cram.crai" ]; then
    # Use gatk to build (works for BAM/CRAM; writes the right index type)
    gatk BuildBamIndex -I "!{sid}.markdup.cram" || true
  fi

  # If a .bai were produced for some reason, normalize the name
  if [ ! -f "!{sid}.markdup.cram.crai" ] && [ -f "!{sid}.markdup.cram.bai" ]; then
    ln -sf "!{sid}.markdup.cram.bai" "!{sid}.markdup.cram.crai"
  fi

  # Final sanity (fail clearly if still missing)
  ls -l "!{sid}.markdup.cram" "!{sid}.markdup.cram.crai"
  '''
}

process REINDEX_CRAM {
  tag { sid }
  container 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'
  cpus 3
  time '1h'
  stageInMode 'link'

  input:
    tuple val(sid), path(cram)

  output:
    tuple val(sid), path(cram), path("${cram}.crai")

  shell:
  '''
  set -euo pipefail
  samtools index -@ !{task.cpus} -c "!{cram}"
  '''
}


process CALL_VARIANTS_CRAM {
  tag { sid }
  container "${IMG_BCFTOOLS}"
  cpus 3
  memory '6 GB'
  time '8h'
  stageInMode 'link'
  publishDir "${params.outdir}/vcf", mode:'copy', overwrite:true

  input:
    tuple val(sid), path(cram), path(crai)
    tuple path(ref_fa), path(ref_fai)

  output:
    tuple val(sid), path("${sid}.vcf.gz"), path("${sid}.vcf.gz.tbi")

  shell:
  '''
  set -euo pipefail
  set -x

  # Show inputs
  ls -lh !{cram} !{crai} !{ref_fa} !{ref_fai} >&2 || true

  # Ensure the exact filenames htslib looks for exist in THIS dir.
  # Only copy if the existing path isn't the same inode.
  if [ ! -e "!{cram}.crai" ] || [ "$(readlink -f "!{cram}.crai" 2>/dev/null || echo NO)" != "$(readlink -f "!{crai}" 2>/dev/null || echo NO)" ]; then
    cp -fL "!{crai}" "!{cram}.crai"
  fi
  if [ ! -e "!{ref_fa}.fai" ] || [ "$(readlink -f "!{ref_fa}.fai" 2>/dev/null || echo NO)" != "$(readlink -f "!{ref_fai}" 2>/dev/null || echo NO)" ]; then
    cp -fL "!{ref_fai}" "!{ref_fa}.fai"
  fi

  # Optional args from params
  ''' +
  (params.region    ? "ARGS_R='-r ${params.region}'\n"    : "ARGS_R=\n") +
  (params.panel_bed ? "ARGS_T='-R ${params.panel_bed}'\n" : "ARGS_T=\n") +
  '''
  # mpileup -> BCF (log stderr, fail loudly if empty)
  bcftools mpileup -Ob -f "!{ref_fa}" ${ARGS_R:+$ARGS_R} ${ARGS_T:+$ARGS_T} \
      -q 20 -Q 20 -a AD,DP,SP "!{cram}" -o mpileup.bcf 2> mpileup.log
  rc=$?
  if [ $rc -ne 0 ] || [ ! -s mpileup.bcf ]; then
    echo "---- mpileup.log ----" >&2
    sed -n '1,200p' mpileup.log >&2 || true
    exit ${rc:-1}
  fi

  # call + index (also surface errors)
  bcftools call -mv -Oz -o "!{sid}.vcf.gz" mpileup.bcf 2> call.log || {
    echo "---- call.log ----" >&2
    sed -n '1,200p' call.log >&2 || true
    exit 1
  }
  bcftools index -t "!{sid}.vcf.gz" 2> index.log || {
    echo "---- index.log ----" >&2
    sed -n '1,200p' index.log >&2 || true
    exit 1
  }
  '''
}




process QC_AGGREGATE {
  tag { sid }
  container null
  cpus 4; memory '4 GB'; time '2h'
  publishDir "${params.outdir}/qc_panel", mode:'copy', overwrite:true

  input:
    tuple val(sid), path(summary_txt)
    tuple val(sid), path(markdup_metrics)
    tuple val(sid), path(fastp_json)
    tuple val(sid), path(cram), path(crai)
    path  panel_bed
    tuple path(ref_fa), path(ref_fai)
    tuple val(sid), path(thresholds_bedgz)
    tuple val(sid), path(regions_bedgz)

  output:
    tuple val(sid), path("panel_qc.json")

  shell:
  '''
  set -euo pipefail
  SAM="!{ params.samtools_bin ?: 'samtools' }"

  # total mapped via idxstats (fast; index-based)
  total_mapped=$("$SAM" idxstats "!{cram}" | awk '{m+=$3} END{print (m?m:0)}' || echo 0)
  # on-target mapped (primary only), multi-threaded
  ont=$("$SAM" view -@ !{task.cpus} -c -F 260 -L "!{panel_bed}" "!{cram}" || echo 0)

  on_target_rate=$(awk -v t="$total_mapped" -v o="$ont" 'BEGIN{ if (t>0) printf "%.1f", (o/t)*100; else print "0.0" }')
  export ON_TARGET_RATE="$on_target_rate"

  python3 - << 'PY'
import gzip, json, statistics as stats, os

summary = "!{summary_txt}"
markdup = "!{markdup_metrics}"
fastp   = "!{fastp_json}"
thresh  = "!{thresholds_bedgz}"
regions = "!{regions_bedgz}"

def _to_int(x):
    try: return int(x)
    except: return None

def _to_float(x):
    try: return float(x)
    except: return None

# --- mean coverage from regions (weighted) ---
total_len = 0
sum_cov   = 0.0
region_means = []
with gzip.open(regions, "rt") as f:
    for line in f:
        p = line.rstrip().split("\t")
        if len(p) < 4:
            continue
        start, end = _to_int(p[1]), _to_int(p[2])
        if start is None or end is None:
            continue
        length = max(0, end - start)

        # mosdepth --by may add a region name column; mean is then at p[4]
        mean = _to_float(p[3])
        if mean is None and len(p) > 4:
            mean = _to_float(p[4])
        if mean is None:
            continue

        total_len += length
        sum_cov   += mean * length
        region_means.append(mean)

mean_cov   = round(sum_cov/total_len, 1) if total_len else None
median_cov = round(stats.median(region_means), 1) if region_means else None


# --- breadth from thresholds ---
# We now run mosdepth with --thresholds 20,50,100 (in that order).
ge = {20:0, 50:0, 100:0}
bases_total = 0

with gzip.open(thresh, "rt") as f:
    for line in f:
        p = line.rstrip().split("\t")
        if len(p) < 4:
            continue

        start, end = _to_int(p[1]), _to_int(p[2])
        if start is None or end is None:
            continue
        L = max(0, end - start)

        # counts start at p[3] if numeric, otherwise skip a name column (start at p[4])
        first_idx = 3
        if _to_int(p[3]) is None and len(p) > 4 and _to_int(p[4]) is not None:
            first_idx = 4

        vals = [ _to_int(x) for x in p[first_idx:first_idx+3] ]  # 20,50,100
        if L > 0:
            bases_total += L
            if len(vals) > 0 and vals[0] is not None: ge[20]  += vals[0]
            if len(vals) > 1 and vals[1] is not None: ge[50]  += vals[1]
            if len(vals) > 2 and vals[2] is not None: ge[100] += vals[2]

def pct(x):
    return round(100.0 * x / bases_total, 1) if bases_total and x is not None else None

ge20, ge50, ge100 = pct(ge[20]), pct(ge[50]), pct(ge[100])


# --- duplication from MarkDuplicates metrics ---
dup_pct = None
try:
    with open(markdup) as fh:
        for line in fh:
            if line.startswith("LIBRARY"):
                header = line.strip().split("\t")
                idx = header.index("PERCENT_DUPLICATION")
                row = next(fh).strip().split("\t")
                dup_pct = round(float(row[idx]) * 100.0, 1)
                break
except Exception:
    pass

# --- Q30 from fastp ---
q30 = None
try:
    j = json.load(open(fastp))
    q30 = round(float(j["summary"]["after_filtering"]["q30_rate"]) * 100.0, 1)
except Exception:
    pass

on_target_rate = float(os.getenv("ON_TARGET_RATE", "0.0"))

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
  '''
}









/* ---------- Optional annotation + clinical ---------- */




process PHARMCAT {
  tag { "pharmcat:${sid}" }
  // Pin a recent tag so behavior is stable across runs
  container "pgkb/pharmcat:3.0.1"
  cpus 4; memory '8 GB'; time '3h'
  publishDir "${params.outdir}/pharmcat", mode: 'copy', overwrite: true

  input:
    tuple val(sid), path(vcf)   // (sid, vcf[.gz/.bgz]) from CALL_VARIANTS or VEP

  output:
    tuple val(sid), path("${sid}.pharmcat.report.html"), emit: pharmcat_html
    tuple val(sid), path("${sid}.pharmcat.report.json"), emit: pharmcat_json

  script:
  // The PharmCAT pipeline accepts plain, gz, or bgz VCFs.
  // We call the wrapper and explicitly request both HTML+JSON (stable on 3.x).
  """
  set -euo pipefail

  IN="input.vcf"
  case "!{vcf}" in
    *.gz|*.bgz)  gunzip -c "!{vcf}" > "\$IN" ;;
    *)           cp "!{vcf}" "\$IN" ;;
  esac

  # Run PharmCAT (Docker image bundles all deps)
  pharmcat_pipeline "\$IN" \
    -bf "!{sid}" \
    -reporterHtml -reporterJson

  # Normalize filenames for downstream
  mv "!{sid}.report.html"  "!{sid}.pharmcat.report.html"
  mv "!{sid}.report.json"  "!{sid}.pharmcat.report.json" || \
    cp "!{sid}".*.json     "!{sid}.pharmcat.report.json" 2>/dev/null || true
  """
}


process MAKE_JSON_SUMMARY {
  tag "summary_json"
  container "${IMG_PY}"
  cpus 6; memory '24 GB'; time '6h'
  publishDir "${params.outdir}/json", mode: 'copy', overwrite: true

  input:
    tuple val(sid), path(vep_vcf,       stageAs: 'vcf_in')
    tuple val(sid), path(pharmcat_json, stageAs: 'pharmcat.json')
    path kb_tsv,        stageAs: 'knowledge.tsv'   // broadcast singleton
    tuple val(sid), path(panel_qc_json, stageAs: 'panel_qc.json')  // <-- per-sample

  output:
    tuple val(sid), path("${sid}.report.json"),       emit: report_json
    tuple val(sid), path("${sid}.mcode_bundle.json"), emit: mcode_json

  shell:
  '''
  set -euo pipefail

  # pass knobs via environment (bash), then start Python
  export REPORT_CODING_ONLY='!{ params.report_coding_only ? "1" : "0" }'
  export REPORT_KEEP_IMPACT='!{ params.report_keep_impact ?: "HIGH,MODERATE" }'
  export REPORT_MAX_PER_GENE='!{ (params.report_max_per_gene ?: 50).toString() }'
  export REPORT_MAX_ROWS='!{ (params.report_max_rows ?: 100000).toString() }'
  export REPORT_SEARCH_LIMIT='!{ (params.report_search_limit ?: 0).toString() }'
  # intentionally NO export REPORT_EXPORT_TSV here (defaults to '1' in Python)
  export ENABLE_EVIDENCE='!{ params.enable_evidence ? "1" : "0" }'
  export ONCOKB_TOKEN='!{ params.oncokb_token ?: "" }'

  python3 - << 'PY'
import os, json, csv, gzip, re, unicodedata
from pathlib import Path

# -------- tunable report filters (env overrides allowed) --------
MIN_QUAL     = float(os.getenv('REPORT_MIN_QUAL', '30'))
MIN_DP       = int(os.getenv('REPORT_MIN_DP',   '10'))
REQUIRE_PASS = os.getenv('REPORT_REQUIRE_PASS', '1') == '1'
APPLY_AB     = os.getenv('REPORT_APPLY_AB',     '0') == '1'
AB_MIN       = float(os.getenv('REPORT_AB_MIN', '0.25'))
AB_MAX       = float(os.getenv('REPORT_AB_MAX', '0.75'))
CODING_ONLY  = os.getenv('REPORT_CODING_ONLY',  '0') == '1'
KEEP_IMPACT   = {s.strip().upper() for s in os.getenv('REPORT_KEEP_IMPACT','').split(',') if s.strip()}
MAX_PER_GENE  = int(os.getenv('REPORT_MAX_PER_GENE','0') or '0')
MAX_ROWS      = int(os.getenv('REPORT_MAX_ROWS','0') or '0')
SEARCH_LIMIT  = int(os.getenv('REPORT_SEARCH_LIMIT','0') or '0')
EXPORT_TSV    = os.getenv('REPORT_EXPORT_TSV','1') == '1'  # default ON

sid        = "!{sid}"
vcf_path   = Path("vcf_in")          # staged VCF (possibly bgzip)
pharm_path = Path("pharmcat.json")
kb_path    = Path("knowledge.tsv")
panel_path = Path("panel_qc.json")

def exists_file(p: Path) -> bool:
    try:    return p.exists() and p.is_file() and p.stat().st_size > 0
    except: return False

def open_vcf(p: Path):
    # gzip magic check (works regardless of extension)
    with p.open('rb') as fh:
        head = fh.read(2)
    if head == b'\\x1f\\x8b':
        return gzip.open(p, 'rt')
    return p.open('rt')

def norm(s): return unicodedata.normalize('NFKC', s) if isinstance(s,str) else s

# ---------- optional PharmCAT ----------
phc_obj = None
if exists_file(pharm_path):
    try:    phc_obj = json.loads(pharm_path.read_text())
    except: phc_obj = {"_error":"Failed to parse PharmCAT JSON"}

# ---------- optional knowledge TSV ----------
kb = {}
if exists_file(kb_path):
    with kb_path.open(newline='') as f:
        for r in csv.DictReader(f, delimiter='\\t'):
            g  = norm((r.get('gene') or '').upper())
            vk = norm(r.get('variant_key') or '')
            kb.setdefault((g, vk), []).append({
                "therapy":  r.get('therapy')  or '',
                "evidence": r.get('evidence') or '',
                "source":   r.get('source')   or '',
                "disease":  r.get('disease')  or '',
                "tier":     r.get('tier')     or '',
                "notes":    r.get('notes')    or ''
            })

# ---------- optional panel QC ----------
qc_panel = None
qc_summary_block = None
if exists_file(panel_path):
    try:    qc_panel = json.loads(panel_path.read_text())
    except: qc_panel = {"_error":"Failed to parse panel_qc.json"}
if qc_panel and isinstance(qc_panel, dict):
    qc_summary_block = {
        "Mean Coverage":      f"{qc_panel.get('mean_coverage')}×" if qc_panel.get('mean_coverage') is not None else None,
        "Median Coverage":    f"{qc_panel.get('median_coverage')}×" if qc_panel.get('median_coverage') is not None else None,
        "≥20×":               f"{qc_panel.get('ge20_pct')}%" if qc_panel.get('ge20_pct') is not None else None,
        "≥50×":               f"{qc_panel.get('ge50_pct')}%" if qc_panel.get('ge50_pct') is not None else None,
        "≥100×":              f"{qc_panel.get('ge100_pct')}%" if qc_panel.get('ge100_pct') is not None else None,
        "On-target Rate":     f"{qc_panel.get('on_target_rate_pct')}%" if qc_panel.get('on_target_rate_pct') is not None else None,
        "Contamination":      f"{qc_panel.get('contamination_pct')}%" if qc_panel.get('contamination_pct') is not None else None,
        "Duplication":        f"{qc_panel.get('duplication_pct')}%" if qc_panel.get('duplication_pct') is not None else None,
        "Insert Size (median)": (f"{qc_panel.get('insert_size_median_bp')} bp" if qc_panel.get('insert_size_median_bp') is not None else None),
        "Q30":                f"{qc_panel.get('q30_pct')}%" if qc_panel.get('q30_pct') is not None else None
    }

# ---------- helpers ----------
LOF_TERMS = {"frameshift_variant","stop_gained","splice_acceptor_variant",
             "splice_donor_variant","start_lost","stop_lost","transcript_ablation"}
BRCA_SET  = {"BRCA1","BRCA2"}

CODING_TERMS = {
  "missense_variant","stop_gained","stop_lost","start_lost","frameshift_variant",
  "inframe_insertion","inframe_deletion","splice_acceptor_variant","splice_donor_variant",
  "protein_altering_variant"
}

def is_lof_annotation(a):
    cons = a.get('Consequence','') or ''
    if (a.get('LoF') or '').upper() in {'HC','LC'}: return True
    if any(term in cons for term in LOF_TERMS):     return True
    hgvsp = a.get('HGVSp') or a.get('HGVSp_VEP') or ''
    return bool(re.search(r'fs', hgvsp, re.I) or re.search(r'(Ter|\\*)\\d*$', hgvsp))

def kb_lookup(gene, hgvsp, all_ann, alt, info):
    res=[]; gU=(gene or '').upper()
    if gene and hgvsp: res.extend(kb.get((gU, f"{gene}:{hgvsp}"), []))
    if gU in BRCA_SET and any(is_lof_annotation(a) for a in all_ann):
        for vk in (f"{gene}:LoF", f"{gene}:p.* (LoF)", f"{gene}:loss_of_function"):
            res.extend(kb.get((gU, vk), []))
    is_sv = ('[' in alt or ']' in alt or 'SVTYPE=' in info)
    if is_sv:
        genes = sorted({ (a.get('SYMBOL') or '') for a in all_ann if a.get('SYMBOL') })
        for a_gene in genes:
            for b_gene in genes:
                if a_gene != b_gene:
                    for vk in (f"{a_gene}-{b_gene}", f"{a_gene}:FUSION"):
                        res.extend(kb.get((a_gene.upper(), vk), []))
    return res

def parse_info(s):
    d={}
    for item in (s or '').split(';'):
        if not item: continue
        if '=' in item:
            k,v=item.split('=',1); d[k]=v
        else:
            d[item]=True
    return d

def parse_fmt_sample(fmt, sample):
    m={}
    if not fmt or not sample: return m
    ks=fmt.split(':'); vs=sample.split(':')
    for i,k in enumerate(ks):
        if i < len(vs): m[k]=vs[i]
    return m

def first_allelic_balance(ad_str):
    # AD like "ref,alt1[,alt2...]"; use first ALT vs REF
    if not ad_str: return None
    try:
        parts=[int(x) for x in ad_str.split(',') if x!='.']
        if len(parts) >= 2:
            ref,alt = parts[0], parts[1]
            denom = ref + alt
            return (alt/denom) if denom>0 else None
    except: pass
    return None

def passes_filters(q, filt, dp, gt, ab, anns):
    if REQUIRE_PASS and (filt or '').upper() not in ('PASS','.'): return False
    if q is None or q < MIN_QUAL: return False
    if dp is None or dp < MIN_DP: return False
    if APPLY_AB and gt in {'0/1','1/0'} and ab is not None:
        if not (AB_MIN <= ab <= AB_MAX): return False
    if CODING_ONLY:
        # keep only if at least one annotation contains a coding-ish consequence
        cons_all = ' '.join([a.get('Consequence','') or '' for a in anns])
        if not any(term in cons_all for term in CODING_TERMS): return False
    if KEEP_IMPACT:
        # keep only if any annotation's IMPACT is in the allowed set
        imp_vals = {(a.get('IMPACT') or '').upper() for a in anns if a.get('IMPACT')}
        if not imp_vals.intersection(KEEP_IMPACT): return False

    return True

def tier_of(v):
    tiers=[(s.get('tier') or '') for s in (v.get('treatments') or [])]
    if any('T1' in t for t in tiers): return 'T1'
    if any('T2' in t for t in tiers): return 'T2'
    return 'T3+'

from urllib.parse import quote_plus

def search_links(v):
    g = (v.get('gene') or '').strip()
    p = (v.get('hgvsp') or '').strip()
    cons = (v.get('consequence') or '').strip()

    terms = [t for t in [g, (p or cons), 'therapy', 'label OR guideline OR trial', '"2019..2025"'] if t]
    goog  = "https://www.google.com/search?q=" + "+".join(map(quote_plus, terms))
    pubm  = "https://pubmed.ncbi.nlm.nih.gov/?term=" + quote_plus(" ".join([t for t in [g, (p or cons), 'therapy'] if t])) + "&filter=years.2019-2025"

    oncokb = f"https://www.oncokb.org/search?query={quote_plus(g+' '+p)}" if g and p else (f"https://www.oncokb.org/gene/{quote_plus(g)}" if g else None)
    civic  = f"https://civicdb.org/search/variants?q={quote_plus((g+' '+p).strip() or g)}" if g else None
    trials = "https://clinicaltrials.gov/search?cond=&term=" + quote_plus(" ".join([t for t in [g, (p or cons)] if t]))

    return {"google": goog, "pubmed": pubm, "oncokb": oncokb, "civic": civic, "trials": trials}

# ---- consequence classification helpers (must be defined before use) ----
VEP_PRIORITY = [
  "transcript_ablation","splice_acceptor_variant","splice_donor_variant",
  "stop_gained","frameshift_variant","stop_lost","start_lost",
  "inframe_insertion","inframe_deletion","missense_variant","protein_altering_variant",
  "synonymous_variant","coding_sequence_variant",
  "5_prime_UTR_variant","3_prime_UTR_variant",
  "non_coding_transcript_exon_variant","non_coding_transcript_variant",
  "intron_variant","upstream_gene_variant","downstream_gene_variant",
  "regulatory_region_variant","TF_binding_site_variant","intergenic_variant"
]



TYPE_MAP = {
  "missense_variant": "missense",
  "frameshift_variant": "frameshift",
  "stop_gained": "nonsense",
  "stop_lost": "stop_lost",
  "start_lost": "start_lost",
  "inframe_insertion": "inframe_ins",
  "inframe_deletion": "inframe_del",
  "protein_altering_variant": "protein_altering",
  "synonymous_variant": "synonymous",
  "splice_acceptor_variant": "splice_acceptor",
  "splice_donor_variant": "splice_donor",
  "5_prime_UTR_variant": "UTR5",
  "3_prime_UTR_variant": "UTR3",
  "non_coding_transcript_exon_variant": "ncRNA_exon",
  "non_coding_transcript_variant": "ncRNA",
  "intron_variant": "intron",
  "upstream_gene_variant": "upstream",
  "downstream_gene_variant": "downstream",
  "regulatory_region_variant": "regulatory",
  "TF_binding_site_variant": "regulatory",
  "intergenic_variant": "intergenic",
}

def _terms_from_cons(cons: str):
    for t in (cons or "").replace(",", "&").split("&"):
        t = t.strip()
        if t:
            yield t

def _prefer_canonical(anns):
    canon = [a for a in anns if (a.get("CANONICAL") or "").upper() == "YES"]
    return canon if canon else anns

def friendly_type(anns):
    idx = {t:i for i,t in enumerate(VEP_PRIORITY)}
    best = None; best_i = 10**9
    for a in _prefer_canonical(anns):
        for t in _terms_from_cons(a.get("Consequence") or ""):
            i = idx.get(t, 10**6)
            if i < best_i:
                best_i, best = i, t
    return TYPE_MAP.get(best, best or "other")

def is_coding_like(anns):
    if any((a.get("BIOTYPE") or "") == "protein_coding" for a in anns):
        return True
    codingish = {
        "missense_variant","frameshift_variant","stop_gained","stop_lost","start_lost",
        "inframe_insertion","inframe_deletion","protein_altering_variant","synonymous_variant"
    }
    return any(t in codingish for a in anns for t in _terms_from_cons(a.get("Consequence") or ""))
# --- add right after the current imports ---
import urllib.request, urllib.error
from urllib.parse import quote_plus

# ----- Online evidence (CIViC + optional OncoKB) -----
def _http_json_post(url, payload, headers=None, timeout=25):
    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(url, data=data,
                                 headers={"Content-Type": "application/json", **(headers or {})})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode("utf-8"))

def civic_query_variants(gene: str, aa_change: str):
    """Query CIViC v2 GraphQL for predictive evidence for gene+AA change."""
    if not gene:
        return None
    term = f"{gene} {aa_change}".strip()
    q = """
    query VAR($gene: String!, $term: String!) {
      search(query: $term, entityType: VARIANTS, genes: [$gene], first: 10) {
        edges {
          node {
            ... on Variant {
              id
              name
              evidenceItems(evidenceType: PREDICTIVE, first: 20) {
                nodes {
                  clinicalSignificance
                  disease { name }
                  therapies { name }
                  evidenceLevel
                  source { url }
                }
              }
            }
          }
        }
      }
    }"""
    try:
        out = _http_json_post("https://civicdb.org/api/graphql",
                              {"query": q, "variables": {"gene": gene, "term": term}})
        edges = (((out or {}).get("data") or {}).get("search") or {}).get("edges") or []
        nodes = []
        for e in edges:
            v = (e.get("node") or {})
            nodes += (((v.get("evidenceItems") or {}).get("nodes")) or [])
        return {"evidenceItems": {"nodes": nodes}}
    except Exception:
        return None

def oncokb_by_protein(gene: str, hgvsp: str, token: str):
    """OncoKB annotate API by protein change. Token optional; skip if missing."""
    aa = (hgvsp or "").replace("p.", "").replace("P.", "")
    if not gene or not aa or not token:
        return None
    url = ("https://www.oncokb.org/api/v1/annotate/mutations/byProteinChange"
           f"?hugoSymbol={quote_plus(gene)}&alteration={quote_plus(aa)}")
    # Try Bearer first, fall back to X-API-KEY (covers old/new deployments)
    for hdr in ({"Authorization": f"Bearer {token}", "Accept": "application/json"},
                {"X-API-KEY": token, "Accept": "application/json"}):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=hdr), timeout=25) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception:
            pass
    return None

# put near the other helpers in MAKE_JSON_SUMMARY
AA3 = {
  'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Gln':'Q','Glu':'E','Gly':'G',
  'His':'H','Ile':'I','Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P','Ser':'S',
  'Thr':'T','Trp':'W','Tyr':'Y','Val':'V','Sec':'U','Pyl':'O'
}
import re

def _aa_from_hgvsp(hgvsp: str):
    """Return both formats: ('C13*', 'Cys13Ter') when possible."""
    s = hgvsp or ''
    # keep only the part after any ":" and the leading 'p.'
    m = re.search(r':\s*[pP]\.(.+)$', s)
    core = m.group(1) if m else re.sub(r'^[pP]\.', '', s)

    # 3-letter → 1-letter (Cys13Ter -> C13*)
    m3 = re.match(r'([A-Za-z]{3})(\d+)([A-Za-z]{3}|Ter|\*)$', core)
    if m3:
        a1 = AA3.get(m3.group(1).title(), m3.group(1)[0].upper())
        a2raw = m3.group(3)
        a2 = '*' if a2raw in ('Ter','*') else AA3.get(a2raw.title(), a2raw[0].upper())
        return (f"{a1}{m3.group(2)}{a2}", f"{m3.group(1).title()}{m3.group(2)}{a2raw}")

    # already 1-letter? (e.g., C13*, R140H)
    m1 = re.match(r'([A-Z\*])(\d+)([A-Z\*])$', core)
    if m1:
        return (''.join(m1.groups()), core)

    # fallback: return core twice
    return (core, core)


def enrich_rows_with_online(report_dict: dict):
    token = (os.getenv("ONCOKB_TOKEN") or "").strip()

    # (optional) only enrich a limited, high-value subset
    LIMIT = int(os.getenv("REPORT_SEARCH_LIMIT", "200") or "0")
    candidates = [r for r in report_dict.get("rows", [])
                  if (r.get("coding") and (r.get("impact","").upper() in {"HIGH","MODERATE"}))]
    rows_to_enrich = candidates[:LIMIT] if LIMIT else candidates

    for r in rows_to_enrich:
        g = (r.get("gene") or "").strip()
        p = (r.get("hgvsp") or "").strip()
        if not g or not p:
            continue

        aa1, aa3 = _aa_from_hgvsp(p)  # <- use the new normalizer

        # CIViC: try 1-letter first, then 3-letter
        civ = civic_query_variants(g, aa1) or civic_query_variants(g, aa3)

        # OncoKB: use 1-letter form only
        okb = oncokb_by_protein(g, aa1, token) if token else None

        if civ:
            items = ((civ.get("evidenceItems") or {}).get("nodes")) or []
            treats = []
            assoc  = r.get("association")
            for e in items:
                cs  = (e.get("clinicalSignificance") or "").strip()
                dis = ((e.get("disease") or {}).get("name")) or ""
                lvl = (e.get("evidenceLevel") or "").strip()
                url = ((e.get("source") or {}).get("url")) or ""
                for t in (e.get("therapies") or []):
                    nm = t.get("name")
                    if nm:
                        treats.append({"therapy": nm, "evidence": f"{cs or 'Association'} · CIViC {lvl}".strip(" ·"), "source_url": url})
                if not assoc and (cs or dis):
                    assoc = f"{cs or 'Association'} — {dis}".strip(" —")
            if treats:
                r.setdefault("treatments", []).extend(treats)
            if assoc and not r.get("association"):
                r["association"] = assoc
            if items and not r.get("known_label"):
                r["known_label"] = "Known"

        if okb:
            if okb.get("hotspot"):
                r["known_label"] = "Hotspot"
            if okb.get("oncogenic"):
                r["oncogenic"] = okb.get("oncogenic")
            for tx in (okb.get("treatments") or []):
                for d in (tx.get("drugs") or []):
                    nm = d.get("drugName") or d.get("name")
                    if nm:
                        r.setdefault("treatments", []).append({
                            "therapy": nm,
                            "evidence": f"OncoKB {tx.get('level','')}".strip(),
                            "source_url": "https://www.oncokb.org/"
                        })

    # refresh tiers after enrichment
    for r in report_dict.get("rows", []):
        r["tier"] = tier_of(r)


# ---------- stream VCF and build filtered rows ----------
rows = []
csq_fields = []

# ---- stream VCF and build filtered rows ----
rows = []
csq_fields = []
TAB = chr(9)  # literal tab, avoids '\t' getting mangled by editors

with open_vcf(vcf_path) as fh:
    for line in fh:
        if line.startswith('##INFO=<ID=CSQ'):
            # robustly capture the field order from the header
            m = re.search(r'Format:\s*([^">]+)', line)
            if m:
                csq_fields = [f.strip() for f in m.group(1).split('|')]
            continue
        if line.startswith('#'):
            continue

        # Trim newline(s) safely; then split by TAB (not spaces!)
        parts = line.rstrip().split(TAB)
        if len(parts) < 8:
            continue

        chrom, pos, _id, ref, alt, qual_s, filt, info = parts[:8]
        fmt  = parts[8] if len(parts) > 8 else ''
        samp = parts[9] if len(parts) > 9 else ''

        info_d = parse_info(info)
        fmt_d  = parse_fmt_sample(fmt, samp)

        # numeric fields
        try:
            q = None if qual_s == '.' else float(qual_s)
        except Exception:
            q = None

        # depth: prefer sample DP, else INFO/DP, else AD sum
        dp = None
        if fmt_d.get('DP') and fmt_d['DP'] != '.':
            try: dp = int(fmt_d['DP'])
            except Exception: dp = None
        if dp is None and info_d.get('DP') and info_d['DP'] != '.':
            try: dp = int(info_d['DP'])
            except Exception: dp = None
        if dp is None and fmt_d.get('AD'):
            try:
                nums = [int(x) for x in fmt_d['AD'].split(',') if x != '.']
                dp = sum(nums) if nums else None
            except Exception:
                dp = None

        gt = fmt_d.get('GT') or ''
        ab = first_allelic_balance(fmt_d.get('AD'))

        # ----- metrics (inside the loop) -----
        AD_list = None
        if fmt_d.get('AD'):
            try:
                AD_list = [int(x) for x in fmt_d['AD'].split(',') if x != '.']
            except Exception:
                AD_list = None

        AF = None
        if 'AF' in info_d and info_d['AF'] not in ('.', ''):
            try: AF = float(info_d['AF'].split(',')[0])
            except Exception: AF = None
        elif 'AC' in info_d and 'AN' in info_d:
            try:
                ac0 = int(info_d['AC'].split(',')[0]); an = int(info_d['AN'])
                AF = (ac0 / an) if an else None
            except Exception:
                AF = None

        VAF = None
        if AD_list and len(AD_list) >= 2:
            denom = AD_list[0] + AD_list[1]
            VAF = (AD_list[1] / denom) if denom > 0 else None

        # parse CSQ annotations
        anns = []
        csq_val = None
        for it in info.split(';'):
            if it.startswith('CSQ='):
                csq_val = it[4:]
                break
        if csq_val and csq_fields:
            for rec in csq_val.split(','):
                vals = rec.split('|')
                anns.append({k: (vals[i] if i < len(vals) else '') for i, k in enumerate(csq_fields)})

        # basic gene/cons/hgvsp pick
        gene = ''; cons = ''; hgvsp = ''
        if anns:
            a0 = anns[0]
            gene  = a0.get('SYMBOL', '') or ''
            cons  = a0.get('Consequence', '') or ''
            hgvsp = a0.get('HGVSp', '') or a0.get('HGVSp_VEP', '') or ''

        # apply report filters
        if not passes_filters(q, filt, dp, gt, ab, anns):
            continue

        # build the row dict
        e = {
            "chrom": chrom, "pos": int(pos), "ref": ref, "alt": alt,
            "gene": gene, "consequence": cons, "hgvsp": hgvsp,
            "qual": q, "filter": filt, "depth": dp
        }
        if ab is not None: e["ab"] = round(ab, 4)
        if AF is not None: e["af"] = round(AF, 6)

        metrics = {}
        if dp is not None:  metrics["DP"] = dp
        if AD_list:         metrics["AD"] = AD_list[:2]
        if VAF is not None: metrics["VAF"] = round(VAF, 6)
        if q is not None:   metrics["QUAL"] = q
        if metrics:         e["metrics"] = metrics

        # knowledge lookups (and quick columns)
        tre = kb_lookup(gene, hgvsp, anns, alt, info)
        if tre:
            first = tre[0]
            if (first.get("disease") or "") and not e.get("cancer_type"):
                e["cancer_type"] = first["disease"]
            if (first.get("evidence") or "") and not e.get("association"):
                e["association"] = first["evidence"]
            e["treatments"] = tre

        # classification
        e["type"]   = friendly_type(anns)
        e["coding"] = is_coding_like(anns)
        imp = (anns[0].get("IMPACT") if anns else "") or ""
        if imp: e["impact"] = imp

        # tier + links
        e["tier"] = tier_of({"treatments": e.get("treatments")})
        if (e["tier"] or "").startswith('T3'):
            e["search"] = search_links(e)

        # per-gene cap (needs stable local dict)
        if MAX_PER_GENE:
            _g = (gene or 'NA')
            _cnt = locals().setdefault("_per_gene_counts", {}).get(_g, 0)
            if _cnt >= MAX_PER_GENE:
                continue

        rows.append(e)
        if MAX_PER_GENE:
            locals()["_per_gene_counts"][_g] = locals()["_per_gene_counts"].get(_g, 0) + 1

        # global max rows
        if MAX_ROWS and len(rows) >= MAX_ROWS:
            break


# ---------- end VCF loop ----------
IMPACT_RANK = {'HIGH':0,'MODERATE':1,'LOW':2,'MODIFIER':3,'':4}
TYPE_RANK = {
  'nonsense':0,'frameshift':1,'splice_acceptor':2,'splice_donor':3,
  'start_lost':4,'stop_lost':5,'missense':6,'inframe_del':7,'inframe_ins':8,
  'protein_altering':9,'synonymous':10
}
def _v(r,k): return (r.get('metrics') or {}).get(k, 0.0)
def sort_key(r):
    imp = (r.get('impact') or '').upper()
    typ = r.get('type') or ''
    tier_boost = 0 if (r.get('treatments') or []) else 1  # KB hit first if present
    return (
        tier_boost,
        IMPACT_RANK.get(imp, 4),
        TYPE_RANK.get(typ, 99),
        -_v(r,'VAF'),
        -_v(r,'QUAL')
    )
rows = sorted(rows, key=sort_key)

TSG = {'TP53','BRCA1','BRCA2','PTEN','RB1','ATM','CDKN2A','NF1','APC','SMAD4'}
LOF_TYPES = {'nonsense','frameshift','splice_acceptor','splice_donor','start_lost','stop_lost'}

def _pct(x): 
    return (f"{round(100*x,1)}%" if isinstance(x,(int,float)) and x is not None else None)

def make_comment(r):
    g   = r.get('gene') or ''
    p   = r.get('hgvsp') or r.get('type') or ''
    imp = (r.get('impact') or '').upper()
    typ = r.get('type') or ''
    vaf = (r.get('metrics') or {}).get('VAF')
    bits = []

    if imp in {'HIGH','MODERATE'} and r.get('coding'):
        bits.append(f"{g} {p}: {imp.lower()}-impact {typ.replace('_',' ')} variant")
    if g in TSG and typ in LOF_TYPES:
        bits.append("pattern consistent with loss-of-function in a tumor suppressor")
    pvaf = _pct(vaf)
    if pvaf: bits.append(f"VAF ~{pvaf}")
    if r.get('treatments'):
        bits.append("KB evidence present")

    return "; ".join(bits) + "."



# ---------- simple analytics from rows ----------
import statistics as stats
def _safe_median(vals):
    vals = [x for x in vals if x is not None]
    return (stats.median(vals) if vals else None)

# 1) purity from VAF (heuristic)
vafs = []
for r in rows:
    m = r.get("metrics") or {}
    v = m.get("VAF")
    # prefer likely somatic: coding & not obviously common in population
    if r.get("coding") and (r.get("af") is None or r.get("af") < 0.001):
        if v is not None and 0.02 <= v <= 0.9:
            vafs.append(v)
purity_est = None
if vafs:
    # use the robust median; for noisy data you can switch to a small kernel density / histogram mode
    v_med = _safe_median(vafs)
    if v_med is not None:
        purity_est = max(0.0, min(1.0, 2.0 * v_med))

# 2) sv burden (only if your VCF contains SVTYPE)
sv_burden = 0
for r in rows:
    info = r.get("consequence","")  # placeholder; your VCF parsing could capture SVTYPE earlier if you start calling SVs
# If you later parse SVTYPE from INFO, count per-type here.

# 3) simple counts
snv = sum(1 for r in rows if len(r.get("ref",""))==1 and len(r.get("alt",""))==1)
ind = len(rows) - snv
# Ti/Tv (only on SNVs)
def _is_transition(ref,alt):
    return (ref,alt) in {("A","G"),("G","A"),("C","T"),("T","C")}
titv_num = sum(1 for r in rows if len(r["ref"])==1 and len(r["alt"])==1 and _is_transition(r["ref"],r["alt"]))
titv_den = max(1, snv - titv_num)
titv = (titv_num / titv_den) if snv else None

med_dp   = _safe_median([(r.get("metrics") or {}).get("DP") for r in rows])
med_qual = _safe_median([(r.get("metrics") or {}).get("QUAL") for r in rows])

# seed a summary (leave fields None if unknown)
summary = {
  "tumor_content": (round(purity_est*100,1) if purity_est is not None else None),
  "microbial_species": None,
  "mutation_signatures": [],
  "sv_burden": sv_burden,
  "counts": {"snv": snv, "indel": ind, "titv": titv, "median_dp": med_dp, "median_qual": med_qual}
}
clinical_comments = [make_comment(r) for r in rows
                     if (r.get('impact') or '').upper() in {'HIGH','MODERATE'} and r.get('coding')]
if clinical_comments:
    summary["clinical_comments"] = clinical_comments[:50]  # optional cap

# ---------- assemble report ----------
report={
  "meta":{
    "report_version":"0.0.1",
    "reference":"!{params.vep_assembly}",
    "pipeline":"OpenCare Clinical (panel-enabled)",
    "pipeline_version":"v0.6.2",
    "panel_name": "!{params.panel_name}" if "!{params.panel_name}" else None,
    "assay":"!{params.assay}",
    "filters":{"min_qual":MIN_QUAL,"min_dp":MIN_DP,"require_pass":REQUIRE_PASS,
               "apply_ab":APPLY_AB,"ab_min":AB_MIN,"ab_max":AB_MAX,"coding_only":CODING_ONLY},
    "disclaimer":"Not for clinical decision-making without local validation."
  },
  "patient":{"id":"!{params.patient_id}","report_date":None},"alt_id": sid,
  "qc":({"panel": qc_panel} if qc_panel else {}),
  "summary": summary,
  "rows":rows,
  "pharmcat": phc_obj
}

if qc_summary_block:
    report.setdefault("qc", {}).setdefault("summary", {}).update(qc_summary_block)
from collections import Counter

impact_counts = Counter((r.get('impact') or 'NA').upper() for r in rows)
type_counts   = Counter((r.get('type') or 'other') for r in rows)

top_variants = []
for r in rows[:min(50, len(rows))]:
    top_variants.append({
        "gene": r.get("gene"),
        "hgvsp": r.get("hgvsp"),
        "impact": r.get("impact"),
        "type": r.get("type"),
        "tier": r.get("tier"),
        "vaf": (r.get("metrics") or {}).get("VAF"),
        "dp":  (r.get("metrics") or {}).get("DP"),
        "qual":(r.get("metrics") or {}).get("QUAL")
    })

summary.update({
    "impact_counts": dict(impact_counts),
    "type_counts": dict(type_counts),
    "top_variants": top_variants
})

# --- Enrich rows from the web if enabled ---
# Uses CIViC always; adds OncoKB if ONCOKB_TOKEN is set.
if os.getenv("ENABLE_EVIDENCE", "0") == "1":
    enrich_rows_with_online(report)

# minimal mCODE-ish bundle
bundle={"resourceType":"Bundle","type":"collection","entry":[]}
bundle["entry"].append({"resource":{"resourceType":"Patient","id":"!{params.patient_id}"}})
for i,r in enumerate(rows,1):
    bundle["entry"].append({
      "resource":{
        "resourceType":"Observation","id":f"var-{i}","status":"final",
        "code":{"text":"Genetic variant"},
        "subject":{"reference":"Patient/!{params.patient_id}"},
        "component":[
          {"code":{"text":"Gene"},"valueString":r.get("gene","")},
          {"code":{"text":"HGVSp"},"valueString":r.get("hgvsp","")},
          {"code":{"text":"Consequence"},"valueString":r.get("consequence","")},
          {"code":{"text":"Location"},"valueString":f"{r['chrom']}:{r['pos']} {r['ref']}>{r['alt']}"}
        ]
      }
    })
# write a tab-delimited dump of all kept rows (optional)
if EXPORT_TSV and rows:
    import csv, gzip
    cols = ['chrom','pos','ref','alt','gene','consequence','hgvsp',
            'impact','qual','depth','ab','type','coding','tier']
    with gzip.open(f"{sid}.all_rows.tsv.gz", "wt", newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter=TAB)
        w.writeheader()
        for r in rows:
            out = {k: ('' if r.get(k) is None else r.get(k, '')) for k in cols}
            w.writerow(out)


Path(f"{sid}.report.json").write_text(json.dumps(report,indent=2))
Path(f"{sid}.mcode_bundle.json").write_text(json.dumps(bundle,indent=2))
PY
  '''
}

process FILTER_FOR_HTML {
  tag { sid }
  container "${IMG_BCFTOOLS}"
  cpus 6
  memory '20 GB'
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
TMP="!{sid}.filled.vcf.gz"
OUT="!{sid}.reportready.vcf.gz"

# Hook up the params you defined up top
USE_CSQ="!{ params.report_use_csq ? '1' : '0' }"
CSQ_REGEX='!{params.report_csq_regex}'

# 1) Fill only the tags we actually use (DP is not filled by +fill-tags)
bcftools +fill-tags "$IN" -Oz -o "$TMP" -- -t AN,AC,AF
tabix -f -p vcf "$TMP"

# 2) Liberal pass on genotype + basic depth/qual (no CSQ yet)
BASE_EXPR='(GT="0/1" || GT="1/1" || GT="1|0" || GT="0|1" || GT="1|1") &&
           QUAL>=20 &&
           (FMT/DP[0]>=8 || INFO/DP>=8 || (FMT/AD[0:0]+FMT/AD[0:1])>=8)'

bcftools view -f .,PASS -i "$BASE_EXPR" "$TMP" -Oz -o pass1.vcf.gz
tabix -f -p vcf pass1.vcf.gz
PASS1=$(bcftools view -H pass1.vcf.gz | wc -l || true)

# 3) Optional CSQ gate (quote-safe, outside bcftools)
if [ "$USE_CSQ" = "1" ]; then
  bcftools view -h pass1.vcf.gz > header.vcf
  bgzip -cd pass1.vcf.gz | grep -v '^#' | grep -E 'CSQ=.*('"$CSQ_REGEX"')' \
    | cat header.vcf - | bgzip -c > "$OUT" || :
else
  cp -f pass1.vcf.gz "$OUT"
fi

# 4) Always leave an index so Nextflow is happy
if ! tabix -f -p vcf "$OUT" 2>/dev/null; then : > "$OUT.tbi"; fi

# 5) Quick diagnostics
echo "n_before=${PASS1} n_after=$(bcftools view -H "$OUT" | wc -l || true)" >&2


'''
}



process BUILD_HTML {
  tag "$sample_id"
  container 'python:3.11'
  cpus 6
  memory '20 GB'
  time '30m'

  input:
  tuple val(sample_id), path(report_json)
  path  build_script
  path  pathway_db
  path  gene_domains

  output:
  path "OpenCare_${sample_id.toString().replaceAll(/[^A-Za-z0-9._-]/,'_')}_report.html"
  publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true

  script:
  def PD = pathway_db   && pathway_db.size()   > 0 ? "--pathway-db '${pathway_db}'"     : ''
  def GD = gene_domains && gene_domains.size() > 0 ? "--gene-domains '${gene_domains}'" : ''
  def safe_id = sample_id.toString().replaceAll(/[^A-Za-z0-9._-]/,'_')

  """
  set -euo pipefail
  mkdir -p "modules_${safe_id}" "out_${safe_id}"

 

  #  evidence enrichment (CIViC always; add OncoKB if token present)
  export ENABLE_EVIDENCE='!{ params.enable_evidence ? "1" : "0" }'
  # Only export ONCOKB_TOKEN if provided; leaving it empty = CIViC-only mode.
  if [ -n "!{ params.oncokb_token ?: '' }" ]; then
    export ONCOKB_TOKEN="!{ params.oncokb_token }"
  fi

  python3 "${build_script}" \
    --report "${report_json}" \
    --modules "modules_${safe_id}" \
    --outdir  "out_${safe_id}" \
    ${PD} ${GD}

  cp "out_${safe_id}/OpenCARE2_report.html" "OpenCare_${safe_id}_report.html"
  """
}
