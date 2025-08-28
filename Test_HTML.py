# Create a modular build of the OpenCARE PRO preview and generate a responsive HTML
# with your requested changes. We split the HTML into small "modules" that the
# builder stitches together. Chosen language: **Python** (reasoning below in chat).
#
# Outputs:
# - /mnt/data/opencare_modules/* (partials)
# - /mnt/data/opencare_build/opencare_pro_modular.html (final file to open)

import os, json, datetime, urllib.parse
import hashlib
from textwrap import dedent

base_dir = "/mnt/data"
mod_dir = os.path.join(base_dir, "opencare_modules")
out_dir = os.path.join(base_dir, "opencare_build")
os.makedirs(mod_dir, exist_ok=True)
os.makedirs(out_dir, exist_ok=True)
out_file = os.path.join(out_dir, "opencare_pro_modular.html")

# --- Assets / data -----------------------------------------------------------------
svg_logo = dedent("""
<svg xmlns="http://www.w3.org/2000/svg" width="180" height="40" viewBox="0 0 540 120">
  <defs>
    <linearGradient id="g" x1="0" x2="1" y1="0" y2="1">
      <stop offset="0" stop-color="#1aa7b6"/>
      <stop offset="1" stop-color="#0d6e76"/>
    </linearGradient>
  </defs>
  <rect fill="none" width="540" height="120"/>
  <g fill="url(#g)">
    <text x="10" y="85" font-family="Segoe UI, Arial, sans-serif" font-size="88" font-weight="700">OpenCARE</text>
  </g>
  <path d="M62 48c12-16 32-22 48-18-9 10-26 20-48 18z" fill="#12c3d6"/>
</svg>
""").strip()
logo_data_uri = "data:image/svg+xml;utf8," + urllib.parse.quote(svg_logo)

report = {
    "meta": {
        "brand": "OpenCARE",
        "report_version": "0.3.7-demo",
        "pipeline": "NF Clinical Toy",
        "pipeline_version": "v0.6.2",
        "reference": "GRCh38",
        "vep_cache": "VEP 110 (GRCh38)",
        "pharmcat": "PharmCAT 2.3 (toy)",
        "generated_at": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "disclaimer": "Demo only — not for clinical decision-making.",
        "params": {"ui": "ipr"},
        "containers": []
    },
    "patient": {
        "id": "PATIENT01",
        "alt_id": "ALT-2025-001",
        "gender": "Female",
        "age_group": "Adult",
        "report_date": "2025-08-27",
        "biopsy_name": "biop2",
        "biopsy_details": "Ultrasound-guided bx of Liver / VCC",
        "physician": "REDACTED",
        "cancer": "Cholangiocarcinoma"
    },
    "summary": {
        "tumor_content": "58%",
        "microbial_species": "None",
        "mutation_signatures": ["ID1 (DNA replication slippage)", "DBS5 (Platinum)", "SBS31 (Platinum; possible HRD)"],
        "sv_burden": "48 (33%)"
    },
    "qc": {
        "coverage_mean": 132,
        "coverage_median": 127,
        "on_target_rate": 0.91,
        "contamination": 0.011,
        "dup_rate": 0.09,
        "insert_median": 320,
        "read_q30": 0.93,
        "pct_gt20x": 0.96,
        "samples": [
            {"name": "SAMPLE1", "mean": 132, "median": 127, "pct20x": 0.96, "pct50x": 0.85, "pct100x": 0.71}
        ]
    },
    "rows": [
        {"chrom":"7","pos":55249071,"ref":"T","alt":"G","gene":"EGFR","consequence":"missense_variant","hgvsp":"p.L858R","cancer_type":"NSCLC","association":"Predictive—Benefit","tier":"AMP T1","metrics":{"DP":178,"AD":[130,48],"VAF":0.27,"QUAL":583},
         "links":{"ClinVar":"https://www.ncbi.nlm.nih.gov/clinvar/?term=EGFR+L858R","OncoKB":"https://www.oncokb.org/gene/EGFR","IGV":"#"},
         "treatments":[
            {"therapy":"Osimertinib","evidence":"NCCN guideline","source_url":"https://www.nccn.org/guidelines/guidelines-detail?category=1&id=1450","drug_label":"https://www.accessdata.fda.gov/drugsatfda_docs/label/2020/208065s013lbl.pdf","pmid":"29151359"},
            {"therapy":"Gefitinib","evidence":"FDA label","source_url":"https://www.accessdata.fda.gov/drugsatfda_docs/label/2015/206995Orig1s000lbl.pdf","drug_label":"https://www.accessdata.fda.gov/drugsatfda_docs/label/2015/206995Orig1s000lbl.pdf","pmid":"12629580"}
         ]},
        {"chrom":"12","pos":25398285,"ref":"G","alt":"T","gene":"KRAS","consequence":"missense_variant","hgvsp":"p.G12C","cancer_type":"NSCLC","association":"Predictive—Benefit","tier":"AMP T1","metrics":{"DP":142,"AD":[86,56],"VAF":0.39,"QUAL":421},
         "links":{"ClinVar":"https://www.ncbi.nlm.nih.gov/clinvar/?term=KRAS+G12C","OncoKB":"https://www.oncokb.org/gene/KRAS","IGV":"#"},
         "treatments":[{"therapy":"Sotorasib","evidence":"FDA label","source_url":"https://www.accessdata.fda.gov/drugsatfda_docs/label/2023/214665s005lbl.pdf","drug_label":"https://www.accessdata.fda.gov/drugsatfda_docs/label/2023/214665s005lbl.pdf","pmid":"34133884"}]},
        {"chrom":"12","pos":25398284,"ref":"G","alt":"A","gene":"KRAS","consequence":"missense_variant","hgvsp":"p.G12D","cancer_type":"CRC","association":"Predictive—Benefit","tier":"AMP T2","metrics":{"DP":121,"AD":[77,44],"VAF":0.36,"QUAL":337},
         "links":{"ClinVar":"https://www.ncbi.nlm.nih.gov/clinvar/?term=KRAS+G12D","OncoKB":"https://www.oncokb.org/gene/KRAS","IGV":"#"},
         "treatments":[{"therapy":"Adagrasib (investigational)","evidence":"Phase II","source_url":"https://clinicaltrials.gov/search?cond=Colorectal+Cancer&term=adagrasib+KRAS+G12D","drug_label":"https://www.fda.gov/","pmid":"38338988"}]},
        {"chrom":"17","pos":430456,"ref":"C","alt":"A","gene":"BRCA1","consequence":"stop_gained","hgvsp":"p.E23*","cancer_type":"Breast/Ovarian","association":"Predictive—Benefit","tier":"AMP T2","metrics":{"DP":98,"AD":[60,38],"VAF":0.39,"QUAL":290},
         "links":{"ClinVar":"https://www.ncbi.nlm.nih.gov/clinvar/?term=BRCA1+loss+of+function","OncoKB":"https://www.oncokb.org/gene/BRCA1","IGV":"#"},
         "treatments":[{"therapy":"Olaparib","evidence":"ESMO guideline","source_url":"https://www.esmo.org/guidelines","drug_label":"https://www.ema.europa.eu/en/documents/product-information/lynparza-epar-product-information_en.pdf","pmid":"37285449"}]},
        {"chrom":"10","pos":89692904,"ref":"G","alt":"A","gene":"PTEN","consequence":"missense_variant","hgvsp":"p.R130Q","cancer_type":"Solid tumors","association":"Predictive—Resistance","tier":"Biological (non-actionable)","metrics":{"DP":110,"AD":[95,15],"VAF":0.14,"QUAL":210},
         "links":{"ClinVar":"https://www.ncbi.nlm.nih.gov/clinvar/?term=PTEN+R130Q","OncoKB":"https://www.oncokb.org/gene/PTEN","IGV":"#"},
         "treatments":[{"therapy":"AKT inhibitor (class)","evidence":"Phase II","source_url":"https://clinicaltrials.gov/search?term=AKT+inhibitor+PTEN","drug_label":"https://clinicaltrials.gov/search?term=AKT+inhibitor+PTEN","pmid":"34819606"}]},
        {"chrom":"2","pos":29443694,"ref":"N","alt":"EML4-ALK","gene":"ALK","consequence":"translocation","hgvsp":"EML4-ALK","cancer_type":"Cholangiocarcinoma","association":"Predictive—Benefit","tier":"AMP T1","metrics":{"DP":0,"AD":[0,0],"VAF":None,"QUAL":None},
         "links":{"ClinVar":"https://www.ncbi.nlm.nih.gov/clinvar/?term=ALK+EML4","OncoKB":"https://www.oncokb.org/gene/ALK","IGV":"#"},
         "treatments":[
            {"therapy":"Alectinib","evidence":"NCCN guideline","source_url":"https://www.nccn.org/guidelines/guidelines-detail?category=1&id=1450","drug_label":"https://www.accessdata.fda.gov/drugsatfda_docs/label/2022/208434s018lbl.pdf","pmid":"28668676"},
            {"therapy":"Brigatinib","evidence":"Guideline-supported","source_url":"https://www.nccn.org/guidelines/guidelines-detail?category=1&id=1450","drug_label":"https://www.accessdata.fda.gov/drugsatfda_docs/label/2024/208772s014lbl.pdf","pmid":"28475456"}
         ]}
    ]
}

gene_domains = {
    "EGFR":{"length":1210,"domains":[["LBD",25,645],["TM",646,668],["Kinase",712,978]]},
    "KRAS":{"length":189,"domains":[["G-domain",1,166],["HVR",167,188]]},
    "BRCA1":{"length":1863,"domains":[["RING",24,64],["BRCT",1646,1859]]},
    "ALK":{"length":1620,"domains":[["Kinase",1116,1378]]},
    "PTEN":{"length":403,"domains":[["Phosphatase",14,185],["C2",190,351]]}
}

# --- Partials ----------------------------------------------------------------------
head_css = dedent(r"""
<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>OpenCARE · Pipeline Reports — PRO v3.2</title>
<style>
:root{--brand:#146c74;--brand-strong:#0e4e54;--accent:#0d8ea0;--bg:#ffffff;--fg:#0f172a;--muted:#475569;--card:#fff;--border:#e5e7eb;--th:#f3f4f6;--hover:#f9fafb;--zebra:#fafafa;--ok:#0fb98a;--warn:#f59e0b;--danger:#ef4444;--chip:#e6f7fb;--chipfg:#064e57;--shadow:0 10px 24px rgba(0,0,0,.06)}
[data-theme="dark"]{--brand:#0b3f44;--brand-strong:#072a2e;--accent:#46c2cf;--bg:#081217;--fg:#e5f6f8;--muted:#97c2c8;--card:#0c1b20;--border:#13353a;--th:#0f2429;--hover:#0c1f24;--zebra:#0a191e;--chip:#0e2b31;--chipfg:#c8f0f4;--shadow:0 12px 28px rgba(0,0,0,.4)}
*{box-sizing:border-box}body{margin:0;font-family:system-ui,Segoe UI,Roboto,Arial;background:var(--bg);color:var(--fg)}a{color:var(--accent)}:focus-visible{outline:3px solid color-mix(in oklab,var(--accent) 60%,transparent);outline-offset:2px;border-radius:6px}
.appbar{position:sticky;top:0;z-index:60;background:linear-gradient(90deg,var(--brand),var(--brand-strong));color:#fff;padding:8px 14px;display:flex;gap:16px;align-items:center;justify-content:space-between}
.brand{display:flex;align-items:center;gap:10px}.brand img{height:28px;display:block;filter:drop-shadow(0 1px 0 rgba(255,255,255,.7)) contrast(1.05)}
[data-theme="dark"] .brand img{filter:drop-shadow(0 1px 0 rgba(0,0,0,.9)) brightness(1.25) contrast(1.2)}
.brand .name{font-weight:900;letter-spacing:.2px;font-size:14px}
.appbar .seg{display:flex;border:1px solid rgba(255,255,255,.35);border-radius:12px;overflow:hidden}
.appbar .seg button{border:0;background:transparent;padding:8px 12px;color:#fff;opacity:.9;cursor:pointer}
.appbar .seg button.active{background:rgba(255,255,255,.18);opacity:1}
.appbar .btn{background:transparent;border:1px solid rgba(255,255,255,.35);color:#fff;border-radius:10px;padding:8px 12px;cursor:pointer}
.container{display:block;padding:12px 2vw;max-width:min(2000px,96vw)} /* responsive wide */
.card{background:var(--card);border:1px solid var(--border);border-radius:14px;box-shadow:var(--shadow);padding:14px;margin-right:78px}
.card h2{margin:0 0 8px}.section-title{display:flex;align-items:center;justify-content:space-between}.section-title small{color:var(--muted)}
.kv{display:grid;grid-template-columns:minmax(220px, 22vw) 1fr;gap:8px 14px}.kv .label{color:var(--muted)}.sep{height:1px;background:var(--border);margin:12px 0}
.sidebar{position:fixed;right:8px;top:64px;z-index:70}
.sideicons{display:flex;flex-direction:column;gap:10px}
.sideicons a{position:relative;display:flex;align-items:center;justify-content:center;width:44px;height:44px;border:1px solid var(--border);border-radius:12px;background:var(--card);text-decoration:none;color:var(--fg);transition:transform .18s ease, background-color .18s ease, color .18s ease}
.sideicons a:hover{transform:translateX(-8px);background:var(--accent);color:#fff}
.sideicons a.active{outline:2px solid color-mix(in oklab,var(--accent) 45%,transparent)}
.sideicons a span{font-size:20px}
.sideicons a .label{position:absolute;right:52px;opacity:0;pointer-events:none;background:var(--card);color:var(--fg);border:1px solid var(--border);border-radius:8px;padding:4px 8px;font-size:12px;transform:translateX(8px);transition:opacity .18s ease, transform .18s ease}
.sideicons a:hover .label{opacity:1;transform:translateX(0)}
.controls{display:flex;flex-wrap:wrap;gap:8px;align-items:center;margin:0 0 8px}
input[type="search"]{flex:1;min-width:280px;padding:10px;border:1px solid var(--border);border-radius:10px;background:transparent;color:var(--fg)}
.tablewrap{max-height:65vh;overflow:auto;border:1px solid var(--border);border-radius:12px;background:var(--card);box-shadow:var(--shadow)}
table{border-collapse:collapse;width:100%} /* expand with screen */
th,td{border-bottom:1px solid var(--border);padding:10px 12px;vertical-align:top}thead th{position:sticky;top:0;background:var(--th);z-index:2;cursor:pointer;text-align:left}tbody tr:hover{background:var(--hover)}tbody tr:nth-child(even){background:var(--zebra)}
td:first-child, th:first-child{position:sticky;left:0;z-index:3;background:var(--card)}
.badge{display:inline-block;background:var(--chip);color:var(--chipfg);border-radius:999px;padding:4px 10px;margin:2px 4px 2px 0;font-size:12px}
small.muted{color:var(--muted)}#counts{border:1px solid var(--border);padding:4px 10px;border-radius:999px}
th .arrow{opacity:.4;margin-left:4px}th.sorted .arrow{opacity:1}.alert{color:var(--warn);font-weight:700;margin-right:6px}
.linkify a{font-size:12px;text-decoration:none;border-bottom:1px dotted var(--accent);color:var(--accent)}.linkify a:hover{text-decoration:underline}
.geneBtn{background:transparent;border:0;color:var(--accent);text-decoration:underline;cursor:pointer;padding:0 4px;border-radius:6px}
.inspectIcon{margin-left:6px;cursor:pointer;border:1px solid var(--border);background:var(--card);border-radius:8px;padding:2px 6px}
.anchorlink{opacity:.6;text-decoration:none;margin-left:6px}
.tier{display:inline-flex;align-items:center;gap:6px;border:1px solid var(--border);border-radius:999px;padding:2px 8px}.tier .dot{width:8px;height:8px;border-radius:50%}
.evwrap{display:flex;flex-wrap:nowrap;gap:6px;overflow:auto}
.evgroup{flex:0 0 auto;border:1px dashed var(--border);border-radius:10px;padding:4px 8px;background:var(--zebra)}
.evlabel{font-size:11px;color:var(--muted);display:block;margin-bottom:3px}
.evchips a{margin-right:6px;white-space:nowrap}
.mini{display:inline-flex;gap:6px;flex-wrap:wrap}
.mini span{border:1px solid var(--border);border-radius:8px;padding:2px 6px;font-size:12px;background:var(--zebra)}
#pathwayChips, #qcChips{display:flex;flex-wrap:wrap;gap:8px;margin-bottom:8px}
#pathwayChips .chip, #qcChips .chip{border:1px solid var(--border);border-radius:999px;padding:6px 10px;background:var(--card);cursor:pointer}
#pathwayChips .chip.active, #qcChips .chip.active{background:var(--accent);color:#fff;border-color:transparent}
.pcard{border:1px solid var(--border);border-radius:12px;padding:12px;background:var(--card)}
.pflex{display:grid;grid-template-columns:2fr 1fr;gap:12px} /* widen left on big screens */
.node{cursor:pointer;transition:transform .08s ease}.node:hover{transform:scale(1.04)}
#mechPanel{border:1px dashed var(--border);border-radius:12px;padding:10px}
#inspector{position:fixed;right:0;top:0;bottom:0;width:520px;max-width:96vw;background:var(--card);border-left:1px solid var(--border);box-shadow:var(--shadow);transform:translateX(100%);transition:transform .25s ease;z-index:80}#inspector.open{transform:translateX(0)}.ins-header{display:flex;justify-content:space-between;align-items:center;padding:12px 14px;border-bottom:1px solid var(--border)}.ins-body{padding:12px 14px;overflow:auto;height:calc(100% - 52px)}.card-lite{border:1px dashed var(--border);border-radius:10px;padding:10px;margin-top:6px}.kv2{display:grid;grid-template-columns:140px 1fr;gap:8px}
#clinCards{display:grid;grid-template-columns:repeat(auto-fill,minmax(280px,1fr));gap:10px}
.cCard{border:1px solid var(--border);border-radius:12px;padding:10px;background:var(--card)}.cHead{display:flex;justify-content:space-between;align-items:center;margin-bottom:6px}.cGene{font-weight:700}.cType{background:var(--chip);color:var(--chipfg);border-radius:999px;padding:2px 8px;font-size:12px}
#help{position:fixed;inset:0;background:rgba(0,0,0,.55);display:none;align-items:center;justify-content:center;z-index:90}#help .box{background:var(--card);color:var(--fg);border:1px solid var(--border);border-radius:14px;box-shadow:var(--shadow);padding:16px;max-width:700px;width:96vw}kbd{border:1px solid var(--border);padding:2px 6px;border-radius:6px}
@media (min-width:1600px){
  .container{max-width:min(2200px,96vw)}
  .kv{grid-template-columns:minmax(260px, 26vw) 1fr}
  th,td{padding:12px 16px}
}
@media print{.appbar,.sidebar,.controls,.pagination,#inspector,#help,#geneViewer{display:none!important}.container{padding:0}.card{page-break-inside:avoid;box-shadow:none;border-color:#888}a[href]:after{content:" (" attr(href) ")";font-size:10px;color:#666}}
.anchor{scroll-margin-top:76px}
</style>
</head>
<body>
""")

appbar = dedent("""
<div class="appbar" role="banner">
  <div class="brand">
    <img src="LOGO_URI" alt="OpenCARE logo" />
    <div class="name">OpenCARE · <span style="opacity:.9">Pipeline Reports</span></div>
    <span class="pill" id="metaVer" aria-label="Report version" style="border:1px solid rgba(255,255,255,.25);border-radius:999px;padding:3px 8px;font-size:12px;opacity:.9"></span>
    <div class="seg" role="tablist" aria-label="Theme" style="margin-left:8px">
      <button id="lightBtn" role="tab" aria-selected="true" class="active">Light</button>
      <button id="darkBtn" role="tab" aria-selected="false">Dark</button>
    </div>
    <div class="seg" role="tablist" aria-label="View mode">
      <button id="summaryBtn" role="tab" aria-selected="true" class="active">Summary</button>
      <button id="detailedBtn" role="tab" aria-selected="false">Detailed</button>
    </div>
  </div>
  <div class="right">
    <button class="btn" id="helpBtn" aria-haspopup="dialog" aria-controls="help">?</button>
    <button class="btn" onclick="window.print()" aria-label="Print to PDF">Print/PDF</button>
    <button class="btn" id="csvBtn">CSV</button>
    <button class="btn" id="jsonBtn">JSON</button>
    <button class="btn" id="mcodeBtn">FHIR/mCODE</button>
  </div>
</div>
""")

summary = dedent("""
<div class="container">
  <main>
    <section id="summary" class="card anchor" aria-labelledby="h-summary">
      <div class="section-title">
        <h2 id="h-summary">1. Summary</h2>
        <small id="crumbs"></small>
      </div>

      <div class="kv" role="list">
        <div class="label" role="listitem">Diagnosis</div><div id="diag"></div>
        <div class="label" role="listitem">Patient ID</div><div id="pid"></div>
        <div class="label" role="listitem">Alternate ID</div><div id="altid"></div>
        <div class="label" role="listitem">Report Date</div><div id="rdate"></div>
        <div class="label" role="listitem">Case Type</div><div id="ctype"></div>
        <div class="label" role="listitem">Physician</div><div id="phys"></div>
        <div class="label" role="listitem">Biopsy Name</div><div id="bname"></div>
        <div class="label" role="listitem">Biopsy Details</div><div id="bdet"></div>
        <div class="label" role="listitem">Gender</div><div id="gender"></div>
      </div>
      <div class="sep"></div>

      <div id="clinCommentsSec">
        <h3>Clinical Comments</h3>
        <div id="clinCards"></div>
        <div class="sep"></div>
      </div>

      <h3>Tumour Summary</h3>
      <div class="kv">
        <div class="label">Tumour Content</div><div id="tcontent"></div>
        <div class="label">Microbial Species</div><div id="microb"></div>
        <div class="label">Mutation Signatures</div><div id="mutsig"></div>
        <div class="label">SV Burden</div><div id="svb"></div>
      </div>

      <div id="qcBlock" style="display:none">
        <div class="sep"></div>
        <h3>QC Summary</h3>
        <div class="kv" id="qcSummary"></div>
        <div style="height:6px"></div>
        <div id="qcChips"></div>
        <div id="qcPanel" class="pcard"></div>
      </div>
    </section>
""")

therapeutic = dedent("""
    <section id="therapeutic" class="card anchor" aria-labelledby="h-ther">
      <div class="section-title">
        <h2 id="h-ther">2. Therapeutic Alterations (Detailed)</h2>
        <small id="counts" aria-live="polite"></small>
      </div>
      <div class="controls">
        <input id="q" type="search" placeholder="Filter gene / variant / organ / therapy / evidence…" aria-label="Quick filter"/>
        <button class="btn" id="resetBtn">Reset</button>
      </div>
      <div class="tablewrap" role="region" aria-label="Therapeutic Alterations Table">
        <table id="tbl" role="table" aria-describedby="h-ther">
          <thead>
            <tr role="row">
              <th role="columnheader" data-k="gene">Gene <span class="arrow">↕</span></th>
              <th role="columnheader" data-k="hgvsp">Observed Variant <span class="arrow">↕</span></th>
              <th role="columnheader" data-k="known">Known/Hotspot <span class="arrow">↕</span></th>
              <th role="columnheader" data-k="locus">Locus <span class="arrow">↕</span></th>
              <th role="columnheader" data-k="cancer_type">Cancer Type <span class="arrow">↕</span></th>
              <th role="columnheader" data-k="association">Clinical Association <span class="arrow">↕</span></th>
              <th role="columnheader" data-k="tier">Tier <span class="arrow">↕</span></th>
              <th role="columnheader">Metrics</th>
              <th role="columnheader">Treatment Evidence</th>
            </tr>
          </thead>
          <tbody></tbody>
        </table>
      </div>
    </section>
""")

pathways = dedent("""
    <section id="pathways" class="card anchor" aria-labelledby="h-path">
      <div class="section-title">
        <h2 id="h-path">3. Pathway Analysis</h2>
        <small>Pick an affected pathway, then click a highlighted node for ~50-word mechanism + references</small>
      </div>
      <div id="pathwayChips"></div>
      <div id="pathwayOne" class="pcard pflex"></div>
    </section>
  </main>
</div>
""")

sidebar = dedent("""
<aside class="sidebar" role="navigation" aria-label="Report Sections">
  <div class="sideicons">
    <a href="#summary"><span>🧾</span><em class="label">Summary</em></a>
    <a href="#therapeutic"><span>💊</span><em class="label">Therapeutic</em></a>
    <a href="#pathways"><span>🧬</span><em class="label">Pathways</em></a>
  </div>
</aside>
""")

inspector = dedent("""
<div id="inspector" aria-hidden="true" role="dialog" aria-labelledby="insTitle">
  <div class="ins-header">
    <strong id="insTitle">Variant</strong>
    <button class="btn" id="insClose" aria-label="Close inspector">Close</button>
  </div>
  <div class="ins-body">
    <div class="kv2">
      <div class="label">Gene</div><div id="insGene"></div>
      <div class="label">Protein</div><div id="insProt"></div>
      <div class="label">Consequence</div><div id="insCons"></div>
      <div class="label">Tier</div><div id="insTier"></div>
      <div class="label">Locus</div><div id="insLocus"></div>
      <div class="label">Cancer Type</div><div id="insCancer"></div>
      <div class="label">Association</div><div id="insAssoc"></div>
      <div class="label">Depth (DP)</div><div id="insDP"></div>
      <div class="label">Alleles (AD)</div><div id="insAD"></div>
      <div class="label">VAF</div><div id="insVAF"></div>
      <div class="label">QUAL</div><div id="insQUAL"></div>
    </div>
    <div class="card-lite" id="insLinks" aria-label="External links"></div>
    <h4>Treatments</h4>
    <div id="insTreats"></div>
  </div>
</div>
""")

help_modal = dedent("""
<div id="help" role="dialog" aria-modal="true" aria-labelledby="helpTitle" style="display:none">
  <div class="box">
    <h3 id="helpTitle">Keyboard & Tips</h3>
    <ul><li><kbd>/</kbd> focus search</li><li><kbd>Esc</kbd> close modals</li></ul>
    <div style="text-align:right"><button class="btn" id="helpClose">Close</button></div>
  </div>
</div>
""")

scripts = dedent(r"""
<script>
const report = REPORT_JSON;
const geneDomains = GENE_DOMAINS;

/* Top meta */
function set(id, v){ const el=document.getElementById(id); if(el) el.textContent = v || ''; }
document.getElementById('metaVer').textContent = `${report.meta.report_version} · ${report.meta.reference}`;
set('diag', report.patient.cancer); set('pid', report.patient.id); set('altid', report.patient.alt_id);
set('rdate', report.patient.report_date); set('ctype', report.patient.age_group); set('phys', report.patient.physician);
set('bname', report.patient.biopsy_name); set('bdet', report.patient.biopsy_details); set('gender', report.patient.gender);
document.getElementById('crumbs').textContent = `${report.patient.cancer} · ${report.patient.id}`;
const sum = report.summary || {}; set('tcontent', sum.tumor_content || '—'); set('microb', sum.microbial_species || '—'); set('mutsig', (sum.mutation_signatures||[]).join('; ')); set('svb', sum.sv_burden || '—');

/* Clinician cards (appear only in Summary view) */
const ccWrap = document.getElementById('clinCards');
ccWrap.innerHTML = (report.rows||[]).map(r => {
  const tx = (r.treatments||[]).map(t=>t.therapy).filter(Boolean).join(', ') || '—';
  return `<div class="cCard">
    <div class="cHead"><span class="cGene">${r.gene} ${r.hgvsp || r.consequence}</span>
      <span class="cType">${r.cancer_type || '—'}</span></div>
    <div style="display:flex;gap:8px;align-items:center;flex-wrap:wrap">
      <span class="badge">${r.tier}</span><span>${r.association}</span></div>
    <div style="margin-top:6px"><small class="muted">Suggested:</small> ${tx}</div>
  </div>`;
}).join('');

/* Theme + View (Summary / Detailed) */
const lightBtn = document.getElementById('lightBtn');
const darkBtn = document.getElementById('darkBtn');
const summaryBtn = document.getElementById('summaryBtn');
const detailedBtn = document.getElementById('detailedBtn');
function setTheme(t){ document.body.setAttribute('data-theme', t); localStorage.setItem('theme', t); }
function setView(v){
  const detailed = (v==='detailed');
  document.getElementById('therapeutic').style.display = detailed ? '' : 'none';
  document.getElementById('qcBlock').style.display = detailed ? '' : 'none';
  document.getElementById('clinCommentsSec').style.display = detailed ? 'none' : '';
  localStorage.setItem('viewmode', v);
}
const savedTheme = localStorage.getItem('theme'); if(savedTheme) setTheme(savedTheme);
let savedView = localStorage.getItem('viewmode') || 'summary';
if(savedView==='clinician') savedView='summary';
if(savedView==='analyst') savedView='detailed';
setView(savedView);
if(savedTheme==='dark'){ lightBtn.classList.remove('active'); darkBtn.classList.add('active'); }
if(savedView==='detailed'){ summaryBtn.classList.remove('active'); detailedBtn.classList.add('active'); }
lightBtn.onclick = ()=>{ setTheme('light'); lightBtn.classList.add('active'); darkBtn.classList.remove('active'); };
darkBtn.onclick = ()=>{ setTheme('dark'); darkBtn.classList.add('active'); lightBtn.classList.remove('active'); };
summaryBtn.onclick = ()=>{ setView('summary'); summaryBtn.classList.add('active'); detailedBtn.classList.remove('active'); };
detailedBtn.onclick = ()=>{ setView('detailed'); detailedBtn.classList.add('active'); summaryBtn.classList.remove('active'); };

/* Scroll spy */
const navLinks = Array.from(document.querySelectorAll('.sideicons a')); const sections = navLinks.map(a => document.querySelector(a.getAttribute('href')));
function onScrollSpy(){ let idx=0, y=window.scrollY+100; sections.forEach((sec,i)=>{ if(sec && sec.offsetTop < y) idx=i; }); navLinks.forEach((a,i)=> a.classList.toggle('active', i===idx)); }
window.addEventListener('scroll', onScrollSpy); onScrollSpy();

/* Helpers */
function escapeHtml(s){return (s||'').toString().replace(/[&<>]/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;'}[c]))}
function locusOf(r){ return `${r.chrom}:${r.pos} ${r.ref}>${r.alt}`; }
function assocIcon(a){ if(/Benefit/i.test(a)) return `<span style="color:var(--ok)">●</span> ${escapeHtml(a)}`; if(/Resistance/i.test(a)) return `<span style="color:var(--danger)">●</span> ${escapeHtml(a)}`; return escapeHtml(a||'—'); }
function tierBadge(tier){ const color = /T1/.test(tier)?'var(--ok)':(/T2/.test(tier)?'var(--warn)':'var(--muted)'); return `<span class="tier"><span class="dot" style="background:${color}"></span>${escapeHtml(tier||'—')}</span>`; }
function isLowQuality(m){ if(!m) return false; const dp = +m.DP || 0; const qual = +m.QUAL || 0; return (dp < 30) || (qual && qual < 200); }
const knownMap = {"EGFR:p.L858R":"Hotspot","KRAS:p.G12C":"Hotspot","KRAS:p.G12D":"Known","ALK:EML4-ALK":"Fusion","BRCA1:LoF":"LoF","PTEN:p.R130Q":"Known"};
function knownBadge(r){ const key1 = `${r.gene}:${r.hgvsp||''}`; const key2 = (/stop_gained|frameshift/i.test(r.consequence) && (r.gene==='BRCA1'||r.gene==='BRCA2')) ? `${r.gene}:LoF` : null; const label = knownMap[key1] || (key2 ? knownMap[key2] : null); return label ? `<span class="badge" title="Catalog match">${label}</span>` : '<small class="muted">—</small>'; }
function metricsMini(m){ if(!m) return '<small class="muted">—</small>'; const chips=[]; if(m.DP!=null) chips.push(`<span>DP ${m.DP}</span>`); if(m.VAF!=null) chips.push(`<span>VAF ${Math.round(m.VAF*1000)/10}%</span>`); if(m.QUAL!=null) chips.push(`<span>QUAL ${m.QUAL}</span>`); return `<div class="mini">${chips.join('')}</div>`; }

/* Therapeutic table */
const tbody = document.querySelector('#tbl tbody'); const counts = document.getElementById('counts'); const qInput = document.getElementById('q'); let state = { q:'', sortK:'tier', sortDir:-1 };
function applyFilters(){
  const q=state.q.toLowerCase(), f=v=>(v||'').toLowerCase().includes(q);
  let rows = report.rows.filter(r => (!q || f(r.gene)||f(r.hgvsp)||f(r.cancer_type)||f(r.association)||f(r.tier)||f(locusOf(r))||(r.treatments||[]).some(t=>f(t.therapy)||f(t.evidence))));
  if(state.sortK){
    rows.sort((a,b)=>{ const av=(state.sortK==='locus')? `${a.chrom}:${a.pos}` : (state.sortK==='known'? (knownMap[`${a.gene}:${a.hgvsp||''}`]||'') : (a[state.sortK]||'')); const bv=(state.sortK==='locus')? `${b.chrom}:${b.pos}` : (state.sortK==='known'? (knownMap[`${b.gene}:${b.hgvsp||''}`]||'') : (b[state.sortK]||'')); return av.localeCompare(bv,undefined,{numeric:true})*state.sortDir; });
  } return rows;
}
function evidenceBlock(r){
  const groups = {Guideline:[], Label:[], PubMed:[], Trial:[]};
  (r.treatments||[]).forEach(t=>{
    const therapy = escapeHtml(t.therapy||'');
    if(/NCCN|ESMO|Guideline/i.test(t.evidence||'')){ groups.Guideline.push(`<a href="${t.source_url||'#'}" target="_blank" rel="noopener" title="${escapeHtml(t.evidence||'Guideline')}">${therapy}</a>`); }
    if(/FDA|EMA|label/i.test((t.evidence||'') + (t.drug_label||''))){ groups.Label.push(`<a href="${t.drug_label||t.source_url||'#'}" target="_blank" rel="noopener" title="Drug label">${therapy}</a>`); }
    if(t.pmid){ const pm = (t.pmid||'').replace(/[^\\d]/g,''); groups.PubMed.push(`<a href="https://pubmed.ncbi.nlm.nih.gov/${pm}/" target="_blank" rel="noopener" title="PubMed ${escapeHtml(t.pmid)}">${therapy}</a>`); }
    if(/Phase|trial/i.test(t.evidence||'')){ groups.Trial.push(`<a href="${t.source_url||'#'}" target="_blank" rel="noopener" title="${escapeHtml(t.evidence)}">${therapy}</a>`); }
  });
  return `<div class="evwrap">`+
    (groups.Guideline.length? `<div class="evgroup"><span class="evlabel">Guideline</span><span class="evchips">${groups.Guideline.join('')}</span></div>`:'')+
    (groups.Label.length? `<div class="evgroup"><span class="evlabel">Label</span><span class="evchips">${groups.Label.join('')}</span></div>`:'')+
    (groups.PubMed.length? `<div class="evgroup"><span class="evlabel">PubMed</span><span class="evchips">${groups.PubMed.join('')}</span></div>`:'')+
    (groups.Trial.length? `<div class="evgroup"><span class="evlabel">Trial</span><span class="evchips">${groups.Trial.join('')}</span></div>`:'')+
  `</div>`;
}
function renderTable(){
  const all = applyFilters(); counts.textContent = `${all.length} variants`; tbody.innerHTML='';
  all.forEach((r,i)=>{
    const id = `var-${i+1}`; const alert = isLowQuality(r.metrics) ? `<span class="alert" title="Low-quality metrics (toy thresholds)">!</span>` : '';
    const tr=document.createElement('tr'); tr.id=id;
    tr.innerHTML = `
      <td><button class="geneBtn" data-gene="${escapeHtml(r.gene)}" data-hgvsp="${escapeHtml(r.hgvsp||'')}" title="Open gene map">${alert}${escapeHtml(r.gene||'')}</button>
          <button class="inspectIcon" data-idx="${i}" title="Inspect variant">🔎</button>
          <a href="#${id}" class="anchorlink" title="Deep link">#</a></td>
      <td>${escapeHtml(r.hgvsp||'')}</td>
      <td>${knownBadge(r)}</td>
      <td>${escapeHtml(locusOf(r))}</td>
      <td>${escapeHtml(r.cancer_type||'')}</td>
      <td>${assocIcon(r.association||'')}</td>
      <td>${tierBadge(r.tier||'')}</td>
      <td>${metricsMini(r.metrics)}</td>
      <td>${evidenceBlock(r)}</td>`;
    tbody.appendChild(tr);
  });
  tbody.querySelectorAll('.inspectIcon').forEach(btn=>{ btn.addEventListener('click',(e)=> openInspector(applyFilters()[parseInt(e.currentTarget.dataset.idx,10)]) ); });
  tbody.querySelectorAll('.geneBtn').forEach(btn=>{ btn.addEventListener('click', ()=> openGeneViewer(btn.dataset.gene, btn.dataset.hgvsp) ); });
}
document.getElementById('resetBtn').onclick = ()=>{ state.q=''; qInput.value=''; renderTable(); };
qInput.addEventListener('input', ()=>{ state.q=qInput.value.trim(); renderTable(); });
for(const th of document.querySelectorAll('thead th[data-k]')){ th.addEventListener('click', ()=>{ const k = th.dataset.k; if(state.sortK===k){ state.sortDir*=-1; } else { state.sortK=k; state.sortDir=1; } document.querySelectorAll('thead th').forEach(x=>x.classList.remove('sorted')); th.classList.add('sorted'); th.querySelector('.arrow').textContent = state.sortDir===1?'↑':'↓'; renderTable(); }); }
renderTable();

/* Inspector */
const insp = document.getElementById('inspector'); document.getElementById('insClose').onclick = ()=>{ insp.classList.remove('open'); insp.setAttribute('aria-hidden','true'); };
function openInspector(r){
  document.getElementById('insTitle').textContent = `${r.gene} · ${r.hgvsp || r.consequence}`;
  document.getElementById('insGene').textContent = r.gene||'';
  document.getElementById('insProt').textContent = r.hgvsp||'';
  document.getElementById('insCons').textContent = r.consequence||'';
  document.getElementById('insTier').innerHTML = `${tierBadge(r.tier||'')}`;
  document.getElementById('insLocus').textContent = `${r.chrom}:${r.pos} ${r.ref}>${r.alt}`;
  document.getElementById('insCancer').textContent = r.cancer_type||'';
  document.getElementById('insAssoc').innerHTML = assocIcon(r.association||'');
  const m = r.metrics||{};
  document.getElementById('insDP').textContent = (m.DP==null?'—':m.DP);
  document.getElementById('insAD').textContent = (m.AD? m.AD.join(', '):'—');
  document.getElementById('insVAF').textContent = (m.VAF==null?'—': (Math.round(m.VAF*1000)/10)+'%');
  document.getElementById('insQUAL').textContent = (m.QUAL==null?'—':m.QUAL);
  const L = r.links||{}; document.getElementById('insLinks').innerHTML = ['ClinVar','OncoKB','IGV'].map(k=> L[k] ? `<a href="${L[k]}" target="_blank" rel="noopener">${k}</a>` : '').filter(Boolean).join(' · ') || '<span class="muted">No links</span>';
  document.getElementById('insTreats').innerHTML = evidenceBlock(r);
  insp.classList.add('open'); insp.setAttribute('aria-hidden','false');
}

/* Protein domain viewer */
function openGeneViewer(gene, hgvsp){
  const info = geneDomains[gene]; if(!info){ alert('No domain map for '+gene); return; }
  const m = /[pP]\.[A-Za-z*]?(\d+)/.exec(hgvsp||''); const aa = m ? parseInt(m[1],10) : null;
  const W = 680, H = 120, pad = 20, len = info.length;
  const scale = (x)=> pad + x/len*(W-2*pad);
  let svg = `<svg viewBox="0 0 ${W} ${H}" width="100%" height="140">`;
  svg += `<rect x="${pad}" y="${50}" width="${W-2*pad}" height="${12}" rx="6" fill="#dbeafe" stroke="#93c5fd"></rect>`;
  const colors = ["#d1fae5","#fde68a","#fecaca","#fbcfe8","#99f6e4","#c7d2fe"];
  info.domains.forEach((d,i)=>{ const x = scale(d[1]); const w = Math.max(12, scale(d[2]) - scale(d[1])); svg += `<rect x="${x}" y="46" width="${w}" height="20" rx="6" fill="${colors[i%colors.length]}" stroke="#6b7280"></rect><text x="${x+4}" y="40" font-size="11" fill="#64748b">${d[0]} (${d[1]}-${d[2]})</text>`; });
  if(aa){ const x = scale(aa); const acc = getComputedStyle(document.body).getPropertyValue('--accent'); svg += `<line x1="${x}" y1="${30}" x2="${x}" y2="${90}" stroke="${acc}" stroke-width="2"></line><circle cx="${x}" cy="${90}" r="4" fill="${acc}"></circle><text x="${x+6}" y="${28}" font-size="12" fill="${acc}">aa ${aa}</text>`; }
  svg += `<text x="${pad}" y="${110}" font-size="11" fill="#64748b">1</text><text x="${W-pad-10}" y="${110}" font-size="11" fill="#64748b">${len}</text></svg>`;
  const w = window.open('', '_blank','width=820,height=300'); w.document.write('<!doctype html><title>'+gene+' protein</title>'+svg); w.document.close();
}

/* Pathways: compute affected-only, mark mutated nodes in red */
const gene2path = { EGFR: ['EGFR signaling','MAPK','PI3K–AKT'], KRAS: ['RAS–MAPK','MAPK','PI3K–AKT'], BRCA1: ['Homologous recombination (HR)'], ALK: ['ALK fusion signaling','MAPK','PI3K–AKT'], PTEN: ['PI3K–AKT','mTOR'] };
function computePathways(){ const map = new Map(); (report.rows||[]).forEach(r=> (gene2path[r.gene]||[]).forEach(p=> { if(!map.has(p)) map.set(p, new Set()); map.get(p).add(r.gene); })); return map; }
function pNode(id,label,x,y,active=false,mut=false){ const fill = active ? 'var(--chip)' : '#e5e7eb'; const stroke = mut? 'var(--danger)' : '#94a3b8'; const sw = mut? 2 : 1; return `<g class="node" data-id="${id}" transform="translate(${x},${y})"><rect rx="6" ry="6" width="70" height="26" fill="${fill}" stroke="${stroke}" stroke-width="${sw}"></rect><text x="35" y="17" text-anchor="middle" font-size="11" fill="#111">${label}</text></g>`; }
function mechText(gene, pathway){
  const T = {
    EGFR: {"EGFR signaling":"EGFR L858R heightens autophosphorylation and adaptor recruitment, enabling ligand-independent output. Downstream RAS–RAF–MEK–ERK and PI3K–AKT remain persistently active, sustaining proliferation and survival. Clinically, EGFR tyrosine kinase inhibitors exploit this dependence, though resistance frequently arises through secondary alterations or pathway bypass.",
            "MAPK":"The L858R substitution promotes continuous RAS activation and ERK flux independent of ligand. This enforces proliferation programs, often conferring sensitivity to EGFR-directed agents. Reactivation of ERK through downstream mutations or adaptive signaling can diminish response to single-agent therapy.",
            "PI3K–AKT":"EGFR phosphorylation augments PI3K recruitment, elevating PIP3 and AKT activation to sustain growth and survival. Crosstalk with MAPK underlies incomplete responses; dual-pathway inhibition is sometimes explored to prevent escape."},
    KRAS: {"RAS–MAPK":"KRAS G12X impairs GAP-mediated GTP hydrolysis, locking KRAS in its active GTP-bound state. Persistent RAF–MEK–ERK signaling drives proliferation and transcriptional programs of growth. Allele-selective KRAS inhibitors or downstream MEK/ERK blockade can attenuate this signaling depending on context.",
            "MAPK":"Constitutively active KRAS seeds continuous MAPK signaling downstream of RTKs. Upstream EGFR inhibition has limited effect when KRAS remains active, whereas KRAS-selective or MEK/ERK inhibitors may provide benefit in certain tumors.",
            "PI3K–AKT":"KRAS activates PI3K as a parallel effector, promoting PIP3/AKT signaling for survival and metabolism. Combined MAPK and PI3K pathway therapies are explored clinically to counter compensatory crosstalk."},
    ALK: {"ALK fusion signaling":"EML4–ALK creates a constitutively dimerized kinase that signals without ligand, activating MAPK and PI3K pathways. ALK inhibitors suppress this oncogenic driver; resistance emerges via kinase-domain mutations or bypass signaling, motivating next-generation inhibitors and combinations.",
           "MAPK":"Constitutive ALK activity triggers the RAF–MEK–ERK cascade. Pharmacologic ALK inhibition dampens ERK signaling, but adaptive reactivation through alternative RTKs can occur, suggesting benefit from rational combinations.",
           "PI3K–AKT":"ALK output engages PI3K–AKT to support survival and metabolic programs. Combinations with PI3K or mTOR inhibitors can delay resistance in preclinical and clinical settings."},
    PTEN: {"PI3K–AKT":"PTEN loss-of-function reduces PIP3 dephosphorylation, increasing AKT activation and mTOR signaling. This can blunt responses to upstream RTK inhibitors and favors AKT/mTOR-directed strategies depending on co-mutations and tumor lineage.",
            "mTOR":"Elevated PIP3 and AKT augment mTOR activity, promoting growth and protein synthesis. Therapeutic approaches may include AKT or mTOR inhibition; efficacy varies by context and compensatory signaling."},
    BRCA1: {"Homologous recombination (HR)":"Truncating BRCA1 disrupts homologous recombination for double-strand break repair. HR deficiency increases reliance on PARP-mediated repair, sensitizing to PARP inhibitors and platinum therapies, though reversion mutations or context-specific factors influence benefit."}
  };
  return (T[gene] && T[gene][pathway]) ? T[gene][pathway] : "No mechanism text available for this gene–pathway combination in the demo.";
}
function renderPathway(name, genes){
  const W=540,H=220; let svg = `<svg viewBox="0 0 ${W} ${H}" width="100%" height="220">`;
  svg += `<defs><marker id="arr" markerWidth="8" markerHeight="6" refX="6" refY="3" orient="auto"><path d="M0,0 L0,6 L6,3 z" fill="#94a3b8"/></marker><marker id="arrR" markerWidth="8" markerHeight="6" refX="6" refY="3" orient="auto"><path d="M0,0 L0,6 L6,3 z" fill="var(--danger)"/></marker></defs>`;
  const gset = new Set(genes);
  const mut = (g)=> gset.has(g);
  function arrow(x1,y1,x2,y2,hot){ const col = hot? 'var(--danger)' : '#94a3b8'; const mk = hot? 'url(#arrR)' : 'url(#arr)'; return `<line x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}" stroke="${col}" marker-end="${mk}"/>`; }
  if(/RAS|MAPK/.test(name)){
    svg += pNode('RTK','RTK',20,40,mut('EGFR')||mut('ALK'),mut('EGFR')||mut('ALK'));
    svg += pNode('RAS','RAS',130,40,mut('KRAS'),mut('KRAS'));
    svg += pNode('RAF','RAF',240,40,false,false);
    svg += pNode('MEK','MEK',350,40,false,false);
    svg += pNode('ERK','ERK',460,40,false,false);
    svg += arrow(90,53,130,53, mut('EGFR')||mut('ALK'));
    svg += arrow(200,53,240,53, mut('KRAS'));
    svg += arrow(310,53,350,53, mut('KRAS'));
    svg += arrow(420,53,460,53, mut('KRAS'));
  } else if(/PI3K/.test(name)){
    svg += pNode('RTK','RTK',20,40,mut('EGFR')||mut('ALK'),mut('EGFR')||mut('ALK'));
    svg += pNode('PI3K','PI3K',130,40,false,false);
    svg += pNode('AKT','AKT',240,40,false,false);
    svg += pNode('mTOR','mTOR',350,40,mut('PTEN'),mut('PTEN'));
    svg += arrow(90,53,130,53, mut('EGFR')||mut('ALK'));
    svg += arrow(200,53,240,53, mut('EGFR')||mut('ALK'));
    svg += arrow(310,53,350,53, mut('EGFR')||mut('ALK'));
    if(mut('PTEN')){ svg += `<line x1="300" y1="90" x2="240" y2="66" stroke="var(--danger)" marker-end="url(#arrR)"/>`; svg += `<text x="302" y="88" font-size="10" fill="var(--danger)">PTEN</text>`; }
  } else if(/ALK fusion/.test(name)){
    svg += pNode('ALK','EML4–ALK',20,40,mut('ALK'),mut('ALK'));
    svg += pNode('MAPK','MAPK',180,20,true,false);
    svg += pNode('PI3K','PI3K–AKT',180,80,true,false);
    svg += arrow(90,53,180,33, mut('ALK'));
    svg += arrow(90,53,180,93, mut('ALK'));
  } else if(/Homologous recombination/.test(name)){
    svg += pNode('BRCA','BRCA1/2',20,40,mut('BRCA1')||mut('BRCA2'),mut('BRCA1')||mut('BRCA2'));
    svg += pNode('HR','HR repair',150,40,false,false);
    svg += pNode('PARP','PARP',280,40,false,false);
    svg += arrow(90,53,150,53, mut('BRCA1')||mut('BRCA2'));
    svg += arrow(220,53,280,53, mut('BRCA1')||mut('BRCA2'));
  } else if(/EGFR signaling/.test(name)){
    svg += pNode('EGFR','EGFR',20,40,mut('EGFR'),mut('EGFR'));
    svg += pNode('MAPK','MAPK',140,20,true,false);
    svg += pNode('PI3K','PI3K–AKT',140,80,true,false);
    svg += arrow(90,53,140,33, mut('EGFR'));
    svg += arrow(90,53,140,93, mut('EGFR'));
  } else if(/mTOR/.test(name)){
    svg += pNode('AKT','AKT',20,40,false,false);
    svg += pNode('mTOR','mTOR',150,40,mut('PTEN'),mut('PTEN'));
    svg += arrow(90,53,150,53, mut('PTEN'));
  }
  svg += `</svg>`; return svg;
}
function buildPathwayUI(){
  const map = computePathways();
  const chips = document.getElementById('pathwayChips');
  const panel = document.getElementById('pathwayOne');
  chips.innerHTML = '';
  let first = null;
  map.forEach((genes, name)=>{
    const b = document.createElement('button');
    b.className = 'chip'; b.textContent = name; b.title = `Genes: ${Array.from(genes).join(', ')}`;
    b.onclick = ()=>{
      chips.querySelectorAll('.chip').forEach(x=>x.classList.remove('active'));
      b.classList.add('active');
      const g = Array.from(genes);
      const left = `<div>${renderPathway(name,g)}</div>`;
      const gene = g[0];
      const mech = mechText(gene, name);
      const refs = [`https://pubmed.ncbi.nlm.nih.gov/?term=${encodeURIComponent(gene+" "+name+" clinical")}&filter=years.2019-2025`];
      const right = `<aside id="mechPanel"><strong>Mechanism</strong><div id="mechText" style="margin:6px 0 8px"></div><strong>References</strong><ul id="mechRefs" style="padding-left:16px;margin:6px 0"></ul><small class="muted">Click a highlighted node for details.</small></aside>`;
      panel.innerHTML = `<div class="pflex">${left}${right}</div>`;
      document.getElementById('mechText').textContent = mech;
      document.getElementById('mechRefs').innerHTML = refs.map(u=> `<li><a href="${u}" target="_blank" rel="noopener">${u.replace('https://pubmed.ncbi.nlm.nih.gov/','PubMed: ')}</a></li>`).join('');
      panel.querySelectorAll('.node').forEach(nd=>{
        nd.addEventListener('click', ()=>{
          const label = nd.querySelector('text').textContent;
          let gene = null;
          if(/EGFR/i.test(label)) gene='EGFR';
          if(/RAS/.test(label) && g.includes('KRAS')) gene='KRAS';
          if(/ALK/.test(label)) gene='ALK';
          if(/BRCA/.test(label)) gene='BRCA1';
          if(/mTOR/.test(label) && g.includes('PTEN')) gene='PTEN';
          if(!gene){ document.getElementById('mechText').textContent='No gene-specific mechanism here.'; document.getElementById('mechRefs').innerHTML=''; return; }
          const txt = mechText(gene, name);
          document.getElementById('mechText').textContent = txt;
          const q = encodeURIComponent(gene+" "+name+" clinical"); 
          document.getElementById('mechRefs').innerHTML = [`https://pubmed.ncbi.nlm.nih.gov/?term=${q}&filter=years.2019-2025`].map(u=> `<li><a href="${u}" target="_blank" rel="noopener">${u.replace('https://pubmed.ncbi.nlm.nih.gov/','PubMed: ')}</a></li>`).join('');
        });
      });
    };
    chips.appendChild(b);
    if(!first){ first=b; }
  });
  if(first){ first.click(); }
}
buildPathwayUI();

/* QC SUMMARY (Detailed only) */
function buildQC(){
  const qc = report.qc || {};
  const s = document.getElementById('qcSummary');
  s.innerHTML = `
    <div class="label">Mean Coverage</div><div>${qc.coverage_mean||'—'}×</div>
    <div class="label">Median Coverage</div><div>${qc.coverage_median||'—'}×</div>
    <div class="label">≥20×</div><div>${qc.pct_gt20x!=null? Math.round(qc.pct_gt20x*100)+'%':'—'}</div>
    <div class="label">On-target Rate</div><div>${qc.on_target_rate!=null? Math.round(qc.on_target_rate*100)+'%':'—'}</div>
    <div class="label">Contamination</div><div>${qc.contamination!=null? (Math.round(qc.contamination*1000)/10)+'%':'—'}</div>
    <div class="label">Duplication</div><div>${qc.dup_rate!=null? Math.round(qc.dup_rate*100)+'%':'—'}</div>
    <div class="label">Insert Size (median)</div><div>${qc.insert_median||'—'} bp</div>
    <div class="label">Q30</div><div>${qc.read_q30!=null? Math.round(qc.read_q30*100)+'%':'—'}</div>`;
  const chips = document.getElementById('qcChips');
  const panel = document.getElementById('qcPanel');
  const tabs = [
    {key:'coverage', label:'Coverage'},
    {key:'ontarget', label:'On-target'},
    {key:'contam', label:'Contamination'},
    {key:'reads', label:'Read Quality'}
  ];
  chips.innerHTML='';
  tabs.forEach(t=>{
    const b=document.createElement('button'); b.className='chip'; b.textContent=t.label;
    b.onclick=()=>{
      chips.querySelectorAll('.chip').forEach(x=>x.classList.remove('active'));
      b.classList.add('active');
      if(t.key==='coverage'){
        const s = qc.samples && qc.samples[0] || {};
        panel.innerHTML = `<div><strong>Coverage distribution (toy)</strong><div class="mini"><span>mean ${s.mean||'—'}×</span><span>median ${s.median||'—'}×</span><span>≥20× ${s.pct20x!=null? Math.round(s.pct20x*100)+'%':'—'}</span><span>≥50× ${s.pct50x!=null? Math.round(s.pct50x*100)+'%':'—'}</span><span>≥100× ${s.pct100x!=null? Math.round(s.pct100x*100)+'%':'—'}</span></div></div>`;
      } else if(t.key==='ontarget'){
        panel.innerHTML = `<div><strong>On-target</strong><div class="mini"><span>rate ${Math.round((qc.on_target_rate||0)*100)}%</span><span>dup ${Math.round((qc.dup_rate||0)*100)}%</span></div></div>`;
      } else if(t.key==='contam'){
        panel.innerHTML = `<div><strong>Contamination</strong><div class="mini"><span>estimate ${(Math.round((qc.contamination||0)*1000)/10)}%</span></div><small class="muted">Toy value for demo</small></div>`;
      } else {
        panel.innerHTML = `<div><strong>Read quality</strong><div class="mini"><span>Q30 ${Math.round((qc.read_q30||0)*100)}%</span><span>insert median ${qc.insert_median||'—'} bp</span></div></div>`;
      }
    };
    chips.appendChild(b);
  });
  const first = chips.querySelector('.chip'); if(first) first.click();
}
buildQC();

/* Help modal */
const help = document.getElementById('help');
document.getElementById('helpBtn').onclick = ()=>{ help.style.display='flex'; };
document.getElementById('helpClose').onclick = ()=>{ help.style.display='none'; };
</script>
""")

# Write partials
partials = {
    "1_head.html": head_css,
    "2_appbar.html": appbar.replace("LOGO_URI", logo_data_uri),
    "3_summary.html": summary,
    "4_therapeutic.html": therapeutic,
    "5_pathways.html": pathways,
    "6_sidebar.html": sidebar,
    "7_inspector.html": inspector,
    "8_help.html": help_modal,
    "9_scripts.html": scripts.replace("REPORT_JSON", json.dumps(report)).replace("GENE_DOMAINS", json.dumps(gene_domains)),
    "z_close.html": "\n</body>\n</html>\n",
}
for fname, content in partials.items():
    with open(os.path.join(mod_dir, fname), "w", encoding="utf-8") as f:
        f.write(content)

# Stitch them together
with open(out_file, "w", encoding="utf-8") as out:
    for fname in sorted(os.listdir(mod_dir)):
        with open(os.path.join(mod_dir, fname), "r", encoding="utf-8") as f:
            out.write(f.read())

print("Built:", out_file)

# Compute quick stats and checksum
def md5sum(path, chunk=65536):
    h = hashlib.md5()
    with open(path, "rb") as fp:
        while True:
            b = fp.read(chunk)
            if not b: break
            h.update(b)
    return h.hexdigest()

size = os.path.getsize(out_file)
checksum = md5sum(out_file)
modules = sorted(os.listdir(mod_dir))

print("Modules included (%d):" % len(modules))
for m in modules:
    print(" -", m)
print("Output size:", size, "bytes")
print("MD5:", checksum)
print("Ready at:", out_file)