#!/usr/bin/env python3
"""
Build OpenCARE PRO modular HTML wired to Nextflow outputs

- Reads a required report.json and optional mcode_bundle.json
- Keeps the original UI/sections; only swaps in real data + search chips + panel QC
- Outputs: <outdir>/OpenCARE2_report.html

Optional extras
- --pathway-db points to a JSON file mapping genes -> pathways (fallback to bundled resources)
- --gene-domains points to a JSON with {GENE: {length:int, domains:[[name,beg,end],...]} }
- If ENABLE_ONLINE=1 in the environment, we’ll fetch missing gene domain maps from UniProt
"""

from pathlib import Path
from textwrap import dedent
import argparse
import json
import os
import time
import urllib.parse
import urllib.request

# --- Pathway DB helpers (permanent embed + legacy fix) ---
import re as _re

def _read_json_smart(p):
    if not p:
        return None
    p = Path(p)
    return json.loads(p.read_text(encoding="utf-8")) if p.is_file() else None

# --- Pathway DB helpers (robust, explicit-file-first) ---
def load_pathway_db_from_args(args):
    """
    Prefer an explicit --pathway-db file if provided and readable.
    Otherwise fall back to resources; else return [].
    """
    # 1) Explicit argument
    p = getattr(args, "pathway_db", None)
    if p:
        p = Path(p)
        try:
            if p.is_file():
                obj = json.loads(p.read_text(encoding="utf-8"))
                print(f"[pathway] loaded explicit file: {p} "
                      f"type={type(obj).__name__} "
                      f"size={(len(obj) if isinstance(obj,(list,dict)) else 'NA')}")
                return obj
            else:
                print(f"[pathway] WARN: explicit file not found: {p}")
        except Exception as e:
            print(f"[pathway] WARN: failed to read explicit file {p}: {e}")

    # 2) Fall back to bundled resources
    for c in (SCRIPT_DIR / "resources" / "pathway_gene_map_v1.json",
              Path("resources/pathway_gene_map_v1.json")):
        try:
            if Path(c).is_file():
                obj = json.loads(Path(c).read_text(encoding="utf-8"))
                print(f"[pathway] loaded resource file: {c} "
                      f"type={type(obj).__name__} "
                      f"size={(len(obj) if isinstance(obj,(list,dict)) else 'NA')}")
                return obj
        except Exception as e:
            print(f"[pathway] WARN: failed to read resource {c}: {e}")

    print("[pathway] no DB found; using empty []")
    return []

def patch_pathway_placeholder(html_path: Path, pathway_db):
    """
    For legacy HTML that still contains 'const pathwayDB = [];',
    replace it once with the real JSON.
    """
    html_path = Path(html_path)
    s = html_path.read_text(encoding="utf-8")
    if _re.search(r"const\s+pathwayDB\s*=\s*\[\];", s):
        s2 = _re.sub(
            r"const\s+pathwayDB\s*=\s*\[\];",
            "const pathwayDB   = " + json.dumps(pathway_db, ensure_ascii=False) + ";",
            s,
            count=1
        )
        html_path.write_text(s2, encoding="utf-8")
        try:
            size = len(pathway_db) if isinstance(pathway_db, list) else len(pathway_db.keys())
        except Exception:
            size = "?"
        print(f"Patched {html_path} with pathway DB (size: {size})")

# ------------------ Static assets ------------------
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


# ------------------ HTML partials ------------------
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
.brand{display:flex;align-items:center;gap:10px}
.brand img{height:28px;display:block;filter:drop-shadow(0 1px 1px rgba(0,0,0,.55)) brightness(1.2) contrast(1.35) saturate(1.08)}
[data-theme="dark"] .brand img{filter:drop-shadow(0 1px 1px rgba(0,0,0,.85)) brightness(1.25) contrast(1.25)}
.brand .name{font-weight:900;letter-spacing:.2px;font-size:14px}
.appbar .seg{display:flex;border:1px solid rgba(255,255,255,.35);border-radius:12px;overflow:hidden}
.appbar .seg button{border:0;background:transparent;padding:8px 12px;color:#fff;opacity:.9;cursor:pointer}
.appbar .seg button.active{background:rgba(255,255,255,.18);opacity:1}
.appbar .btn{background:transparent;border:1px solid rgba(255,255,255,.35);color:#fff;border-radius:10px;padding:8px 12px;cursor:pointer}

.container{display:block;padding:12px 2vw;max-width:min(2000px,96vw)}
.card{background:var(--card);border:1px solid var(--border);border-radius:14px;box-shadow:var(--shadow);padding:14px;margin-right:78px}
.card h2{margin:0 0 8px}.section-title{display:flex;align-items:center;justify-content:space-between}.section-title small{color:var(--muted)}
.kv{display:grid;grid-template-columns:minmax(220px, 22vw) 1fr;gap:8px 14px}.kv .label{color:var(--muted)}.sep{height:1px;background:var(--border);margin:12px 0}
/* QC advice tiles */
.qcboxes{
  display:grid; grid-template-columns:repeat(auto-fill,minmax(260px,1fr));
  gap:12px; margin-top:8px;
}
.qcbox{
  position:relative;
  border:1px solid var(--border);
  border-radius:12px;
  padding:14px 10px 10px;
  background:var(--card);

  /* use shared vars so ring + dot always match */
  --qc-base: var(--ok);                     /* default; overridden per state */
  --qc-ring: color-mix(in oklab, var(--qc-base) 58%, transparent);
  box-shadow: 0 0 0 2px var(--qc-ring) inset;
}

/* states define the base color only */
.qcbox.good { --qc-base: var(--ok); }
.qcbox.acc  { --qc-base: var(--warn); }
.qcbox.bad  { --qc-base: var(--danger); }
/* put Acceptable at top-right; others can stay as-is */
.qcbox.acc .qcflag { right:10px; left:auto; }   /* Acceptable → right */
.qcbox.good .qcflag,
.qcbox.bad  .qcflag  { left:10px; right:auto; } /* keep Good/Concerning on the left */

/* add the dot (if your pill doesn't have it yet) */
.qcflag .dot{
  width:7px; height:7px; border-radius:50%; background:currentColor; display:inline-block;
  margin-right:6px;
}
/* --- QC pill placement: force states --- */
.qcboxes .qcbox.acc  > .qcflag { right:10px !important; left:auto !important; }  /* Acceptable → top-right */
.qcboxes .qcbox.good > .qcflag,
.qcboxes .qcbox.bad  > .qcflag { left:10px !important; right:auto !important; }  /* Good/Concerning → top-left */
/* --- QC pill placement: force states --- */
.qcboxes .qcbox.acc  > .qcflag { right:10px !important; left:auto !important; }  /* Acceptable → top-right */
.qcboxes .qcbox.good > .qcflag,
.qcboxes .qcbox.bad  > .qcflag { left:10px !important; right:auto !important; }  /* Good/Concerning → top-left */

/* Small pill-style QC chip (like the T1/T2 badge) */
.qcflag{
  position:absolute;
  top:8px;
  left:10px;          /* default left; overridden per state below */
  display:inline-flex;
  align-items:center;
  gap:6px;
  padding:2px 8px;
  border-radius:999px;
  font-size:11px;
  font-weight:600;
  line-height:1.1;
  letter-spacing:.2px;
  border:1px solid currentColor;
  background: color-mix(in oklab, var(--card) 92%, transparent);
}

/* state colors + side placement */
.qcbox.good .qcflag { color: var(--ok);    left:10px; right:auto; }
.qcbox.acc  .qcflag { color: var(--warn);  right:10px; left:auto; }   /* Acceptable → top-right */
.qcbox.bad  .qcflag { color: var(--danger);left:10px; right:auto; }

/* the little status dot */
.qcflag .dot{
  width:7px; height:7px; border-radius:50%;
  background:currentColor; display:inline-block; margin-right:6px;
}


/* body sits under the floating pill */
.qcbox .body{ margin-top:0; }

/* --- Two-sample QC layout (scoped) --- */
#qcPanel .qc2grid{
  display:grid;
  grid-template-columns:repeat(auto-fit,minmax(320px,1fr));
  gap:12px;
}
#qcPanel .qcCard{
  border:1px solid var(--border);
  border-radius:12px;
  padding:12px;
  background:var(--card);
  box-shadow:var(--shadow);
}
#qcPanel .qcCard h4{
  margin:0 0 8px;
  display:flex;align-items:center;justify-content:space-between;
}
#qcPanel .qcMini{
  display:grid;grid-template-columns:1fr 1fr;gap:6px 10px;font-size:13px;margin-bottom:8px;
}
#qcPanel .qcMini .label{color:var(--muted)}
/* reuse your existing .qcboxes, .qcbox, .qcflag from above */


/* Clinical comment body formatting */
.cBody{border-left:3px solid color-mix(in oklab,var(--accent) 65%,transparent);padding-left:10px}
.cBody>div{margin:4px 0}

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
table{border-collapse:collapse;width:100%}
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

/* Evidence layout */
.evwrap{display:flex;flex-wrap:wrap;gap:10px;overflow:visible}
.evgroup{flex:1 1 260px;min-width:240px;border:1px dashed var(--border);border-radius:10px;padding:8px 10px;background:var(--zebra)}
.evlabel{display:block;margin-bottom:6px;font-size:11px;color:var(--muted);letter-spacing:.4px}
.evchips{display:flex;flex-wrap:wrap;gap:8px}
.evchips a{display:inline-block;white-space:nowrap}
.evgroup .subtabs{display:flex;flex-wrap:wrap;gap:8px}

.mini{display:inline-flex;gap:6px;flex-wrap:wrap}
.mini span{border:1px solid var(--border);border-radius:8px;padding:2px 6px;font-size:12px;background:var(--zebra)}

#pathwayChips,#qcChips{display:flex;flex-wrap:wrap;gap:8px;margin-bottom:8px}
#pathwayChips .chip,#qcChips .chip{border:1px solid var(--border);border-radius:999px;padding:6px 10px;background:var(--card);color:var(--fg);cursor:pointer}
#pathwayChips .chip.active,#qcChips .chip.active{background:var(--accent);color:#fff;border-color:transparent}

.subtabs{display:flex;gap:8px;flex-wrap:wrap}
.subtabs .subtab{border:1px solid var(--border);border-radius:10px;padding:6px 10px;background:var(--card);color:var(--fg);text-decoration:none;line-height:1.2}
.subtabs .subtab:hover{background:var(--accent);color:#fff;border-color:transparent}

/* Segmented control */
.seg{display:flex;gap:0;border:1px solid var(--border);background:var(--card);border-radius:12px;overflow:hidden}
.seg button{border:0;background:transparent;padding:8px 14px;color:var(--fg);cursor:pointer;opacity:.9}
.seg button + button{border-left:1px solid var(--border)}
.seg button.active{background:var(--accent);color:#fff;opacity:1}
.appbar .seg{background:transparent;border:1px solid rgba(255,255,255,.35)}
.appbar .seg button{color:#fff}
.appbar .seg button + button{border-left:1px solid rgba(255,255,255,.35)}
.appbar .seg button.active{background:var(--accent);color:#fff}

/* Dark chip contrast */
[data-theme="dark"] #pathwayChips .chip,[data-theme="dark"] #qcChips .chip{background:#fff;color:#000;border-color:#fff}
[data-theme="dark"] #pathwayChips .chip.active,[data-theme="dark"] #qcChips .chip.active{background:#fff;color:#000;border-color:var(--accent);box-shadow:0 0 0 2px var(--accent) inset}

.pcard{border:1px solid var(--border);border-radius:12px;padding:12px;background:var(--card)}
.pflex{display:grid;grid-template-columns:2fr 1fr;gap:12px}

#inspector{position:fixed;right:0;top:0;bottom:0;width:520px;max-width:96vw;background:var(--card);border-left:1px solid var(--border);box-shadow:var(--shadow);transform:translateX(100%);transition:transform .25s ease;z-index:80}
#inspector.open{transform:translateX(0)}
.ins-header{display:flex;justify-content:space-between;align-items:center;padding:12px 14px;border-bottom:1px solid var(--border)}
.ins-body{padding:12px 14px;overflow:auto;height:calc(100% - 52px)}
.card-lite{border:1px dashed var(--border);border-radius:10px;padding:10px;margin-top:6px}
.kv2{display:grid;grid-template-columns:140px 1fr;gap:8px}


/* Clinical comment boxes */
#clinCards{display:grid;grid-template-columns:repeat(auto-fill,minmax(340px,1fr));gap:12px}
.cCard{border:1px solid var(--border);border-radius:12px;padding:12px;background:var(--card)}
.cHead{display:flex;justify-content:space-between;align-items:center;margin-bottom:8px;gap:8px}
.cGene{font-weight:700}

/* use the same tier badge component as the table */
.cType{background:transparent;color:inherit;border-radius:0;padding:0;font-size:12px}

/* subtle tier accents (inner ring + colored body rule) */
.cCard.t1{ box-shadow:0 0 0 2px color-mix(in oklab,var(--ok) 60%,transparent) inset; }
.cCard.t2{ box-shadow:0 0 0 2px color-mix(in oklab,var(--warn) 60%,transparent) inset; }
.cCard.t3{ box-shadow:0 0 0 1px color-mix(in oklab,var(--muted) 45%,transparent) inset; }

.cBody{border-left:3px solid color-mix(in oklab,var(--accent) 65%,transparent);padding-left:12px}
.cCard.t1 .cBody{ border-left-color: color-mix(in oklab,var(--ok) 70%,transparent); }
.cCard.t2 .cBody{ border-left-color: color-mix(in oklab,var(--warn) 70%,transparent); }



#help{position:fixed;inset:0;background:rgba(0,0,0,.55);display:none;align-items:center;justify-content:center;z-index:90}
#help .box{background:var(--card);color:var(--fg);border:1px solid var(--border);border-radius:14px;box-shadow:var(--shadow);padding:16px;max-width:700px;width:96vw}
kbd{border:1px solid var(--border);padding:2px 6px;border-radius:6px}

@media (min-width:1600px){
  .container{max-width:min(2200px,96vw)}
  .kv{grid-template-columns:minmax(260px, 26vw) 1fr}
  th,td{padding:12px 16px}
}
@media print{
  .appbar,.sidebar,.controls,.pagination,#inspector,#help,#geneViewer{display:none!important}
  .container{padding:0}
  .card{page-break-inside:avoid;box-shadow:none;border-color:#888}
  a[href]:after{content:" (" attr(href) ")";font-size:10px;color:#666}
}
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
    <button class="btn" id="tbBtn" title="Download tumor board (HTML)">Tumor board</button>            
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
        <div class="label" role="listitem">Clinical Diagnosis</div><div id="diag"></div>
        <div class="label" role="listitem">Patient ID</div><div id="pid"></div>
        <div class="label" role="listitem">Alternate ID</div><div id="altid"></div>
        <div class="label" role="listitem">Report Date</div><div id="rdate"></div>
        <div class="label" role="listitem">Case Type</div><div id="ctype"></div>
        <div class="label" role="listitem">Physician</div><div id="phys"></div>
        <div class="label" role="listitem">Biopsy Name</div><div id="bname"></div>
        <div class="label" role="listitem">Biopsy Details</div><div id="bdet"></div>
        <div class="label" role="listitem">Gender</div><div id="gender"></div>
        <div class="label">Arm-level CNAs</div><div id="armcn"></div>  
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
      <div id="pathwayFilterBar" class="seg" role="tablist" aria-label="Pathway filter" style="margin:6px 0 10px">
        <button id="pfAll" class="active" aria-selected="true">All</button>
        <button id="pfAffected" aria-selected="false">Affected</button>
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


# ------------------ Client-side JS (kept robust; no duplicate blocks) ------------------
scripts = dedent(r"""
<script>
// Full dataset is only for download; we render a filtered view
const report = {
  ...reportFull,
  rows: (reportFull.rows || []).filter(r => {
    const codingOk = (r.coding !== false) && /^(HIGH|MODERATE)$/i.test(String(r.impact || ''));
    const tierOk   = /^T[12]/i.test(String(r.tier || ''));
    const hasTx    = Array.isArray(r.treatments) && r.treatments.length > 0;
    return codingOk || tierOk || hasTx;
  })
};


// Expose to the window so late-loaded snippets can read them
window.reportFull = reportFull;
window.report     = report;

// Wire the JSON download button to the full JSON file on disk
const JSON_DL = "__JSON_DL__";
const jsonBtn = document.getElementById('jsonBtn');
if (jsonBtn) jsonBtn.onclick = () => { window.location.href = JSON_DL; };

const MCODE       = __MCODE_JSON__;
const geneDomains = __GENE_DOMAINS__;
const pathwayDB   = __PATHWAY_DB_JSON__;
// Genes to highlight in pathways = only T1/T2 rows
function genesForPathways(){
  return new Set((report.rows||[])
    .filter(r => /^T[12]/i.test(String(r.tier||'')))
    .map(r => String(r.gene||'').toUpperCase())
    .filter(Boolean));
}

/* ---------- Helpers ---------- */
function set(id, v){ const el=document.getElementById(id); if(el) el.textContent = v ?? ''; }
function escapeHtml(s){return (s??'').toString().replace(/[&<>]/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;'}[c]))}
function locusOf(r){ return `${r.chrom}:${r.pos} ${r.ref}>${r.alt}`; }
function assocIcon(a){ if(/Benefit/i.test(a)) return `<span style="color:var(--ok)">●</span> ${escapeHtml(a)}`; if(/Resistance/i.test(a)) return `<span style="color:var(--danger)">●</span> ${escapeHtml(a)}`; return escapeHtml(a||'—'); }
function tierBadge(tier){ const t=(tier||'').toUpperCase(); const color = /T1/.test(t)?'var(--ok)':(/T2/.test(t)?'var(--warn)':'var(--muted)'); return `<span class="tier"><span class="dot" style="background:${color}"></span>${escapeHtml(t||'—')}</span>`; }
function isLowQuality(m){ if(!m) return false; const dp = +m.DP || 0; const qual = +m.QUAL || 0; return (dp < 30) || (qual && qual < 200); }
function tierRank(t){
  const s = String(t||'').toUpperCase();
  if (s.startsWith('T1')) return 1;
  if (s.startsWith('T2')) return 2;
  if (s.startsWith('T3')) return 3;   // treats T3 and T3+ together
  if (s.startsWith('T4')) return 4;
  return 99;
}
function cancerTypeFor(row){
  if (!row) return null;
  const t = Array.isArray(row.treatments)
    ? row.treatments.find(x => x.cancer_type || x.disease)
    : null;
  return row.cancer_type || (t && (t.cancer_type || t.disease)) || null;
}

function cancerTypePreferT1(gene){
  const G = String(gene||'').toUpperCase();
  const rows = (report.rows||[]).filter(r => String(r.gene||'').toUpperCase() === G);
  // Prefer a T1 row that actually has a cancer type
  const pick = rows.find(r => /^T1/i.test(String(r.tier||'')) && cancerTypeFor(r))
             ||  rows.find(r => cancerTypeFor(r));
  return pick ? cancerTypeFor(pick) : null;
}

/* ---------- Clinical Comments (robust) ---------- */

function renderClinicalComments(){
  const sec  = document.getElementById('clinCommentsSec');
  const wrap = document.getElementById('clinCards');
  if (!sec || !wrap) return;

  // index variants by (GENE|HGVSP) to fetch tier/therapies later
  const rowIndex = new Map();
  (report.rows||[]).forEach(r=>{
    const g = (r.gene||'').toUpperCase();
    const p = (r.hgvsp||'').replace(/\s/g,'').toUpperCase();
    if (g) rowIndex.set(`${g}|${p}`, r);
    if (g && !rowIndex.has(`${g}|`)) rowIndex.set(`${g}|`, r);
  });

  const sources = [
    report.clinical_comments,
    report.summary?.clinical_comments,
    report.interpretation?.comments,
    report.clinical?.comments,
    report.comments
  ];

  const firstNonEmpty = (...vals)=>{ for(const v of vals){ if(v!=null && String(v).trim()) return String(v).trim(); } return ''; };
  const findVariantStr = (o, text) =>
    firstNonEmpty(o.hgvsp, o.variant, o.protein_change,
                  (text||'').match(/p\.[A-Za-z][a-z]{2}\d+(?:[A-Za-z*]+)?/i)?.[0]);

  const guessGeneFromText = (text) => {
    const m = /^\s*([A-Z0-9]{2,12})(?=[\s:])/i.exec(text||''); // e.g. "AADACL3 ENSP..." or "EGFR: ..."
    return m ? m[1].toUpperCase() : '';
  };

  const items = [];
  const push = (o, hintGene) => {
    if (!o || typeof o !== 'object') return;
    const text = firstNonEmpty(o.text,o.comment,o.summary,o.message,o.body);
    if (!text) return;

    // gene & variant
    let gene = firstNonEmpty(o.gene, o.gene_symbol, o.symbol, hintGene);
    if (!gene) gene = guessGeneFromText(text) || '—';
    const variant = findVariantStr(o, text);

    // match tier/therapies from main variants table if possible
    const key = `${gene.toUpperCase()}|${(variant||'').replace(/\s/g,'').toUpperCase()}`;
    const row  = rowIndex.get(key) || rowIndex.get(`${gene.toUpperCase()}|`);
    const tier = firstNonEmpty(o.tier, o.level, row?.tier);
    // Cancer type: only for T1 cards; prefer a T1 row for this gene, but fall back if needed
    // now: always compute a cancer type if we can
    const ctype =
      (row ? cancerTypeFor(row) : null)
      || cancerTypePreferT1(gene)
      || (report.patient?.cancer || null);


    // Effect (very light NLP, plus association fallback)
    let effect = '';
    const effM = (text.match(/\b(activation|activating|loss[- ]?of[- ]?function|LoF|gain[- ]?of[- ]?function|rewiring|fusion[- ]?driven|amplification|overexpression)\b/i) || [])[0];
    effect = effM ? effM.replace(/-/g,' ') : (row?.association || '');
    
    // Therapy: pick explicit mention or from treatments
    let therapy = (text.match(/Therapy:\s*([^;.]+)/i)||[])[1] || '';
    if (!therapy && row?.treatments?.length) {
      therapy = row.treatments.map(t=>t.therapy).filter(Boolean).join(', ');
    }
   

    // Impact/type: pull the bit that describes the variant nature (e.g., “high-impact nonsense variant”)
    let impact = '';
    (text.split(/[.;]/).map(s=>s.trim())).some(s=>{
      if (/(variant|frameshift|nonsense|missense|splice|fusion|amplification|deletion|truncating)/i.test(s)) { impact = s; return true; }
      return false;
    });
    
   
             
    // VAF
    const vaf = (text.match(/VAF[^0-9%~]*~?\s*([\d.]+%)/i)||[])[1] || '';

    // references/PMIDs
    const pmids = []
      .concat(o.pmids || o.PMIDs || o.references || o.refs || [])
      .map(String).filter(Boolean);

    items.push({ gene, variant, tier, effect, therapy, impact, vaf, pmids, cancer: ctype });
  };

  const flatten = (x, hint) => {
    if (!x) return;
    if (typeof x === 'string'){ push({text:x}, hint); return; }
    if (Array.isArray(x))     { x.forEach(v => flatten(v, hint)); return; }
    if (typeof x === 'object'){
      if (x.text || x.comment || x.summary || x.message || x.body) push(x, hint);
      if (Array.isArray(x.items)) flatten(x.items, x.gene || x.gene_symbol || hint);
      Object.entries(x).forEach(([k,v])=>{
        if (/^(gene|gene_symbol|symbol|text|comment|summary|message|body|pmids?|PMIDs?|references|refs|tier|level|variant|hgvsp|protein_change)$/i.test(k)) return;
        flatten(v, x.gene || x.gene_symbol || k);
      });
    }
  };
  sources.forEach(s => flatten(s));

  // Fallback: if no curated comments, auto-compose from top T1/T2
  if (!items.length) {
  const MAX = Number(new URLSearchParams(location.search).get('cards')) || 12;
  const top = (report.rows || [])
  .filter(r => /T1|T2/i.test(String(r.tier||'')))
  .slice(0, MAX);
 

  top.forEach(r => {
    items.push({
      gene: r.gene,
      variant: r.hgvsp || r.consequence || '',
      tier: r.tier || '',
      effect: r.association || r.consequence || '',
      therapy: (r.treatments||[]).map(t => t.therapy).filter(Boolean).join(', '),
      impact: r.known_label || r.type || '',
      vaf: (r.metrics && r.metrics.VAF!=null)
             ? (Math.round(r.metrics.VAF*1000)/10)+'%'
             : '',
      pmids: (r.treatments||[]).map(t => t.pmid).filter(Boolean),
      cancer: r.cancer_type || (report.patient?.cancer || '')           
    });
  });
}
if (!items.length){ sec.style.display='none'; return; }
                 
 // de-dupe by gene+variant+impact text so near-identical lines don’t repeat
const uniq = Array.from(new Map(
  items.map(c => [`${c.gene}|${c.variant}|${c.impact}|${c.vaf}`, c])
).values());

// sort: T1 → T2 → T3/4, then gene, then variant
uniq.sort((a,b) =>
  tierRank(a.tier) - tierRank(b.tier) ||
  (a.gene||'').localeCompare(b.gene||'', undefined, {numeric:true, sensitivity:'base'}) ||
  (a.variant||'').localeCompare(b.variant||'', undefined, {numeric:true, sensitivity:'base'})
);
                 
wrap.innerHTML = uniq.map(it => {
  const tier = (it.tier || '').toUpperCase();
  const cardTierClass = tier.startsWith('T1') ? 't1' : (tier.startsWith('T2') ? 't2' : 't3');

  const refs = it.pmids?.length ? `<div class="linkify" style="margin-top:6px">
    ${it.pmids.map(p => `<a target="_blank" rel="noopener" href="https://pubmed.ncbi.nlm.nih.gov/${String(p).replace(/\D/g,'')}/">PMID ${p}</a>`).join(' · ')}
  </div>` : '';

  return `
    <div class="cCard ${cardTierClass}">
      <div class="cHead">
        <div class="cGene">${escapeHtml(it.gene || '—')}</div>
        <div class="cType">${tierBadge(tier || '—')}</div>
      </div>
      <div class="cBody">
        <div class="cRow"><strong>Variant:</strong> ${escapeHtml(it.variant || '—')}</div>
        <div class="cRow"><strong>Effect:</strong> ${escapeHtml(it.effect || '—')}</div>
        ${ it.therapy
            ? `<div class="cRow"><strong>Therapy:</strong> ${escapeHtml(it.therapy)}</div>`
            : '' }                 
        <div class="cRow"><strong>Impact/type:</strong> ${escapeHtml(it.impact || '—')}</div>
        <div class="cRow"><strong>VAF:</strong> ${escapeHtml(it.vaf || '—')}</div>
        ${ it.cancer
            ? `<div class="cRow"><strong>Cancer type:</strong> ${escapeHtml(it.cancer)}</div>`
            : '' }
      </div>
      ${refs}
    </div>`;
}).join('');
}

  
renderClinicalComments();



/* ---------- Tumor board (HTML download) ---------- */
function buildTumorBoard(){
  const pt  = report.patient || {};
  const css = `<style>
    body{font-family:system-ui,Segoe UI,Roboto,Arial;margin:0;padding:18px;color:#0f172a;line-height:1.45}
    h1{margin:0 0 8px} .muted{color:#64748b} .sec{margin-top:18px}
    table{border-collapse:collapse;width:100%} th,td{border:1px solid #e5e7eb;padding:8px;text-align:left;vertical-align:top}
    th{background:#f3f4f6}
  </style>`;

  const key = (report.rows||[]).filter(r=>/T1|T2/i.test(r.tier||''));
  const keyRows = key.map(r=>{
    const ev = (r.treatments||[]).map(t=>t.therapy||t.evidence||'').filter(Boolean).join(', ');
    return `<tr>
      <td>${escapeHtml(r.gene||'')}</td>
      <td>${escapeHtml(r.hgvsp || r.consequence || '')}</td>
      <td>${escapeHtml(r.cancer_type||'')}</td>
      <td>${escapeHtml(r.association||'')}</td>
      <td>${escapeHtml(r.tier||'')}</td>
      <td>${escapeHtml(ev||'—')}</td>
    </tr>`;
  }).join('');
(() => {
  // CSS for pill + left/right placement
  const st = document.createElement('style');
  st.textContent = `
  .qcflag{
    position:absolute; top:8px;
    display:inline-flex; align-items:center; gap:6px;
    padding:2px 8px; border:1px solid currentColor; border-radius:999px;
    font-size:11px; font-weight:600; line-height:1.1;
    background: color-mix(in oklab, var(--card) 92%, transparent);
  }
  .qcflag .dot{ width:7px; height:7px; border-radius:50%; background:currentColor; display:inline-block; margin-right:6px }
  .qcbox.good .qcflag{ color:var(--ok);    left:10px; right:auto; }
  .qcbox.bad  .qcflag{ color:var(--danger);left:10px; right:auto; }
  .qcbox.acc  .qcflag{ color:var(--warn);  right:10px; left:auto; }
  `;
  document.head.appendChild(st);

  // Create the pills if missing
  document.querySelectorAll('.qcbox').forEach(box => {
    if (box.querySelector('.qcflag')) return;
    box.style.position = 'relative';
    const pill = document.createElement('div');
    pill.className = 'qcflag';
    const label = box.classList.contains('acc') ? 'Acceptable'
                : box.classList.contains('good') ? 'Good'
                : 'Concerning';
    pill.innerHTML = `<span class="dot"></span>${label}`;
    box.prepend(pill);
  });

  console.log('Added', document.querySelectorAll('.qcflag').length, 'qcflag pills');
})();


  const recs = Array.from(new Set(
    key.flatMap(r => (r.treatments||[]))
       .filter(t => /FDA|label|Guideline|NCCN|ESMO/i.test(t.evidence||'') || /label/i.test(t.source_url||''))
       .map(t => (t.therapy||'').trim()).filter(Boolean)
  ));

  const pathways = computePathwaysFull().filter(p=>p.affected).map(p=>p.name).join(', ') || '—';

  const qc = document.getElementById('qcSummary');
  const qcLine = qc ? qc.textContent.replace(/\s+/g,' ').trim() : '—';

  const pmids = new Set();
  (report.rows||[]).forEach(r => (r.treatments||[]).forEach(t => { if (t.pmid) pmids.add(String(t.pmid).replace(/\D/g,'')); }));
  const refsHtml = pmids.size
    ? '<ol>'+Array.from(pmids).map(p=>`<li><a target="_blank" rel="noopener" href="https://pubmed.ncbi.nlm.nih.gov/${p}/">PMID ${p}</a></li>`).join('')+'</ol>'
    : '<div class="muted">None</div>';

  return `<!doctype html><meta charset="utf-8">
    <title>Tumor board – ${escapeHtml(pt.id||'')}</title>${css}
    <h1>Tumor board</h1>
    <div class="muted">${escapeHtml(pt.cancer||'')} · ${escapeHtml(pt.id||'')} · ${escapeHtml(report.meta?.generated_on||'')}</div>

    <div class="sec">
      <h2>Key actionable variants (Tier 1/2)</h2>
      <table>
        <thead><tr><th>Gene</th><th>Protein/Type</th><th>Context</th><th>Association</th><th>Tier</th><th>Evidence/Therapy</th></tr></thead>
        <tbody>${keyRows || '<tr><td colspan="6" class="muted">None</td></tr>'}</tbody>
      </table>
    </div>

    <div class="sec">
      <h2>Recommendations</h2>
      ${recs.length ? '<ul>'+recs.map(r=>`<li>${escapeHtml(r)}</li>`).join('')+'</ul>' : '<div class="muted">No label/guideline therapies in data.</div>'}
    </div>

    <div class="sec"><h2>Pathways affected</h2><div>${escapeHtml(pathways)}</div></div>
    <div class="sec"><h2>Quality & suitability</h2><div>${qcLine}</div></div>
    <div class="sec"><h2>References</h2>${refsHtml}</div>`;
}

const tbBtn = document.getElementById('tbBtn');
if (tbBtn) tbBtn.onclick = () => dl('tumor_board.html', buildTumorBoard(), 'text/html');


/* ---------- Top meta ---------- */
const meta = report.meta || {};
document.getElementById('metaVer').textContent = `${meta.report_version || ''} · ${meta.reference || ''}`;
const patient = report.patient || {};
set('diag', patient.cancer || '');
set('pid', patient.id || '');
function pickAlternateId(rep){
  const p   = rep.patient || {};
  const s   = rep.sample  || {};
  const sp  = rep.specimen|| {};
  const sum = rep.summary || {};
  const mt  = rep.meta    || {};
  const firstNonEmpty = (...vals)=>{ for(const v of vals){ const t = v==null?'':String(v).trim(); if(t) return t; } return ''; };
  return firstNonEmpty(
    rep.alt_id,                        // ← add this
    p.alt_id, p.alternate_id, p.SID, p.sid, p.sample_id, p.sampleId, p.specimen_id,
    s.sample_id, s.id, s.name,
    sp.specimen_id, sp.id, sp.name,
    rep.sample_id, sum.sample_id, mt.sample_id
  );
}

set('altid', (report.alt_id || pickAlternateId(report) || '').trim() || '—');
set('ctype', patient.age_group || '');
set('phys', patient.physician || '');
set('bname', patient.biopsy_name || '');
set('bdet', patient.biopsy_details || '');
set('gender', patient.gender || '');
document.getElementById('crumbs').textContent = `${(patient.cancer||'').trim()} · ${(patient.id||'').trim()}`;

/* ---------- Summary block ---------- */
const sum = report.summary || {};
set('tcontent', sum.tumor_content || '—');
set('microb', sum.microbial_species || '—');
set('mutsig', (sum.mutation_signatures||[]).join('; '));
set('svb', sum.sv_burden || '—');
set('armcn', sum.arm_cn_summary || '—');

/* ---------- Theme & view toggles ---------- */
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
let savedView = localStorage.getItem('viewmode') || 'detailed';
if(savedView==='clinician') savedView='summary';
if(savedView==='analyst') savedView='detailed';
setView(savedView);
if(savedTheme==='dark'){ lightBtn.classList.remove('active'); darkBtn.classList.add('active'); }
if(savedView==='detailed'){ summaryBtn.classList.remove('active'); detailedBtn.classList.add('active'); }
lightBtn.onclick = ()=>{ setTheme('light'); lightBtn.classList.add('active'); darkBtn.classList.remove('active'); };
darkBtn.onclick = ()=>{ setTheme('dark'); darkBtn.classList.add('active'); lightBtn.classList.remove('active'); };
summaryBtn.onclick = ()=>{ setView('summary'); summaryBtn.classList.add('active'); detailedBtn.classList.remove('active'); };
detailedBtn.onclick = ()=>{ setView('detailed'); detailedBtn.classList.add('active'); summaryBtn.classList.remove('active'); };

/* ---------- Scroll spy ---------- */
const navLinks = Array.from(document.querySelectorAll('.sideicons a'));
const sections = navLinks.map(a => document.querySelector(a.getAttribute('href')));
function onScrollSpy(){ let idx=0, y=window.scrollY+100; sections.forEach((sec,i)=>{ if(sec && sec.offsetTop < y) idx=i; }); navLinks.forEach((a,i)=> a.classList.toggle('active', i===idx)); }
window.addEventListener('scroll', onScrollSpy); onScrollSpy();


/* ---------- Metrics mini ---------- */
function metricsMini(m){
  if(!m) return '<small class="muted">—</small>';
  const chips=[];
  if(m.DP!=null) chips.push(`<span>DP ${m.DP}</span>`);
  if(m.VAF!=null) chips.push(`<span>VAF ${Math.round(m.VAF*1000)/10}%</span>`);
  if(m.QUAL!=null) chips.push(`<span>QUAL ${m.QUAL}</span>`);
  return `<div class="mini">${chips.join('')}</div>`;
}

/* ---------- Evidence rendering ---------- */
/* ---------- Evidence rendering ---------- */
function evidenceBlock(r){
  const tx = (r.treatments || []).map(t => {
    const name = escapeHtml(t.therapy || t.name || '');
    const ev   = escapeHtml(t.evidence || t.level || '');
    const src  = escapeHtml(t.source || '');
    const pmid = t.pmid ? String(t.pmid).replace(/\D/g,'') : '';
    const pm   = pmid ? `<a target="_blank" rel="noopener" href="https://pubmed.ncbi.nlm.nih.gov/${pmid}/">PMID</a>` : '';
    const url  = t.source_url ? `<a target="_blank" rel="noopener" href="${t.source_url}">${src || 'link'}</a>` : (src ? `<span>${src}</span>` : '');
    const chips = [name, ev].filter(Boolean).join(' · ');
    return `<span class="badge">${chips || 'Evidence'}</span> ${[pm,url].filter(Boolean).join(' · ')}`;
  });

  const g = encodeURIComponent(r.gene || '');
  const p = encodeURIComponent(r.hgvsp || '');
  const lookups = [
    `<a class="subtab" target="_blank" rel="noopener" href="https://pubmed.ncbi.nlm.nih.gov/?term=${g}%20${p}%20cancer&filter=years.5">PubMed (5y)</a>`,
    `<a class="subtab" target="_blank" rel="noopener" href="https://www.clinicaltrials.gov/search?term=${g}%20${p}">Trials</a>`
  ];

  return `<div class="evwrap">
    <div class="evgroup"><span class="evlabel">Evidence</span>
      <div class="evchips">${tx.join(' ') || '<small class="muted">—</small>'}</div>
    </div>
    <div class="evgroup"><span class="evlabel">Look up</span>
      <div class="subtabs">${lookups.join('')}</div>
    </div>
  </div>`;
}

/* ---------- Therapeutic table ---------- */
const tbody = document.querySelector('#tbl tbody');
const counts = document.getElementById('counts');
const qInput = document.getElementById('q');
let state = { q:'', sortK:'tier', sortDir:1 };
const knownMap = Object.fromEntries(
  (report.rows || []).map(r => [`${r.gene}:${r.hgvsp || ''}`, (r.known_label || '').toString()])
);
// mark Tier header as sorted on load (so the arrow matches)
const thTier = document.querySelector('thead th[data-k="tier"]');
if (thTier) {
  thTier.classList.add('sorted');
  thTier.querySelector('.arrow').textContent = state.sortDir === 1 ? '↑' : '↓';
}
function applyFilters(){
  const q = state.q.toLowerCase();
  const f = v => (v||'').toLowerCase().includes(q);

  const rows = (report.rows||[]).filter(r =>
    !q || f(r.gene)||f(r.hgvsp)||f(r.cancer_type)||f(r.association)||f(r.tier)||
    f(locusOf(r)) || (r.treatments||[]).some(t => f(t.therapy)||f(t.evidence))
  );

  const k   = state.sortK;
  const dir = state.sortDir; // 1 = asc, -1 = desc

  rows.sort((a,b) => {
    if (k === 'tier') {
      const d = tierRank(a.tier) - tierRank(b.tier);
      if (d) return dir * d; // T1→T2→T3
      const g = (a.gene||'').localeCompare(b.gene||'', undefined, {numeric:true, sensitivity:'base'});
      if (g) return dir * g;
      const p = (a.hgvsp||'').localeCompare(b.hgvsp||'', undefined, {numeric:true, sensitivity:'base'});
      return dir * p;
    }
    const va = (a[k] ?? '').toString();
    const vb = (b[k] ?? '').toString();
    const d  = va.localeCompare(vb, undefined, {numeric:true, sensitivity:'base'});
    return dir * d;
  });

  return rows;
}



function renderTable(){
  const all = applyFilters(); counts.textContent = `${all.length} variants`; tbody.innerHTML='';
  all.forEach((r,i)=>{
    const id = `var-${i+1}`; const alert = isLowQuality(r.metrics) ? `<span class="alert" title="Low-quality metrics">!</span>` : '';
    const tr=document.createElement('tr'); tr.id=id;
    tr.innerHTML = `
      <td><button class="geneBtn" data-gene="${escapeHtml(r.gene)}" data-hgvsp="${escapeHtml(r.hgvsp||'')}" title="Open gene map">${alert}${escapeHtml(r.gene||'')}</button>
          <button class="inspectIcon" data-idx="${i}" title="Inspect variant">🔎</button>
          <a href="#${id}" class="anchorlink" title="Deep link">#</a></td>
      <td>${escapeHtml(r.hgvsp||'')}</td>
      <td>${
        r.known_label
          ? `<span class="badge" title="From web annotation">${escapeHtml(r.known_label)}</span>`
          : '<small class="muted">—</small>'
      }</td>
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
for(const th of document.querySelectorAll('thead th[data-k]')){
  th.addEventListener('click', ()=>{
    const k = th.dataset.k;
    if(state.sortK===k){ state.sortDir*=-1; } else { state.sortK=k; state.sortDir=1; }
    document.querySelectorAll('thead th').forEach(x=>x.classList.remove('sorted'));
    th.classList.add('sorted');
    th.querySelector('.arrow').textContent = state.sortDir===1?'↑':'↓';
    renderTable();
  });
}
renderTable();

/* ---------- Inspector ---------- */
const insp = document.getElementById('inspector');
document.getElementById('insClose').onclick = ()=>{ insp.classList.remove('open'); insp.setAttribute('aria-hidden','true'); };
function openInspector(r){
  document.getElementById('insTitle').textContent = `${r.gene} · ${r.hgvsp || r.consequence || ''}`;
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
  const L = r.links||{};
  let linkHtml = ['ClinVar','OncoKB','IGV'].map(k=> L[k] ? `<a href="${L[k]}" target="_blank" rel="noopener">${k}</a>` : '').filter(Boolean).join(' · ');
  if (r.search){
    const add = [];
    if (r.search.pubmed) add.push(`<a href="${r.search.pubmed}" target="_blank" rel="noopener">PubMed (5y)</a>`);
    if (r.search.trials) add.push(`<a href="${r.search.trials}" target="_blank" rel="noopener">ClinicalTrials</a>`);
    if (r.search.oncokb) add.push(`<a href="${r.search.oncokb}" target="_blank" rel="noopener">OncoKB</a>`);
    if (r.search.civic)  add.push(`<a href="${r.search.civic}"  target="_blank" rel="noopener">CIViC</a>`);
    if (r.search.google) add.push(`<a href="${r.search.google}" target="_blank" rel="noopener">Google</a>`);
    linkHtml = (linkHtml? linkHtml + ' · ' : '') + add.join(' · ');
  }
  document.getElementById('insLinks').innerHTML = linkHtml || '<span class="muted">No links</span>';
  document.getElementById('insTreats').innerHTML = evidenceBlock(r);
  insp.classList.add('open'); insp.setAttribute('aria-hidden','false');
}

/* ---------- Gene domain viewer ---------- */
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

/* ---------- Pathways (DB-backed, tolerant) ---------- */
const ALL_PATHWAYS = [
  'EGFR signaling','RAS–MAPK','PI3K–AKT','mTOR','ALK fusion signaling',
  'Homologous recombination (HR)','WNT/β-catenin','TGF-β/SMAD','JAK/STAT',
  'NOTCH','Hedgehog','Hippo/YAP','NRF2/KEAP1','Telomerase/TERT','Chromatin/SWI–SNF',
  'Angiogenesis (VEGF)','Oncometabolism (IDH1/2)'
];
const pathwayFilter = { mode: 'all' };

function buildPathwayIndex(){
  const idx = new Map();
  const pick = (obj, ...keys) => { for (const k of keys) if (obj && obj[k]!=null) return obj[k]; return null; };

  if (Array.isArray(pathwayDB)) {
    pathwayDB.forEach(r => {
      const g = String(pick(r, 'gene_symbol','GENE','Gene','gene') || '').toUpperCase().trim();
      const raw = pick(r, 'pathways','pathway(s)','PATHWAYS','Pathways','pathway');
      const ps  = Array.isArray(raw) ? raw : String(raw||'').split(/[;,]\s*/).filter(Boolean);
      ps.forEach(p => { if(!idx.has(p)) idx.set(p, new Set()); idx.get(p).add(g); });
    });
  } else if (pathwayDB && typeof pathwayDB === 'object') {
    for (const [gene, val] of Object.entries(pathwayDB)) {
      const g = String(gene).toUpperCase().trim();
      const ps = Array.isArray(val) ? val : String(val||'').split(/[;,]\s*/).filter(Boolean);
      ps.forEach(p => { if(!idx.has(p)) idx.set(p, new Set()); idx.get(p).add(g); });
    }
  }
  return idx;
}

function pathwayDBGeneSet(){
  const idx = buildPathwayIndex();
  const s = new Set();
  for (const genes of idx.values()) {
    (genes || new Set()).forEach(g => s.add(String(g).toUpperCase().trim()));
  }
  return s;
}

function t12DBIntersect(){
  const t12 = genesForPathways();    // already uppercased T1/T2 set
  const db  = pathwayDBGeneSet();
  const inter   = new Set([...t12].filter(g => db.has(g)));
  const onlyT12 = [...t12].filter(g => !db.has(g));
  return { inter, onlyT12, t12Count: t12.size, dbCount: db.size };
}
// Expose a one-shot debugger
window.debugPathwayOverlap = () => {
  const { inter, onlyT12, t12Count, dbCount } = t12DBIntersect();
  console.info('[Pathways overlap]',
    { t12Count, dbCount, overlap: inter.size,
      overlapFirst: [...inter].sort().slice(0,20),
      t12NotInDBFirst: onlyT12.sort().slice(0,20) });
  alert(`T1/T2=${t12Count} · DB=${dbCount} · overlap=${inter.size}`);
};

function computePathwaysWithGenes(geneSet){
  const idx   = buildPathwayIndex();
  const names = Array.from(new Set([...(ALL_PATHWAYS||[]), ...idx.keys()])).sort();
  return names.map(name => {
    const genesInPath = idx.get(name) || new Set();
    const hits = Array.from(genesInPath).filter(g => geneSet.has(g));
    return { name, genes: hits, affected: hits.length > 0 };
  });
}

function computePathwaysFull(){            // T1/T2 only
  return computePathwaysWithGenes(genesForPathways());
}

function buildPathwayUI(){
  const chips = document.getElementById('pathwayChips');
  const panel = document.getElementById('pathwayOne');
  chips.innerHTML = '';
  panel.innerHTML = '';

  // info strip (insert once)
  let info = document.getElementById('pathwayInfo');
  if (!info) {
    info = document.createElement('div');
    info.id = 'pathwayInfo';
    info.style.margin = '6px 0';
    const section = document.getElementById('pathways');
    section.insertBefore(info, chips);
  }
  info.innerHTML = '';

  // Intersection check: T1/T2 vs pathway DB
  const { inter, onlyT12, t12Count } = t12DBIntersect();

  if (t12Count === 0) {
    info.innerHTML = `<div class="pcard"><small class="muted">
      No Tier 1/2 variants in this report — pathway view is empty.
    </small></div>`;
    return;
  }

  if (inter.size === 0) {
    const preview = onlyT12.slice(0,5).join(', ');
    info.innerHTML = `<div class="pcard"><small class="muted">
      No overlap between current pathway database and Tier 1/2 genes.
    </small></div>`;
    return;  // don’t render empty chips
  }

  const data = computePathwaysFull()
    .filter(p => pathwayFilter.mode === 'affected' ? p.affected : true);

  let first = null;
  data.forEach(p => {
    const b = document.createElement('button');
    b.className = 'chip';
    b.innerHTML = p.affected ? `<strong>${p.name}</strong>` : p.name;
    b.onclick = () => {
      panel.innerHTML = `
        <h4 style="margin:0">${p.name}</h4>
        <div>${p.affected ? `Affected genes: ${p.genes.join(', ')}` : 'No hits in sample'}</div>`;
      chips.querySelectorAll('button').forEach(x => x.classList.remove('active'));
      b.classList.add('active');
    };
    chips.appendChild(b);
    if (!first && p.affected) first = b;
  });
  if (first) first.click();
}

buildPathwayUI();
document.getElementById('pfAll').onclick = ()=>{
  pathwayFilter.mode='all';
  document.getElementById('pfAll').classList.add('active');
  document.getElementById('pfAffected').classList.remove('active');
  buildPathwayUI();
};
document.getElementById('pfAffected').onclick = ()=>{
  pathwayFilter.mode='affected';
  document.getElementById('pfAffected').classList.add('active');
  document.getElementById('pfAll').classList.remove('active');
  buildPathwayUI();
};
/* ---------- QC metrics (tolerant) ---------- */
/* ---------- QC metrics (tolerant) ---------- */
function buildQC(){
  const block = document.getElementById('qcBlock');
  const s     = document.getElementById('qcSummary');
  const chips = document.getElementById('qcChips');
  const panel = document.getElementById('qcPanel');
  if (!block || !s || !chips || !panel) return;

  block.style.display = '';      // ensure visible in Detailed view
  s.innerHTML = ''; chips.innerHTML = ''; panel.innerHTML = '';

  // ---- sample/model detection (return up to 2) ----
  function detectQcModels(rep){
    const qc = rep.qc || {};
    const picks = [];

    const tumor  = qc.tumor || qc.exome_tumor || qc.exome_tumour || null;
    const normal = qc.normal|| qc.exome_norm  || qc.germline     || null;
    if (tumor  && typeof tumor  === 'object') picks.push({key:'tumor',  label:'Tumor',  root:tumor});
    if (normal && typeof normal === 'object') picks.push({key:'normal', label:'Normal', root:normal});

    // If nothing explicit, try generic children/arrays
    if (!picks.length){
      if (Array.isArray(qc.samples) && qc.samples.length){
        qc.samples.slice(0,2).forEach((s,i)=>{
          picks.push({key:`sample${i+1}`, label:s.name||s.id||`Sample ${i+1}`, root:s});
        });
      } else {
        const kids = Object.keys(qc).filter(k => qc[k] && typeof qc[k]==='object');
        kids.slice(0,2).forEach((k,i)=>{
          const nice = k.replace(/_/g,' ').replace(/\b\w/g,c=>c.toUpperCase());
          picks.push({key:k, label:nice || `Sample ${i+1}`, root:qc[k]});
        });
      }
    }
    // Fallback: single model = the qc root
    if (!picks.length) picks.push({key:'single', label:'QC', root: qc});
    return picks.slice(0,2);
  }

  const models = detectQcModels(report);

  // Brief listing in the (now minimal) summary line
  s.innerHTML = `<div class="label">Samples</div><div>${
    models.map(m => m.label).join(', ')
  }</div>`;

  // ---- tolerant accessors (scoped) ----
  const isNullishStr = v => typeof v === 'string' && /^(?:na|n\/a|null|none|nan|undefined)$/i.test(v.trim());
  const num  = v => { if (v==null || isNullishStr(v) || String(v).trim()==='') return null;
                      const n = +String(v).replace(/[, ]/g,''); return Number.isFinite(n)? n : null; };
  const frac = v => { const n=num(v); if (n==null) return null; return n>1 ? n/100 : n; };

  // grade + advice text (same semantics you used)
  const grade = (val, {good, acc}) => {
    if (val == null) return {cls:'bad', label:'Concerning'};
    if (val >= good[0] && val < good[1]) return {cls:'good', label:'Good'};
    if (val >= acc[0]  && val < acc[1])  return {cls:'acc',  label:'Acceptable'};
    return {cls:'bad', label:'Concerning'};
  };

  // normalize a single model into metrics + advice boxes
  function normalizeOne(root){
    const qcRoots = [root, report.summary?.qc || {}, report.metrics || {}];
    const panelRoots = [
      ...(qcRoots.map(q => q.panel || {})),
      report.panel || {},
      root?.panel_coverage || {},
      report.qc?.panel_coverage || {}
    ];

    const pick = (roots, keys) => {
      for (const r of roots) for (const k of keys){
        if (r && r[k]!=null && !isNullishStr(r[k])) return r[k];
      } return null;
    };

    const mean   = num (pick([...panelRoots,...qcRoots], ['coverage_mean','mean_coverage','avg_coverage','average_coverage','mean_cov']));
    const median = num (pick([...panelRoots,...qcRoots], ['coverage_median','median_coverage','med_coverage','median_cov']));
    const pct20  = frac(pick([...panelRoots,...qcRoots], ['pct_gt20x','ge20_pct','pct_20x','gt20x_pct','pct_ge_20x','cov_ge20_frac']));
    const pct50  = frac(pick([...panelRoots,...qcRoots], ['pct_gt50x','ge50_pct','pct_50x','gt50x_pct','pct_ge_50x','cov_ge50_frac']));
    const pct100 = frac(pick([...panelRoots,...qcRoots], ['pct_gt100x','ge100_pct','pct_100x','gt100x_pct','pct_ge_100x','cov_ge100_frac']));
    const ontr   = frac(pick([...panelRoots,...qcRoots], ['on_target_rate','on_target_rate_pct','pct_on_target','ontarget','on_target']));
    const dup    = frac(pick([...panelRoots,...qcRoots], ['dup_rate','duplication_rate','duplication_pct','pct_duplication','dup_frac']));
    const contam = frac(pick([...panelRoots,...qcRoots], ['contamination','contamination_pct','pct_contamination','freemix','freemix_pct']));
    const insert = num (pick([...panelRoots,...qcRoots], ['insert_size_median_bp','insert_median','median_insert_size_bp','ins_size_median']));
    const q30    = frac(pick([...panelRoots,...qcRoots], ['read_q30','q30','q30_pct','pct_q30','fraction_q30']));

    const boxes = [
      {title:'Contamination', val:contam, fmt:v=>v!=null?(Math.round(v*1000)/10)+'%':'—',
       g:grade(contam,{good:[0,0.02], acc:[0.02,0.05]}),
       reason:contam==null?'Not reported':(contam<0.02?'Low cross-sample signal':(contam<0.05?'Mild contamination; verify key calls':'Elevated contamination; orthogonal confirmation suggested'))},
      {title:'Mean coverage', val:mean, fmt:v=>v!=null? v+'×':'—',
       g:grade(mean,{good:[300,Infinity], acc:[150,300]}),
       reason:mean==null?'Not reported':(mean>=300?'Ample depth for somatic calling':(mean>=150?'Adequate depth; borderline for low-VAF':'Insufficient depth for reliable low-VAF detection'))},
      {title:'On-target rate', val:ontr, fmt:v=>v!=null? Math.round(v*100)+'%':'—',
       g:grade(ontr,{good:[0.80,Infinity], acc:[0.65,0.80]}),
       reason:ontr==null?'Not reported':(ontr>=0.80?'Efficient capture':(ontr>=0.65?'Acceptable capture; some off-target reads':'Poor capture; coverage uneven'))},
      {title:'Read quality (Q30)', val:q30, fmt:v=>v!=null? Math.round(v*100)+'%':'—',
       g:grade(q30,{good:[0.90,Infinity], acc:[0.85,0.90]}),
       reason:q30==null?'Not reported':(q30>=0.90?'Base quality high':(q30>=0.85?'Moderate base quality':'Sub-optimal base quality; variant filtering impacted'))},
    ];

    return {
      mean, median, pct20, pct50, pct100, ontr, dup, contam, insert, q30, boxes
    };
  }

  // ---- render per-sample cards (2-col responsive) ----
  const grid = document.createElement('div');
  grid.className = 'qc2grid';
  panel.appendChild(grid);

  models.forEach(m => {
    const M = normalizeOne(m.root);

    const card = document.createElement('div');
    card.className = 'qcCard';
    card.innerHTML = `
      <h4>${m.label}</h4>
      <div class="qcMini">
        <div class="label">Mean</div><div>${M.mean!=null? M.mean+'×':'—'}</div>
        <div class="label">Median</div><div>${M.median!=null? M.median+'×':'—'}</div>
        <div class="label">≥20×</div><div>${M.pct20!=null?  Math.round(M.pct20*100)+'%':'—'}</div>
        <div class="label">≥50×</div><div>${M.pct50!=null?  Math.round(M.pct50*100)+'%':'—'}</div>
        <div class="label">≥100×</div><div>${M.pct100!=null? Math.round(M.pct100*100)+'%':'—'}</div>
        <div class="label">On-target</div><div>${M.ontr!=null?   Math.round(M.ontr*100)+'%':'—'}</div>
        <div class="label">Contamination</div><div>${M.contam!=null? (Math.round(M.contam*1000)/10)+'%':'—'}</div>
        <div class="label">Duplication</div><div>${M.dup!=null?     Math.round(M.dup*100)+'%':'—'}</div>
        <div class="label">Insert (med)</div><div>${M.insert!=null? M.insert+' bp':'—'}</div>
        <div class="label">Q30</div><div>${M.q30!=null? Math.round(M.q30*100)+'%':'—'}</div>
      </div>
      <div class="qcboxes">
        ${M.boxes.map(it => `
          <div class="qcbox ${it.g.cls}">
            <div class="qcflag"><span class="dot"></span>${it.g.label}</div>
            <div class="body">
              <div><strong>${it.title}:</strong> ${it.fmt(it.val)}</div>
              <div><small class="muted">${it.reason}</small></div>
            </div>
          </div>
        `).join('')}
      </div>
    `;
    grid.appendChild(card);
  });

  // Optional: simple tabs reused as “info chips” (no-op content; kept for parity)
  chips.innerHTML = `
    <button class="chip active" disabled>Coverage</button>
    <button class="chip" disabled>On-target</button>
    <button class="chip" disabled>Contamination</button>
    <button class="chip" disabled>Read quality</button>`;
}

/* ---------- Downloads & Help ---------- */
function dl(filename, text, mime='application/json'){
  const a = document.createElement('a');
  a.href = URL.createObjectURL(new Blob([text], {type:mime}));
  a.download = filename; a.click(); setTimeout(()=> URL.revokeObjectURL(a.href), 1000);
}
document.getElementById('jsonBtn').onclick = ()=> dl('report.json', JSON.stringify(report,null,2));
document.getElementById('mcodeBtn').onclick = ()=>{ if(!MCODE){ alert('No FHIR/mCODE bundle in this run.'); return; } dl('mcode_bundle.json', JSON.stringify(MCODE,null,2)); };
document.getElementById('csvBtn').onclick = ()=>{
  const header = ['gene','hgvsp','consequence','chrom','pos','ref','alt','tier'];
  const rows = (report.rows||[]).map(r=> [r.gene||'',r.hgvsp||'',r.consequence||'',r.chrom||'',r.pos||'',r.ref||'',r.alt||'',r.tier||'']);
  const csv = [header.join(','), ...rows.map(a=> a.map(x=> `"${String(x).replace(/"/g,'""')}"`).join(','))].join('\n');
  dl('variants.csv', csv, 'text/csv');
};
const help = document.getElementById('help');
document.getElementById('helpBtn').onclick = ()=>{ help.style.display='flex'; };
document.getElementById('helpClose').onclick = ()=>{ help.style.display='none'; };
</script>
""")


# ------------------ Python helpers ------------------
SCRIPT_DIR = Path(__file__).resolve().parent
# --- Online evidence helpers (Python; goes in build.py) ---
import re
import urllib.parse, urllib.request, json, os
def _http_json_post(url, payload, headers=None, timeout=25):
    hdrs = {"Content-Type": "application/json"}
    if headers: hdrs.update(headers)
    req = urllib.request.Request(url, data=json.dumps(payload).encode("utf-8"), headers=hdrs)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode("utf-8"))


def civic_query_variants(gene: str, aa_change: str):
    if not gene:
        return None
    term = f"{gene} {aa_change}".strip()
    q1 = """
    query VAR($gene: String!, $term: String!) {
      search(query: $term, entityType: VARIANTS, genes: [$gene], first: 10) {
        edges { node {
          ... on Variant {
            id name
            evidenceItems(first: 30) { nodes {
              clinicalSignificance evidenceLevel
              disease { name } therapies { name }
              source { url }
            } }
          } } }
      }
    }"""
    q2 = """
    query FALLBACK($gene: String!, $term: String!) {
      search(query: $term, entityType: EVIDENCE_ITEMS, genes: [$gene], first: 30) {
        edges { node {
          ... on EvidenceItem {
            clinicalSignificance evidenceLevel
            disease { name } therapies { name }
            source { url } variant { name }
          } } }
      }
    }"""
    try:
        headers={"Content-Type":"application/json","User-Agent":"opencare/0.1"}
        out = _http_json_post("https://civicdb.org/api/graphql",
                              {"query": q1, "variables": {"gene": gene, "term": term}}, headers=headers)
        edges = (((out or {}).get("data") or {}).get("search") or {}).get("edges") or []
        nodes = []
        for e in edges:
            v = (e.get("node") or {})
            nodes += (((v.get("evidenceItems") or {}).get("nodes")) or [])
        if nodes:
            return {"evidenceItems": {"nodes": nodes}}
        # Fallback: search evidence items directly (no type filter)
        out2 = _http_json_post("https://civicdb.org/api/graphql",
                               {"query": q2, "variables": {"gene": gene, "term": term}}, headers=headers)
        edges2 = (((out2 or {}).get("data") or {}).get("search") or {}).get("edges") or []
        nodes2 = [ (e.get("node") or {}) for e in edges2 ]
        if nodes2:
            return {"evidenceItems": {"nodes": nodes2}}
    except Exception:
        pass
    return None

def oncokb_query(gene: str, aa_change: str, token: str):
    """
    Optional: OncoKB REST call (needs token). Keep it simple.
    Returns list of dicts similar to civic_query_variants().
    """
    if not token:
        return []
    try:
        qp = urllib.parse.urlencode({"hugoSymbol": gene, "variant": aa_change})
        req = urllib.request.Request(
            f"https://www.oncokb.org/api/v1/annotate/mutations/byProteinChange?{qp}",
            headers={"Authorization": f"Bearer {token}", "Accept": "application/json"}
        )
        with urllib.request.urlopen(req, timeout=25) as r:
            d = json.loads(r.read().decode("utf-8"))
    except Exception:
        return []

    out = []
    # Very compressed mapping—adapt as needed
    assoc = (d.get("tumorTypes") or [{}])[0].get("levelOfEvidence") or ""
    for t in (d.get("treatments") or []):
        out.append({
            "variant_name": d.get("variant") or "",
            "therapy": t.get("drug") or "",
            "evidence": t.get("level") or assoc,
            "pmid": "",
            "source_url": t.get("references", [{}])[0].get("url") if t.get("references") else "",
            "source": "OncoKB"
        })
    return out

def enrich_rows_with_evidence(report: dict, oncokb_token: str):
    if os.getenv("ENABLE_EVIDENCE", "0") != "1":
        return
    for r in report.get("rows", []):
        gene = (r.get("gene") or "").strip()
        aa   = ""
        # Try to get AA change like p.L858R -> L858R
        if r.get("hgvsp"):
            m = __import__("re").search(r"[pP]\.([A-Za-z\*]?\d+[A-Za-z\*]?)", r["hgvsp"])
            if m: aa = m.group(1)
        if not gene or not aa:
            continue

        # CIViC (always)
        civic = civic_query_variants(gene, aa)
        

        # OncoKB (if token present)
        oncokb = oncokb_query(gene, aa, oncokb_token) if oncokb_token else []

        # Normalize to the schema your table expects (r.treatments / r.search)
        treats = []
        for e in civic + oncokb:
            treats.append({
                "therapy": e.get("therapy") or "",
                "evidence": e.get("evidence") or "",
                "pmid": e.get("pmid") or "",
                "source_url": e.get("source_url") or "",
                "source": e.get("source"),
            })
        # Merge (don’t clobber any pre-existing)
        r.setdefault("treatments", [])
        r["treatments"].extend(treats)

def bom_safe_load_json(path: Path):
    b = path.read_bytes()
    return json.loads(b.decode("utf-8-sig"))

def load_json_first_that_exists(candidates, what: str):
    for p in candidates:
        if p and Path(p).is_file():
            return bom_safe_load_json(Path(p))
    tried = "\n  ".join(str(c) for c in candidates if c)
    raise FileNotFoundError(f"{what} not found. Tried:\n  {tried}")

def _http_json(url, timeout=15):
    req = urllib.request.Request(url, headers={'User-Agent':'OpenCARE/1.0'})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode('utf-8'))
def enrich_rows_with_evidence(report, oncokb_token=None):
    """No-op. Rows were already enriched in MAKE_JSON_SUMMARY when ENABLE_EVIDENCE=1."""
    return

def fetch_domains_uniprot(genes):
    """
    Return {GENE: {length:int, domains:[[name,beg,end],...], source:'UniProt:ACC', isoform:'canonical'}}
    for a list of HGNC symbols (human).
    """
    out = {}
    for g in genes:
        try:
            q = f"gene_exact:{g}+AND+organism_id:9606"
            url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode({
                "query": q, "fields": "accession,reviewed", "size": 1, "format": "json", "sort": "score"
            })
            srch = _http_json(url)
            if not srch.get("results"):
                continue
            acc = srch["results"][0]["primaryAccession"]
            ent = _http_json(f"https://rest.uniprot.org/uniprotkb/{acc}.json")
            length = int(ent["sequence"]["length"])
            domains = []
            for f in ent.get("features", []):
                if f.get("type") in {"DOMAIN","TRANSMEM","TOPO_DOM","REPEAT","REGION"}:
                    loc = f.get("location") or {}
                    b = (loc.get("begin") or {}).get("value")
                    e = (loc.get("end")   or {}).get("value")
                    if b and e:
                        name = f.get("description") or f.get("type").title()
                        domains.append([name, int(b), int(e)])
            if domains:
                out[g] = {
                    "length": length,
                    "domains": sorted(domains, key=lambda x: x[1]),
                    "source": f"UniProt:{acc}",
                    "isoform": "canonical",
                }
            time.sleep(0.2)  # be nice to the API
        except Exception:
            # Ignore failures; missing genes will just lack viewers.
            pass
    return out

import re

LABEL_RE = re.compile(r'(FDA|EMA|\blabel\b|NCCN|ESMO|guideline)', re.I)
CIVIC_T1 = re.compile(r'\b(A|LEVEL[_ ]?A|B|LEVEL[_ ]?B)\b', re.I)
CIVIC_T2 = re.compile(r'\b(C|LEVEL[_ ]?C|D|LEVEL[_ ]?D)\b', re.I)

def promote_tier(row: dict):
    ts = row.get('treatments') or []
    if any(LABEL_RE.search((t.get('evidence') or '') + ' ' + (t.get('source_url') or '')) for t in ts):
        row['tier'] = 'T1'; return
    if any(CIVIC_T1.search(t.get('evidence') or '') for t in ts):
        row['tier'] = 'T1'; return
    if any(CIVIC_T2.search(t.get('evidence') or '') for t in ts):
        row['tier'] = 'T2'; return
def build_single_file(args):
    raw = json.loads(Path(args.report).read_text(encoding="utf-8"))
    mcode_obj  = (json.loads(Path(args.mcode).read_text(encoding="utf-8"))
                  if args.mcode and Path(args.mcode).exists() else None)
    gene_dom   = (json.loads(Path(args.gene_domains).read_text(encoding="utf-8"))
                  if args.gene_domains and Path(args.gene_domains).exists() else {})

# ------------------ Main ------------------
def main():
    ap = argparse.ArgumentParser(description="Build OpenCARE PRO HTML (real data)")
    ap.add_argument('--report', type=Path, required=True, help='path to report.json (real data)')
    ap.add_argument('--modules', type=Path, required=True, help='directory to write HTML partials')
    ap.add_argument('--outdir', type=Path, required=True, help='output directory for final HTML')
    ap.add_argument('--pathway-db', type=Path, default=None, help='optional pathway_gene_map_v1.json')
    ap.add_argument('--gene-domains', type=Path, default=None, help='optional gene_domains.json')
    ap.add_argument('--mcode', type=Path, default=None, help='optional FHIR/mCODE bundle JSON')
    args = ap.parse_args()
    

   

    # A stable download name for the JSON button
    dl_name = (
        f"OpenCare_{Path(args.report).stem}_report.full.json"
        if Path(args.report).name.endswith('.report.json')
        else Path(args.report).name
    )

    # -------- Load REAL inputs --------
    report = bom_safe_load_json(args.report)

    # Load pathway DB (permanent)
    pathway_db = load_pathway_db_from_args(args)

    print("[debug] pathway_db loaded:",
          type(pathway_db).__name__,
          "size=" + (str(len(pathway_db)) if isinstance(pathway_db,(list,dict)) else "NA"),
          "from", args.pathway_db)

    # Ensure dirs exist (unified names)
    args.modules.mkdir(parents=True, exist_ok=True)
    args.outdir.mkdir(parents=True, exist_ok=True)

    # One canonical HTML filename for this run
    report_stem   = Path(args.report).stem
    out_html_name = f"OpenCare_{report_stem}_report.html"
    out_file      = (args.outdir / out_html_name).resolve()

    # One canonical JSON download name (same logic used elsewhere)
    dl_name = (
        f"OpenCare_{report_stem}_report.full.json"
        if args.report.name.endswith('.report.json')
        else args.report.name
    )

    # -------- Load REAL inputs --------
    rep = json.loads(Path(args.report).read_text(encoding="utf-8"))
    bootstrap = (
        '<script id="bootstrap-json">'
        'window.reportFull = ' + json.dumps(rep, ensure_ascii=False) + ';'
        '</script>'
    )

    #  CIViC enrichment
    CIVIC = {}
    civ_path = os.getenv('CIVIC_OFFLINE', '')
    if civ_path and os.path.isfile(civ_path):
        with open(civ_path, 'r', encoding='utf-8') as f:
            CIVIC = json.load(f)
    def aa_key(hgvsp):
        m = re.search(r'[pP]\.([A-Za-z\*]?\d+[A-Za-z\*]?)', hgvsp or '')
        return ('p.'+m.group(1).upper()) if m else None

    for r in report.get('rows', []):
        g = (r.get('gene') or '').upper()
        keys = []
        if r.get('hgvsp'):
            ak = aa_key(r['hgvsp'])
            if ak: keys.append(ak.upper())
        if r.get('hgvsp'): keys.append(r['hgvsp'].upper())
        # also try variant name if you store one
        for k in keys:
            for hit in CIVIC.get((g,k), []):
                r.setdefault('treatments', []).append(hit)
        promote_tier(r)  # from §2

    # -------- Evidence enrichment (CIViC + optional OncoKB) --------
    oncokb_token = os.getenv("ONCOKB_TOKEN", "")

    # pathway DB (fallback to bundled resources; else empty list)
    

    report = bom_safe_load_json(args.report)

    pt_cancer = (report.get('patient') or {}).get('cancer') or ''
    for r in report.get('rows') or []:
        if not r.get('cancer_type'):
            disease = None
            for t in r.get('treatments') or []:
                disease = t.get('disease') or t.get('cancer_type') or disease
                if disease: break
            r['cancer_type'] = disease or pt_cancer

    rep = report  # ensure bootstrap uses enriched dict
    bootstrap = (
        '<script id="bootstrap-json">'
        'window.reportFull = ' + json.dumps(rep, ensure_ascii=False) + ';'
        '</script>'
    )

    # gene domains (fallback to bundled resources; else minimal built-in)
    default_domains = {
        "EGFR": {"length": 1210, "domains":[["TK",712,979],["L-domain",25,645]]},
        "KRAS": {"length": 188,  "domains":[["GTPase",5,167]]},
        "BRCA1":{"length": 1863, "domains":[["RING",24,64],["BRCT",1642,1736]]},
        "PTEN": {"length": 403,  "domains":[["Phosphatase",14,185],["C2",186,351]]},
        "ALK":  {"length": 1620, "domains":[["Kinase",1116,1393]]},
    }
    try:
        gene_domains = load_json_first_that_exists(
            [args.gene_domains, SCRIPT_DIR / "resources" / "gene_domains.json"],
            "Gene domains JSON"
        )
    except FileNotFoundError:
        gene_domains = default_domains

    # Optional FHIR/mCODE
    if args.mcode and args.mcode.is_file():
        try:
            mcode = bom_safe_load_json(args.mcode)
        except Exception:
            mcode = None
    else:
        mcode = None

    # Guardrails for missing keys in report
    report.setdefault("meta", {})
    report.setdefault("patient", {})
    report.setdefault("rows", [])
    report.setdefault("qc", {})

    # -------- Optional UniProt fill (must happen BEFORE injecting JS) --------
    if os.getenv("ENABLE_ONLINE", "0") == "1":
        observed_genes = sorted({
            (r.get("gene") or "").upper()
            for r in report.get("rows", []) if r.get("gene")
        })
        missing = [g for g in observed_genes if g and g not in gene_domains]
        if missing:
            fetched = fetch_domains_uniprot(missing)
            if fetched:
                gene_domains.update(fetched)
                (args.modules / "gene_domains.cache.json").write_text(
                    json.dumps(gene_domains, indent=2), encoding="utf-8"
                )

    # -------- Fill client JS placeholders --------
    scripts_filled = (
        scripts
          .replace("__REPORT_JSON__",      json.dumps(report, ensure_ascii=False))
          .replace("__MCODE_JSON__",       "null" if mcode is None else json.dumps(mcode, ensure_ascii=False))
          .replace("__GENE_DOMAINS__",     json.dumps(gene_domains, ensure_ascii=False))
          .replace("__PATHWAY_DB_JSON__",  json.dumps(pathway_db, ensure_ascii=False))
          .replace("__JSON_DL__",          dl_name)  # <- unified
    )

    # -------- Write partials (deterministic set) --------
    parts = {
      "1_head.html": head_css,
      "2_appbar.html": appbar.replace("LOGO_URI", logo_data_uri),
      "3_summary.html": summary,
      "4_therapeutic.html": therapeutic,
      "5_pathways.html": pathways,
      "6_sidebar.html": sidebar,
      "7_inspector.html": inspector,
      "8_help.html": help_modal,
      "8a_bootstrap.html": bootstrap,   # window.reportFull = <your enriched JSON>
      "9_scripts.html": scripts_filled, # ← the filled JS goes here
      "z_close.html": "\n</body>\n</html>\n",
    }

    for fname, content in parts.items():
        (args.modules / fname).write_text(content, encoding="utf-8")

    # -------- Stitch in a fixed order --------
    order = [
        "1_head.html","2_appbar.html","3_summary.html","4_therapeutic.html",
        "5_pathways.html","6_sidebar.html","7_inspector.html","8_help.html",
        "8a_bootstrap.html","9_scripts.html","z_close.html",
    ]
    # Remove stale partials not in this build
    for stale in args.modules.glob("*.html"):
        if stale.name not in order:
            try:
                stale.unlink()
            except Exception:
                pass

    with open(out_file, "w", encoding="utf-8") as out:
        for name in order:
            out.write((args.modules / name).read_text(encoding="utf-8"))
            out.write("\n")

    # Final safety pass for legacy 'const pathwayDB = [];'
    patch_pathway_placeholder(out_file, pathway_db)

    print(f"Wrote {out_file}")





if __name__ == "__main__":
    main()
