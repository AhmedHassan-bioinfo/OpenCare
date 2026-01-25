#!/usr/bin/env python3
import argparse, gzip, sys
from collections import defaultdict

BINS = [
    (0.00, 0.05, "<5%"),
    (0.05, 0.10, "5–10%"),
    (0.10, 0.20, "10–20%"),
    (0.20, 1.01, "≥20%"),
]

def bin_label(v):
    for lo, hi, lab in BINS:
        if v >= lo and v < hi:
            return lab
    return None

def classify_variant(ref, alt):
    # Match som.py style: only SNVs (1bp->1bp) and INDELs (len differs).
    # Ignore MNVs (same length >1) so SNV/INDEL totals align with som.py.
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) != len(alt):
        return "Indel"
    return None

def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")

def parse_info(info_str):
    d = {}
    if info_str == "." or info_str == "":
        return d
    for part in info_str.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
        else:
            d[part] = True
    return d

def get_truth_vaf(info_str, tag):
    info = parse_info(info_str)
    if tag not in info:
        return None
    try:
        return float(info[tag])
    except Exception:
        return None

def get_query_vaf(fmt, sample_str, vaf_field="AF"):
    # Prefer FORMAT/AF. Fallback to FORMAT/AD.
    keys = fmt.split(":")
    vals = sample_str.split(":")
    if len(keys) != len(vals):
        return None
    m = dict(zip(keys, vals))

    if vaf_field in m and m[vaf_field] not in (".", ""):
        try:
            return float(m[vaf_field])
        except Exception:
            pass

    if "AD" in m and m["AD"] not in (".", ""):
        # AD is typically "ref,alt" for biallelic; use first alt.
        try:
            parts = m["AD"].split(",")
            if len(parts) >= 2:
                ref = float(parts[0])
                alt = float(parts[1])
                denom = ref + alt
                if denom > 0:
                    return alt / denom
        except Exception:
            pass
    return None

def load_truth_variants(path, truth_vaf_tag):
    # returns: dict key->(class, vaf)
    # key = (chrom, pos, ref, alt, class)
    truth = {}
    with open_text(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom, pos, _id, ref, alt, qual, flt, info = cols[:8]
            if "," in alt:
                # normalized VCF should be biallelic; skip multiallelic just in case
                continue
            cls = classify_variant(ref, alt)
            if cls is None:
                continue
            vaf = get_truth_vaf(info, truth_vaf_tag)
            if vaf is None:
                continue
            key = (chrom, int(pos), ref, alt, cls)
            truth[key] = vaf
    return truth

def load_query_variants(path, tumor_sample, query_vaf_field="AF"):
    # returns: dict key->vaf
    query = {}
    sample_index = None
    with open_text(path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                # samples start at col 9 (0-based index 9) => after FORMAT col 8
                if tumor_sample not in header:
                    raise SystemExit(f"[ERROR] Tumor sample '{tumor_sample}' not found in {path}")
                sample_index = header.index(tumor_sample)
                continue
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 10 or sample_index is None:
                continue
            chrom, pos, _id, ref, alt, qual, flt, info, fmt = cols[:9]
            if "," in alt:
                continue
            cls = classify_variant(ref, alt)
            if cls is None:
                continue
            tumor_str = cols[sample_index]
            vaf = get_query_vaf(fmt, tumor_str, vaf_field=query_vaf_field)
            if vaf is None:
                continue
            key = (chrom, int(pos), ref, alt, cls)
            query[key] = vaf
    return query

def f1(p, r):
    return 0.0 if (p + r) == 0 else (2.0 * p * r / (p + r))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scratch", required=True, help="som.py --scratch-prefix directory (e.g. scratch_hc_callable_dp10_for_vafbins)")
    ap.add_argument("--tumor-sample", required=True, help="tumor sample name in query VCF (e.g. Exome_Tumor)")
    ap.add_argument("--truth-vaf-tag", default="TVAF", help="INFO tag in truth VCF (default: TVAF)")
    ap.add_argument("--query-vaf-field", default="AF", help="FORMAT field in query (default: AF)")
    ap.add_argument("--label", default="REGION", help="region label (avoid spaces if you use `column -t`)")
    args = ap.parse_args()

    norm_truth = f"{args.scratch}/normalized_truth.vcf.gz"
    norm_query = f"{args.scratch}/normalized_query.vcf.gz"
    fp_query   = f"{args.scratch}/fp.vcf.gz"

    # Load full normalized truth/query
    truth = load_truth_variants(norm_truth, args.truth_vaf_tag)
    query = load_query_variants(norm_query, args.tumor_sample, args.query_vaf_field)

    truth_keys = set(truth.keys())
    query_keys = set(query.keys())

    tp_keys = truth_keys & query_keys
    fn_keys = truth_keys - query_keys
    fp_keys = query_keys - truth_keys

    # FP VAFs: use fp.vcf.gz (should match fp_keys) for binning; but compute from query too.
    # We’ll bin FP using query VAF values from `query` dict (consistent with fp_keys).
    # (If you prefer to bin using fp.vcf.gz directly, it should be identical after normalization.)
    bins = defaultdict(lambda: {"TP":0, "FP":0, "FN":0})

    # TP/FN binned by truth VAF (TVAF)
    for k in tp_keys:
        cls = k[4]
        b = bin_label(truth[k])
        if b is None:
            continue
        bins[(cls, b)]["TP"] += 1

    for k in fn_keys:
        cls = k[4]
        b = bin_label(truth[k])
        if b is None:
            continue
        bins[(cls, b)]["FN"] += 1

    # FP binned by query AF
    for k in fp_keys:
        cls = k[4]
        b = bin_label(query[k])
        if b is None:
            continue
        bins[(cls, b)]["FP"] += 1

    # Output
    print("Region\tVariantClass\tVAFstratum\tTP\tFP\tFN\tPrecision\tRecall\tF1")
    for cls in ("SNV", "Indel"):
        for _lo,_hi,lab in BINS:
            tp = bins[(cls, lab)]["TP"]
            fp = bins[(cls, lab)]["FP"]
            fn = bins[(cls, lab)]["FN"]
            p = 0.0 if (tp+fp)==0 else tp/(tp+fp)
            r = 0.0 if (tp+fn)==0 else tp/(tp+fn)
            print(f"{args.label}\t{cls}\t{lab}\t{tp}\t{fp}\t{fn}\t{p:.6f}\t{r:.6f}\t{f1(p,r):.6f}")

    # Sanity totals to stderr
    def sum_by(cls, field):
        return sum(bins[(cls, lab)][field] for _lo,_hi,lab in BINS)

    snv_tp, snv_fp, snv_fn = sum_by("SNV","TP"), sum_by("SNV","FP"), sum_by("SNV","FN")
    ind_tp, ind_fp, ind_fn = sum_by("Indel","TP"), sum_by("Indel","FP"), sum_by("Indel","FN")

    print(f"[SANITY] TP/FP/FN totals (SNV):   {snv_tp}/{snv_fp}/{snv_fn}", file=sys.stderr)
    print(f"[SANITY] TP/FP/FN totals (Indel): {ind_tp}/{ind_fp}/{ind_fn}", file=sys.stderr)
    print(f"[SANITY] NOTE: These exclude MNVs (len(ref)==len(alt)>1).", file=sys.stderr)

if __name__ == "__main__":
    main()
