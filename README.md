# OpenCare: Reproducible Nextflow workflow for paired tumor–normal exome interpretation and review reporting

[![Build](https://github.com/AhmedHassan-bioinfo/OpenCare/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/AhmedHassan-bioinfo/OpenCare/actions/workflows/ci.yml)
[![License: Apache-2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18472250.svg)](https://doi.org/10.5281/zenodo.18472250)

<p align="center">
  <img src="assets/OpenCare_Nlogo.jpg" alt="OpenCare Logo" width="550">
</p>

OpenCare is an open-source, containerized **Nextflow DSL2** workflow for **paired tumor–normal (TN) exome** analysis and review-oriented reporting.  
It runs standard processing steps (QC, alignment, duplicate marking, somatic SNV/indel calling, annotation) and produces:

- an interactive **HTML review report**
- structured, versioned **JSON outputs**
- run metadata for reproducibility (logs/trace/timeline/provenance)

This repository is aligned to the manuscript scope and wording.

---

## Scope and intended use

- OpenCare is intended for **research use** and workflow prototyping.
- OpenCare is **not** a regulated diagnostic device or clinical decision system.
- The manuscript’s quantitative benchmarking is demonstrated on a public TN benchmark pair (**HCC1395/HCC1395BL**) under explicitly restricted evaluation regions.

---

## Manuscript-aligned benchmarking summary

The manuscript benchmark uses comparator outputs and harmonized region restriction/normalization for fair cross-pipeline comparison.

Key evaluation conditions (as reported in the manuscript and supplements):

- callable loci defined from **samtools depth** with **DP ≥ 10 in both tumor and matched normal**
- shared evaluation region: **HC ∩ capture targets ∩ callable**
- identical preprocessing before comparison: decomposition + left-normalization
- comparator-based evaluation (som.py / hap.py outputs archived in supplementary artifacts)

For complete evidence bundles, see the manuscript supplementary files (class-resolved outputs and harmonization/provenance package).

---

## What OpenCare reports (and what it does not)

### Included

- QC summaries for review
- somatic SNV/indel interpretation-oriented output
- HTML report for case review/triage
- structured JSON exports for downstream integration

### Important boundaries

- **Indels** may require stricter filtering and/or orthogonal confirmation.
- Current release does **not** infer tumor purity/ploidy.
- Any CNA-style coverage visualization in HTML is **QC-oriented only** and must not be interpreted as a validated CNA caller output.

---

## Quick links

### Prototype video
▶️ https://www.youtube.com/watch?v=jQRYuFybSV4

### Example reports

- Toy demo HTML:  
  https://ahmedhassan-bioinfo.github.io/OpenCare/OpenCARE_demo.html

- Dry-run report:  
  https://ahmedhassan-bioinfo.github.io/OpenCare/OpenCare_ERR194146_report.html

- TN exome demo (HCC1395/HCC1395BL):  
  https://ahmedhassan-bioinfo.github.io/OpenCare/OpenCare_Exome_Tumor_vs_Exome_Norm_report.html

- Benchmark run report (SRR7890849/SRR7890875):  
  https://ahmedhassan-bioinfo.github.io/OpenCare/SRR7890849_SRR7890875_report.html

---

## Installation

### Requirements

- Nextflow (23+ recommended)
- Java 11 or 17
- Docker / Podman / Singularity
- GRCh38 reference FASTA (with index)
- Optional local resources (e.g., VEP cache)

```bash
git clone https://github.com/AhmedHassan-bioinfo/OpenCare
cd OpenCare
```

> **WSL2 note:** keep `work/` and output directories on Linux filesystem (ext4) for better stability/performance.

---

## Quick start (paired TN exome)

```bash
nextflow run main.nf \
  --reads  "$HOME/OpenCare/reads/Exome_*/*_R{1,2}.fastq.gz" \
  --ref_fa "/path/to/refs/hg38.fa" \
  --tumor_id  "Exome_Tumor" \
  --normal_id "Exome_Norm" \
  --patient_id "P01" \
  --outdir "$HOME/OpenCare_out/P01" \
  -w "$HOME/nxf_work" \
  -with-docker -resume
```

---

## Optional components

### Offline-friendly annotation setup

OpenCare can run in restricted environments if required containers/resources are pre-staged locally.

Optional online enrichment steps (if enabled) require internet access and credentials.

### Panel-aware QC / overlays

```bash
nextflow run main.nf \
  --reads "$HOME/OpenCare/reads/<SAMPLE_ID>_R{1,2}.fastq.gz" \
  --ref_fa "/refs/hg38.fa" \
  --panel_bed "/refs/panel_targets.bed" \
  --pathway_db "$PWD/resources/pathway_gene_map_v1.json" \
  --gene_domains "$PWD/resources/gene_domains.json" \
  --outdir "$HOME/OpenCare_out/<SAMPLE_ID>" \
  -w "$HOME/nxf_work" \
  -with-docker -resume
```

---

## Outputs

- HTML report: `results/<SAMPLE_ID>/OpenCare_<SAMPLE_ID>_report.html`
- JSON outputs:
  - report-driving JSON
  - full structured JSON
- Variant / alignment / QC artifacts (depending on mode and options)
- Workflow metadata: trace, timeline, logs, and run artifacts

---

## Reproducibility notes

For manuscript-grade reproduction:

- use the same tagged release / commit
- pin container images used for the run
- use identical reference build/resources
- preserve command lines and Nextflow run metadata

The manuscript supplementary materials provide comparator outputs and harmonization/provenance evidence for the benchmarked analysis.

---

## Limitations

- Benchmark evidence in manuscript is focused on the HCC1395/HCC1395BL TN exome setting.
- Generalization to other sample types/conditions (e.g., low purity, FFPE, differing capture designs) requires separate validation.
- Interpretation remains dependent on local review policies and orthogonal confirmation practice.

---

## Contributing

Contributions are welcome. Please open an issue or PR with a clear description, reproducible steps, and expected behavior.

---

## License

Apache License 2.0

---

## Citation

If you use OpenCare, please cite the manuscript and the archived repository release (DOI badge above).
