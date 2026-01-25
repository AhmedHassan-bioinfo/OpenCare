````markdown
# OpenCare: Reproducible Nextflow workflow for tumor–normal variant calling and review reporting

<!-- Badges -->
[![Build](https://github.com/AhmedHassan-bioinfo/OpenCare/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/AhmedHassan-bioinfo/OpenCare/actions/workflows/ci.yml)
[![License: Apache-2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18365283.svg)](https://doi.org/10.5281/zenodo.18365283)

<p align="center">
  <img src="assets/OpenCare_Nlogo.jpg" alt="OpenCare Logo" width="550"/>
</p>

OpenCare is an open-source, vendor-agnostic **NGS analysis and review-reporting workflow** built with Nextflow. It processes sequencing reads, runs standard somatic calling/annotation tools, and produces a **tumor-board–style HTML report** plus **versioned JSON exports** intended for downstream consumption.

**Design goals:** reproducibility, transparency, portability (HPC, cloud, WSL2), and easy adaptation for research workflows and review-oriented reporting.

> **Scope & safety note**
> - OpenCare is intended for **research and workflow prototyping** and is **not** a regulated diagnostic system.
> - The example evaluation uses the **HCC1395/HCC1395BL** benchmark; performance under low purity, low coverage, or FFPE artifacts requires additional evaluation.
> - Indels often require extra scrutiny/confirmation; see **Limitations** below.
> - The CNA-style plot in the HTML is **QC-only** and must not be interpreted as copy-number calls.

---

## Features

### Reporting and structured exports
- Interactive **HTML** report for case review / tumor-board discussion.
- **Versioned JSON exports** (a slim JSON driving the HTML, plus a fuller JSON for downstream use).
- Optional **FHIR/mCODE-style JSON exports** (lightweight scaffolding for integration experiments; may require local mapping and validation for any specific implementation guide).

### Workflow execution
- **Nextflow DSL2**, containerized execution (Docker/Podman/Singularity).
- Runs on laptops, HPC, cloud, and **WSL2** (Windows).

### Customization
- Tunable reporting thresholds and filters.
- Optional overlays (pathway mapping, gene domains).
- Optional enrichment hooks (see *Offline mode* below).

---

## Offline mode (what it means here)

OpenCare can run in restricted environments once **containers** and **offline resources** (e.g., VEP cache, CIViC snapshot) are staged locally.  
Optional online enrichments (e.g., **OncoKB**) require internet access and credentials and are **not compatible with fully offline execution**.

---

## Workflow

![Workflow Diagram](assets/OpenCare%20Workflow%20diagram.jpg)

---

## Prototype Video

▶️ Watch the prototype: [Demo Video on YouTube](https://www.youtube.com/watch?v=jQRYuFybSV4)

<p align="center">
  <a href="https://www.youtube.com/watch?v=jQRYuFybSV4" target="_blank" rel="noopener noreferrer" aria-label="Watch the demo on YouTube">
    <picture>
      <source srcset="https://i.ytimg.com/vi_webp/jQRYuFybSV4/maxresdefault.webp" type="image/webp">
      <source srcset="https://i.ytimg.com/vi_webp/jQRYuFybSV4/hqdefault.webp" type="image/webp">
      <source srcset="https://img.youtube.com/vi/jQRYuFybSV4/maxresdefault.jpg" type="image/jpeg">
      <img src="https://img.youtube.com/vi/jQRYuFybSV4/hqdefault.jpg" alt="OpenCare demo preview" width="600" loading="lazy" decoding="async"/>
    </picture>
  </a>
</p>

---

## Live Demos

- 🔗 **HTML demo (toy data)** — *Updated: Aug 30, 2025*  
  [OpenCare demo](https://ahmedhassan-bioinfo.github.io/OpenCare/OpenCARE_demo.html)  
  *(Concept visualization using toy data.)*

- 🔗 **Dry run report (real-world dataset)** — *Updated: Sep 24, 2025*  
  [OpenCare ERR194146 report](https://ahmedhassan-bioinfo.github.io/OpenCare/OpenCare_ERR194146_report.html)  
  *(Core functionality; features continue to evolve.)*

- 🔗 **Paired tumor/normal exome (HCC1395/HCC1395BL)** — *Updated: Oct 3, 2025*  
  [OpenCare HCC1395/HCC1395BL report](https://ahmedhassan-bioinfo.github.io/OpenCare/OpenCare_Exome_Tumor_vs_Exome_Norm_report.html)  
  *Includes offline CIViC annotations and a CNA-style coverage ratio plot (QC only).*

<details>
<summary><strong>What is this sample?</strong></summary>

- **Benchmark pair:** HCC1395 (tumor) / HCC1395BL (matched normal; cell-line derived)  
- **Assay:** Paired-end exome FASTQs  
- **Purpose:** Public benchmark demonstration of the paired tumor–normal path  
- **Notes:** Cell-line benchmarks do not represent all clinical sample conditions (purity, FFPE damage, capture variability).

</details>

---

## Install

**Requirements**
- Nextflow ≥ 23
- Java 11 or 17
- Docker (or Podman/Singularity)
- Reference data (e.g., **GRCh38** FASTA + FAI; BWA index auto-built)
- Optional: **VEP cache** for offline annotation
- Optional: **Graphviz** (for DAG rendering)

```bash
git clone https://github.com/AhmedHassan-bioinfo/OpenCare
cd OpenCare
````

> **WSL2/Windows note:** keep the repo, `work/`, and `results/` on Linux (ext4) for file locking and performance. Treat `/mnt/*` as read-only inputs when possible.

---

## Quick start

### Paired tumor/normal exome (Mutect2 path)

```bash
nextflow run main.nf \
  --reads  "$HOME/OpenCare/reads/Exome_*/*_R{1,2}.fastq.gz" \
  --ref_fa "/path/to/refs/hg38.fa" \
  --tumor_id  "Exome_Tumor" \
  --normal_id "Exome_Norm" \
  --patient_id P01 \
  --outdir "$HOME/OpenCare_out/P01" \
  -w "$HOME/nxf_work" \
  -with-docker -resume
```

> **Read pairing & sample IDs:** Nextflow groups `*_R{1,2}.fastq.gz` (or `{1,2}.fastq.gz`). For tumor/normal, set `--tumor_id` and `--normal_id` to the exact sample IDs (the `sid` emitted by read grouping).

### Single-sample mode (experimental)

OpenCare includes a single-sample route for exploratory use. This route is **not evaluated in the manuscript** and should be treated as experimental unless locally validated.

```bash
nextflow run main.nf \
  --reads  "$HOME/OpenCare/reads/<SAMPLE_ID>_R{1,2}.fastq.gz" \
  --ref_fa "/path/to/refs/hg38.fa" \
  --outdir "$HOME/OpenCare_out/<SAMPLE_ID>" \
  -w "$HOME/nxf_work" \
  -with-docker -resume
```

---

## Optional components

### VEP annotation (offline cache inside Docker)

```bash
nextflow run main.nf \
  --reads  "$HOME/OpenCare/reads/<SAMPLE_ID>_R{1,2}.fastq.gz" \
  --ref_fa "/refs/hg38.fa" \
  --vep_cache "/vep_cache" \
  --outdir "$HOME/OpenCare_out/<SAMPLE_ID>" \
  -w "$HOME/nxf_work" \
  -with-docker -resume
```

### Panel QC (mosdepth) and overlays

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

* **Interactive HTML:** `results/<SAMPLE_ID>/OpenCare_<SAMPLE_ID>_report.html`
* **Structured JSON exports:**

  * Slim report JSON (drives HTML rendering)
  * Full report JSON
  * Optional mCODE/FHIR-style bundles (scaffold; validate locally for any target spec)
* **VCF/CRAM/QC:**

  * (VEP-)VCF + index
  * Sorted/mark-dup CRAM + CRAI
  * FastQC/fastp summaries; mosdepth QC (if `--panel_bed`)
* **Execution metadata:** Nextflow trace, logs, timeline, provenance

---

## Configuration tips (especially WSL2)

* Use Linux paths for `-w` and `--outdir`:

  ```bash
  -w "$HOME/nxf_work"  --outdir "$HOME/OpenCare_out"
  ```
* Inputs on `/mnt/*` are fine as read-only. Avoid writing results to Windows mounts to prevent I/O and locking errors.

---

## Common parameters (cheat sheet)

| Param                                      | Purpose                               | Example                     |
| ------------------------------------------ | ------------------------------------- | --------------------------- |
| `--reads`                                  | FASTQ glob with R1/R2 brace expansion | `reads/S1_R{1,2}.fastq.gz`  |
| `--ref_fa`                                 | Reference FASTA (FAI auto-generated)  | `/refs/hg38.fa`             |
| `--tumor_id`, `--normal_id`                | Enable TN path with Mutect2           | `Exome_Tumor`, `Exome_Norm` |
| `--vep_cache`                              | Enable VEP offline cache              | `/vep_cache`                |
| `--panel_bed`                              | Panel regions for mosdepth QC         | `/refs/panel_targets.bed`   |
| `--patient_id`                             | Identifier embedded in JSON exports   | `P01`                       |
| `--assay`                                  | Label: `panel`, `wes`, `wgs`          | `wes`                       |
| `--pathway_db`, `--gene_domains`           | JSON overlays for HTML                | `resources/*.json`          |
| `--region` / `--targets_bed`               | Restrict calling                      | `chr7:55M-56M` / BED        |
| `--max_cpus`, `--max_mem`, `--bwa_threads` | Resource knobs                        | `8`, `24 GB`, `4`           |

> **OncoKB (optional):** set `ONCOKB_TOKEN` to enable where applicable. This requires internet access and a valid token.

---

## Project structure (top-level)

```
OpenCare/
├── main.nf                 # Nextflow workflow (DSL2)
├── HTML/
│   └── build_html.py       # HTML builder
├── resources/
│   ├── pathway_gene_map_v1.json
│   └── gene_domains.json
├── assets/                 # Logos, diagrams, screenshots
├── results/                # Default publish dir (configurable)
├── work/                   # Nextflow work dir (recommend on ext4)
└── test_data/              # Small validation datasets (optional)
```

---

## Troubleshooting

* **WSL2 locking / cross-device link issues**
  Keep `-w` and `--outdir` under Linux home (ext4). Treat `/mnt/*` inputs as read-only.

* **Missing `--reads` / `--ref_fa`**
  Both are required for sequencing input mode; provide valid Linux paths.

---

## Limitations (important)

* **Validation scope:** Example benchmarking focuses on a single public tumor–normal pair (HCC1395/HCC1395BL). Additional datasets and titration experiments are needed to characterize performance across coverage/purity/FFPE conditions.
* **Indels:** Indel calls can carry a higher false-positive burden and may require stricter filters and/or orthogonal confirmation before interpretation.
* **Purity/ploidy:** Tumor purity and ploidy are not estimated in the current release; VAF-based interpretation is therefore unadjusted and context-dependent.
* **CNA-style plot:** The coverage ratio plot is QC-only (unsegmented and uncorrected) and must not be interpreted as copy-number calls.

---

## Roadmap (selected)

* Evaluation on additional public datasets and titration/downsampling experiments.
* Optional integration with dedicated CNA and purity/pIoidy tools.
* Reporting improvements (warnings when report caps are reached; clearer provenance per section).

---

## 📢 Release Announcements

### v0.0.2 (October 3, 2025)

* Integration of the **MSK Cancer Hotspot database** for hotspot-level annotation.

*(See full history in [CHANGELOG.md](CHANGELOG.md))*

---

## Contributing

Contributions welcome — see `CONTRIBUTING.md` for coding style and issue labels.

---

## License

Apache License 2.0

---

## Citation

If you use OpenCare in your work, please cite:

> Ahmed Hassan. **OpenCare: reproducible Nextflow workflow for tumor–normal variant calling and review reporting.** GitHub Repository, 2025.

```
::contentReference[oaicite:0]{index=0}
```

