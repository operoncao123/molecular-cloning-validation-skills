---
name: molecular-cloning-validation
description: Use when validating Gibson assembly cloning projects from synthesized DNA tables, SnapGene/GenBank vector maps, and Sanger sequencing files; especially when generating insert .dna files, cloned plasmid maps, and per-clone sequence validation reports.
---

# Molecular Cloning Validation

## Overview

Use this skill to turn a synthesis table plus a recipient vector into reusable cloning deliverables: synthesized insert maps, cloned plasmid maps, Sanger validation tables, and a visual report. It is designed for workflows where one variable insert fragment is combined with an unchanged second fragment such as a tag, reporter, binding module, or linker.

## Inputs To Look For

- Insert table: `.xlsx`, `.csv`, or `.tsv` with one construct-name column and one DNA-sequence column.
- Recipient vector: preferably `.dna`; GenBank is acceptable. The vector must contain either replacement coordinates or two anchors around the replaced region.
- Optional second fragment: a named feature in the vector, for example `fixed module`, that must be checked against Sanger reads.
- Optional Sanger folders: success and failed report folders with `.ab1` and/or `.seq` files named like `1-2.*`, where the first number maps to the insert table row and the second number is the colony/sample.

## Core Workflow

1. Inspect the table columns and vector features before generating anything.
2. Locate the replacement region with explicit 1-based coordinates or by finding left/right anchors in the vector sequence.
3. Generate per-insert GenBank files and, if SnapGene CLI is available, `.dna` files.
4. Generate cloned plasmid GenBank and `.dna` maps by replacing only the intended vector region.
5. If a second-fragment feature is supplied, extract that feature directly from the original vector and use it as the reference.
6. Align Sanger reads to `insert + second_fragment`, not only to the insert, and classify each colony.
7. Verify generated `.dna` files by reading them back when `snapgene_reader` is available; do not trust SnapGene exit codes alone.
8. Deliver a concise README, CSV tables, Markdown report, HTML report, and a zip archive when useful.

## Reusable Pipeline

Prefer the bundled script for repeatable projects:

```bash
python scripts/cloning_validation_pipeline.py \
  --table synthesis.xlsx \
  --name-column construct_id \
  --dna-column DNA_final \
  --name-map A:targetA,B:targetB \
  --vector recipient_vector.dna \
  --left-anchor LEFT_VECTOR_BOUNDARY_SEQUENCE \
  --right-anchor RIGHT_VECTOR_BOUNDARY_SEQUENCE \
  --left-overlap-seq INSERT_5P_OVERLAP_SEQUENCE \
  --right-overlap-seq INSERT_3P_OVERLAP_SEQUENCE \
  --second-feature-name "fixed module" \
  --sanger-success-dir sanger/success \
  --sanger-failed-dir sanger/failed \
  --out clone_analysis
```

Use `--replace-start` and `--replace-end` instead of anchors when the exact replacement coordinates are known.

The anchor and overlap values are normal project-specific inputs. They must be inferred from the user's vector and insert design; do not reuse example sequences from another project.

## Output Structure

The script writes:

- `01_synthesis_genbank/`: insert GenBank files.
- `02_synthesis_snapgene_dna/`: insert SnapGene files when conversion succeeds.
- `03_cloned_plasmid_genbank/`: full cloned plasmid GenBank files.
- `04_cloned_plasmid_snapgene_dna/`: full cloned plasmid SnapGene files when conversion succeeds.
- `05_sanger_alignment/`: `sanger_results.csv`, `sanger_summary.csv`, `sanger_report.md`, `sanger_report.html`, and second-fragment reference FASTA.
- `generation_summary.csv`: construct lengths and generated plasmid lengths.

## Sanger Calling Rules

- `INSERT_AND_SECOND_FRAGMENT_DNA_EXACT_PASS`: full synthesized insert and second fragment are both exact.
- `CDS_AND_SECOND_FRAGMENT_DNA_EXACT_PASS`: insert CDS and second fragment are exact, but insert noncoding overlap/linker differs.
- `PROTEIN_AND_SECOND_FRAGMENT_DNA_EXACT_PASS`: insert protein is exact and second fragment DNA is exact, but insert CDS DNA differs.
- `PROTEIN_LEVEL_PASS_CONFIRM_DNA_IF_NEEDED`: protein-level pass only; recommend orthogonal confirmation if exact DNA is required.
- `FAIL_OR_INCOMPLETE`: CDS/protein differs, second fragment differs, or useful coverage is incomplete.

## Practical Checks

- Confirm that the second fragment reference comes from the original vector feature, not from memory or a pasted sequence.
- Check that Sanger reads cover the full insert CDS and the full second fragment before calling a clone exact.
- Keep failed-folder reads listed separately; do not use them for pass calls unless the user explicitly asks.
- For `.dna` generation, SnapGene on macOS may write a valid file and then exit nonzero. Verify file existence and sequence readback with `snapgene_reader` instead of relying only on the process exit code.
- If `.dna` conversion is unavailable, still generate GenBank files and tell the user SnapGene conversion was skipped.

## Dependencies

Use existing project environments when possible. Typical Python dependencies:

```bash
pip install pandas openpyxl biopython snapgene-reader
```

Native `.dna` conversion requires SnapGene CLI. Common paths are `/Applications/SnapGene.app/Contents/MacOS/SnapGene` on macOS and `C:\Program Files\SnapGene\SnapGene.exe` on Windows. CLI flags are generally the same; shell quoting and executable paths differ.
