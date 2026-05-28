---
name: molecular-cloning-validation
description: Use when a user needs Gibson assembly cloning validation from synthesized DNA spreadsheets, SnapGene or GenBank vectors, and Sanger sequencing reads, including generation of SnapGene maps and clone pass/fail reports.
---

# Molecular Cloning Validation

## Purpose

Use this for cloning projects where variable synthesized inserts are assembled into a recipient vector and validated by Sanger sequencing. The expected deliverables are insert maps, cloned plasmid maps, and a readable validation report. It also handles a fixed second Gibson fragment, such as a tag, reporter, binding module, or linker, by extracting its reference sequence from the original vector feature.

## Required Context

Before editing or generating files, identify:

- The synthesis table and its construct-name and DNA columns.
- The recipient vector file and whether it is `.dna`, `.gb`, or `.gbk`.
- The replacement region, either exact 1-based coordinates or left/right anchor sequences.
- Any name normalization rules, such as `A -> targetA` and `B -> targetB`.
- The second fragment feature name, if it must also be validated.
- The Sanger success and failed folders, if available.

## Recommended Workflow

1. Read the synthesis table with pandas and inspect columns, row count, and DNA lengths.
2. Read the vector map. For `.dna`, use `snapgene_reader`; for GenBank, use Biopython.
3. Locate the replaced vector region deterministically. Report the 1-based interval and old length.
4. Create insert GenBank files and cloned plasmid GenBank files.
5. Convert GenBank to SnapGene `.dna` when SnapGene CLI is installed.
6. Read back generated `.dna` files with `snapgene_reader` and verify sequences match the expected insert/plasmid sequences.
7. Extract the second fragment from the original vector feature. Do not reconstruct it from memory.
8. Align Sanger reads to `expected_insert + second_fragment_reference` and summarize each colony/sample.
9. Write Markdown, HTML, and CSV reports; zip the output directory if the user needs a transferable package.

## Pipeline Script

The bundled script implements the common flow:

```bash
python scripts/cloning_validation_pipeline.py \
  --table synthesis.xlsx \
  --name-column construct_id \
  --dna-column DNA_final \
  --name-map A:targetA,B:targetB \
  --vector vector.dna \
  --left-anchor AACCGGTTCTAGAGCGCTATCGATGCCACCATG \
  --right-anchor GGGGGGTGGCGGGTCC \
  --left-overlap-seq AACCGGTTCTAGAGCGCTATCGATGCCACC \
  --right-overlap-seq GGGGGTGGCGGGTCC \
  --second-feature-name "fixed module" \
  --sanger-success-dir sanger/报告成功 \
  --sanger-failed-dir sanger/报告失败 \
  --out clone_analysis
```

If coordinates are known, replace the anchors with:

```bash
--replace-start 7252 --replace-end 7551
```

Use `--skip-snapgene` when SnapGene CLI is not available; GenBank outputs and Sanger reports still work.

## Outputs To Produce

- `01_synthesis_genbank/`: GenBank files for synthesized inserts.
- `02_synthesis_snapgene_dna/`: SnapGene `.dna` files for synthesized inserts, when possible.
- `03_cloned_plasmid_genbank/`: full cloned plasmid GenBank maps.
- `04_cloned_plasmid_snapgene_dna/`: full cloned plasmid SnapGene maps, when possible.
- `05_sanger_alignment/sanger_results.csv`: per-read differences.
- `05_sanger_alignment/sanger_summary.csv`: per-construct pass/fail and recommended sample.
- `05_sanger_alignment/sanger_report.md`: Markdown validation report.
- `05_sanger_alignment/sanger_report.html`: simple visual report with status badges.
- `05_sanger_alignment/second_fragment_reference.fasta`: reference sequence extracted from the vector feature.

## Calling Rules

- `INSERT_AND_SECOND_FRAGMENT_DNA_EXACT_PASS`: full insert and second fragment match DNA exactly.
- `CDS_AND_SECOND_FRAGMENT_DNA_EXACT_PASS`: insert CDS and second fragment match DNA exactly, but insert noncoding sequence differs.
- `PROTEIN_AND_SECOND_FRAGMENT_DNA_EXACT_PASS`: insert protein is correct and second fragment DNA is exact, but insert CDS has DNA-level differences.
- `PROTEIN_LEVEL_PASS_CONFIRM_DNA_IF_NEEDED`: proteins are correct but exact DNA is not confirmed.
- `FAIL_OR_INCOMPLETE`: relevant sequence differs or coverage is insufficient.

## Validation Requirements

- Never call a clone exact unless the relevant region is fully covered by the Sanger read and has no differences.
- Keep Sanger reads from failed-report folders separate from pass calls.
- Verify `.dna` files by reading back sequence content. SnapGene CLI may write valid files even when its process exits with warnings or a nonzero code.
- In the final response, state what was verified, where the reports are, and which clones are recommended.

## Dependencies

Typical setup:

```bash
pip install pandas openpyxl biopython snapgene-reader
```

SnapGene native `.dna` conversion requires a local SnapGene installation. On macOS the CLI is often:

```bash
/Applications/SnapGene.app/Contents/MacOS/SnapGene
```
