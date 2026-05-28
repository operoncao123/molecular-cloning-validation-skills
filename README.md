# Molecular Cloning Validation Skills

Reusable AI-agent skills for molecular cloning validation workflows. The package contains one Codex skill and one Claude Code skill, both backed by the same Python pipeline.

- `codex/molecular-cloning-validation/`: Codex skill version.
- `claude-code/molecular-cloning-validation/`: Claude Code skill version.

The pipeline can generate insert maps, cloned plasmid maps, Sanger validation tables, Markdown reports, and HTML reports from a synthesis table, a recipient vector, and optional Sanger sequencing results.

## What the Pipeline Does

1. Reads a construct table from `.xlsx`, `.csv`, or `.tsv`.
2. Reads a recipient vector from SnapGene `.dna`, GenBank, or FASTA.
3. Finds the vector region to replace by either:
   - exact 1-based coordinates with `--replace-start` and `--replace-end`, or
   - project-specific flanking sequences with `--left-anchor` and `--right-anchor`.
4. Replaces that vector region with each synthesized insert.
5. Writes insert GenBank files and cloned plasmid GenBank files.
6. Converts GenBank files to SnapGene `.dna` files when the SnapGene CLI is available.
7. Extracts an optional fixed second fragment directly from a named vector feature.
8. Aligns Sanger reads to `variable insert + fixed module` references.
9. Reports per-sample differences, per-construct pass/fail calls, and recommended colonies.

## Important: Anchors and Overlaps Are Project-Specific

The following parameters are not fixed by the skill and should be different for different users or cloning designs:

- `--left-anchor`: sequence in the original vector at the left boundary of the region to replace.
- `--right-anchor`: sequence in the original vector at the right boundary of the region to replace.
- `--left-overlap-seq`: 5' homology/overlap sequence included in each synthesized insert.
- `--right-overlap-seq`: 3' homology/overlap, linker, or right-boundary sequence included in each synthesized insert.

Changing these values does not affect whether the skill can be used. They are normal input parameters. The agent should inspect each user's vector and synthesis design, then fill these values from that project.

If exact replacement coordinates are already known, prefer coordinates over anchors:

```bash
--replace-start START_POS --replace-end END_POS
```

Coordinates are 1-based and inclusive.

## Typical Command

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

Use `--skip-snapgene` when SnapGene CLI is not installed. GenBank and Sanger reports will still be produced.

## SnapGene CLI on macOS and Windows

The SnapGene CLI is part of the official SnapGene desktop application. The command-line options are effectively the same across platforms, but the executable path and shell quoting differ.

Common macOS path:

```bash
/Applications/SnapGene.app/Contents/MacOS/SnapGene
```

Common Windows paths:

```powershell
C:\Program Files\SnapGene\SnapGene.exe
C:\Program Files (x86)\SnapGene\SnapGene.exe
```

The script attempts to find `SnapGene`, `SnapGene.exe`, and the common macOS/Windows paths automatically. If auto-detection fails, pass the executable explicitly:

```bash
--snapgene-cli /Applications/SnapGene.app/Contents/MacOS/SnapGene
```

PowerShell example:

```powershell
python scripts/cloning_validation_pipeline.py `
  --table synthesis.xlsx `
  --name-column construct_id `
  --dna-column DNA_final `
  --vector recipient_vector.dna `
  --replace-start START_POS `
  --replace-end END_POS `
  --second-feature-name "fixed module" `
  --snapgene-cli "C:\Program Files\SnapGene\SnapGene.exe" `
  --out clone_analysis
```

Some SnapGene builds write valid `.dna` files even when the process exits with warnings or a nonzero code. Always verify generated `.dna` files by reading them back with `snapgene_reader` when possible.

## Output Structure

The pipeline writes:

- `01_synthesis_genbank/`: GenBank files for synthesized inserts.
- `02_synthesis_snapgene_dna/`: SnapGene `.dna` files for synthesized inserts, when conversion succeeds.
- `03_cloned_plasmid_genbank/`: full cloned plasmid GenBank maps.
- `04_cloned_plasmid_snapgene_dna/`: full cloned plasmid SnapGene maps, when conversion succeeds.
- `05_sanger_alignment/sanger_results.csv`: per-read differences.
- `05_sanger_alignment/sanger_summary.csv`: per-construct pass/fail and recommended sample.
- `05_sanger_alignment/sanger_report.md`: Markdown validation report.
- `05_sanger_alignment/sanger_report.html`: HTML visual report.
- `05_sanger_alignment/second_fragment_reference.fasta`: fixed-fragment reference extracted from the vector feature.
- `generation_summary.csv`: construct lengths and generated plasmid lengths.

## Dependencies

```bash
pip install pandas openpyxl biopython snapgene-reader
```

Native SnapGene `.dna` conversion also requires a local SnapGene desktop installation.

## Installation as Skills

For Codex, copy:

```text
codex/molecular-cloning-validation
```

into your Codex skills directory.

For Claude Code, copy:

```text
claude-code/molecular-cloning-validation
```

into your Claude Code skills directory.

## Privacy

This repository contains only reusable workflow instructions and scripts. It does not include project-specific sequence files, Sanger data, or analysis outputs.
