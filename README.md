# Molecular Cloning Validation Skills

Reusable AI-agent skills for two-fragment Gibson assembly validation workflows.

The repository contains two variants of the same workflow:

- `codex/molecular-cloning-validation/`: Codex skill version.
- `claude-code/molecular-cloning-validation/`: Claude Code skill version.

Each skill includes a reusable Python pipeline script that can:

- read synthesized DNA construct tables from `.xlsx`, `.csv`, or `.tsv` files;
- generate insert GenBank files;
- generate cloned plasmid GenBank files from a recipient vector and replacement anchors or coordinates;
- optionally convert GenBank files to SnapGene `.dna` files when the SnapGene CLI is available;
- validate Sanger sequencing reads against `variable insert + fixed module` references;
- produce per-read CSV, construct-summary CSV, Markdown, and HTML reports.

## Typical Command

```bash
python scripts/cloning_validation_pipeline.py \
  --table synthesis.xlsx \
  --name-column construct_id \
  --dna-column DNA_final \
  --vector vector.dna \
  --left-anchor AACCGGTTCTAGAGCGCTATCGATGCCACCATG \
  --right-anchor GGGGGGTGGCGGGTCC \
  --left-overlap-seq AACCGGTTCTAGAGCGCTATCGATGCCACC \
  --right-overlap-seq GGGGGTGGCGGGTCC \
  --second-feature-name "fixed module" \
  --sanger-success-dir sanger/success \
  --sanger-failed-dir sanger/failed \
  --out clone_analysis
```

If exact replacement coordinates are known, use `--replace-start` and `--replace-end` instead of anchors.

## Dependencies

```bash
pip install pandas openpyxl biopython snapgene-reader
```

Native SnapGene `.dna` conversion also requires a local SnapGene installation with the SnapGene CLI.

## Notes

This repository contains only reusable workflow instructions and scripts. It does not include project-specific sequence files, Sanger data, or analysis outputs.
