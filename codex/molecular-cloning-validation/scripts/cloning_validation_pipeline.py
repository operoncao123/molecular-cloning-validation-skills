#!/usr/bin/env python3
"""Reusable Gibson/SnapGene/Sanger pipeline.

Inputs:
- Insert table (.xlsx/.csv/.tsv) with construct names and synthesized DNA.
- Vector map (.dna preferred, GenBank supported) containing the region to replace.
- Optional Sanger success/failed folders containing .ab1 or .seq files named like 1-2.*.

Outputs:
- Insert GenBank and optional SnapGene .dna files.
- Cloned plasmid GenBank and optional SnapGene .dna files.
- Sanger per-read CSV, summary CSV, Markdown report, HTML report.
"""

from __future__ import annotations

import argparse
import csv
import html
import re
import shutil
import subprocess
import sys
import warnings
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

warnings.filterwarnings("ignore", category=DeprecationWarning)

try:
    from snapgene_reader import snapgene_file_to_dict
except Exception:  # pragma: no cover - optional dependency
    snapgene_file_to_dict = None


@dataclass
class FeatureLite:
    start: int
    end: int
    strand: Optional[int]
    type: str
    label: str
    qualifiers: Dict[str, object]


@dataclass
class Construct:
    index: int
    original_name: str
    display_name: str
    insert_dna: str


def clean_dna(value: object) -> str:
    seq = str(value).upper()
    seq = re.sub(r"[^ACGTN]", "", seq)
    return seq


def safe_name(name: str) -> str:
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name.strip())
    name = re.sub(r"_+", "_", name).strip("_")
    return name or "construct"


def parse_name_map(value: str) -> Dict[str, str]:
    mapping = {}
    if not value:
        return mapping
    for item in value.split(","):
        if not item.strip():
            continue
        old, new = item.split(":", 1)
        mapping[old.strip()] = new.strip()
    return mapping


def apply_name_map(name: str, mapping: Dict[str, str]) -> str:
    for old, new in mapping.items():
        if name == old:
            return new
        if name.startswith(old + "_"):
            return new + name[len(old):]
    return name


def read_table(path: Path, name_col: str, dna_col: str, name_map: Dict[str, str]) -> List[Construct]:
    if path.suffix.lower() in {".xlsx", ".xlsm", ".xls"}:
        df = pd.read_excel(path)
    elif path.suffix.lower() == ".tsv":
        df = pd.read_csv(path, sep="\t")
    else:
        df = pd.read_csv(path)
    missing = [c for c in (name_col, dna_col) if c not in df.columns]
    if missing:
        raise ValueError(f"Missing column(s) in {path}: {missing}; available={list(df.columns)}")
    constructs = []
    for i, row in df.iterrows():
        original = str(row[name_col]).strip()
        display = apply_name_map(original, name_map)
        dna = clean_dna(row[dna_col])
        if not dna:
            raise ValueError(f"Empty DNA sequence at row {i + 2} ({original})")
        constructs.append(Construct(i + 1, original, display, dna))
    return constructs


def strand_to_int(value: object) -> Optional[int]:
    if value in {"+", 1, "+1"}:
        return 1
    if value in {"-", -1, "-1"}:
        return -1
    return None


def load_sequence_and_features(path: Path) -> Tuple[str, List[FeatureLite], str]:
    suffix = path.suffix.lower()
    if suffix == ".dna":
        if snapgene_file_to_dict is None:
            raise RuntimeError("Reading .dna requires snapgene_reader. Install with: pip install snapgene-reader")
        data = snapgene_file_to_dict(str(path))
        seq = data["seq"].upper()
        features = []
        for feat in data.get("features", []):
            qualifiers = dict(feat.get("qualifiers", {}) or {})
            label = feat.get("name") or qualifiers.get("label") or feat.get("type") or "feature"
            qualifiers.setdefault("label", label)
            features.append(
                FeatureLite(
                    int(feat["start"]),
                    int(feat["end"]),
                    strand_to_int(feat.get("strand")),
                    feat.get("type") or "misc_feature",
                    str(label),
                    qualifiers,
                )
            )
        topology = data.get("dna", {}).get("topology", "circular")
        return seq, features, topology
    if suffix in {".gb", ".gbk", ".genbank"}:
        rec = SeqIO.read(path, "genbank")
        features = []
        for feat in rec.features:
            start = int(feat.location.start)
            end = int(feat.location.end)
            qualifiers = dict(feat.qualifiers)
            label = first_qualifier(qualifiers, ["label", "gene", "product", "note"]) or feat.type
            features.append(FeatureLite(start, end, feat.location.strand, feat.type, label, qualifiers))
        topology = rec.annotations.get("topology", "circular")
        return str(rec.seq).upper(), features, topology
    rec = SeqIO.read(path, "fasta")
    return str(rec.seq).upper(), [], "linear"


def first_qualifier(qualifiers: Dict[str, object], keys: Iterable[str]) -> str:
    for key in keys:
        value = qualifiers.get(key)
        if value is None:
            continue
        if isinstance(value, list):
            return str(value[0]) if value else ""
        return str(value)
    return ""


def find_feature(features: List[FeatureLite], name: str) -> Optional[FeatureLite]:
    target = name.lower()
    for feat in features:
        haystack = " ".join(
            str(x)
            for x in [
                feat.label,
                feat.type,
                first_qualifier(feat.qualifiers, ["label", "gene", "product", "note"]),
            ]
        ).lower()
        if target in haystack:
            return feat
    return None


def locate_replacement(seq: str, args: argparse.Namespace) -> Tuple[int, int]:
    if args.replace_start and args.replace_end:
        start = args.replace_start - 1
        end = args.replace_end
        if start < 0 or end <= start or end > len(seq):
            raise ValueError("Invalid 1-based replacement coordinates")
        return start, end
    if not args.left_anchor or not args.right_anchor:
        raise ValueError("Provide --replace-start/--replace-end or --left-anchor/--right-anchor")
    left = clean_dna(args.left_anchor)
    right = clean_dna(args.right_anchor)
    start = seq.find(left)
    if start < 0:
        raise ValueError("Left anchor not found in vector")
    right_start = seq.find(right, start + len(left))
    if right_start < 0:
        raise ValueError("Right anchor not found downstream of left anchor")
    return start, right_start + len(right)


def shift_feature(feat: FeatureLite, replace_start: int, replace_end: int, delta: int) -> Optional[FeatureLite]:
    if feat.end <= replace_start:
        return feat
    if feat.start >= replace_end:
        return FeatureLite(feat.start + delta, feat.end + delta, feat.strand, feat.type, feat.label, feat.qualifiers)
    # Drop old features inside or overlapping the replaced region; add new insert features separately.
    return None


def to_seqfeature(feat: FeatureLite) -> SeqFeature:
    qualifiers = dict(feat.qualifiers)
    qualifiers.setdefault("label", feat.label)
    return SeqFeature(FeatureLocation(feat.start, feat.end, strand=feat.strand), type=feat.type, qualifiers=qualifiers)


def make_insert_features(c: Construct, insert_start: int, insert_len: int, args: argparse.Namespace) -> List[SeqFeature]:
    feats = [
        SeqFeature(
            FeatureLocation(insert_start, insert_start + insert_len, strand=1),
            type="misc_feature",
            qualifiers={"label": c.display_name, "note": f"Full synthesized Gibson insert from {c.original_name}."},
        )
    ]
    left_len = len(clean_dna(args.left_overlap_seq)) if args.left_overlap_seq else 0
    right_len = len(clean_dna(args.right_overlap_seq)) if args.right_overlap_seq else 0
    cds_start = insert_start + left_len
    cds_end = insert_start + insert_len - right_len
    if args.left_overlap_seq:
        feats.append(
            SeqFeature(
                FeatureLocation(insert_start, insert_start + left_len, strand=1),
                type="misc_feature",
                qualifiers={"label": "5p_Gibson_overlap"},
            )
        )
    if cds_end > cds_start and (cds_end - cds_start) % 3 == 0:
        cds = c.insert_dna[left_len : insert_len - right_len]
        feats.append(
            SeqFeature(
                FeatureLocation(cds_start, cds_end, strand=1),
                type="CDS",
                qualifiers={
                    "label": f"{c.display_name}_CDS",
                    "product": f"{c.display_name} insert CDS",
                    "translation": str(Seq(cds).translate(to_stop=False)),
                },
            )
        )
    if args.right_overlap_seq:
        feats.append(
            SeqFeature(
                FeatureLocation(insert_start + insert_len - right_len, insert_start + insert_len, strand=1),
                type="misc_feature",
                qualifiers={"label": "C_terminal_linker_right_overlap"},
            )
        )
    return feats


def write_genbank(seq: str, features: List[SeqFeature], path: Path, name: str, topology: str) -> None:
    rec = SeqRecord(Seq(seq), id=safe_name(name)[:16], name=safe_name(name)[:16], description=name)
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"] = topology
    rec.features = features
    path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(rec, path, "genbank")


def find_snapgene_cli(user_value: str) -> Optional[str]:
    if user_value:
        return user_value
    candidates = [
        shutil.which("SnapGene"),
        "/Applications/SnapGene.app/Contents/MacOS/SnapGene",
    ]
    return next((c for c in candidates if c and Path(c).exists()), None)


def convert_to_snapgene(gb_path: Path, dna_path: Path, snapgene_cli: Optional[str]) -> bool:
    if not snapgene_cli:
        return False
    dna_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [snapgene_cli, "-q", "-c", "SnapGene DNA", "-i", str(gb_path), "-o", str(dna_path)]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    # Some SnapGene builds write valid .dna files and then exit nonzero/segfault; verify by file readability later.
    return dna_path.exists() and dna_path.stat().st_size > 0


def read_sanger(path: Path) -> str:
    if path.suffix.lower() == ".ab1":
        return str(SeqIO.read(path, "abi").seq).upper()
    text = path.read_text(errors="ignore")
    lines = [line.strip() for line in text.splitlines() if line.strip() and not line.startswith(">")]
    return clean_dna("".join(lines))


def best_alignment(ref: str, query: str):
    candidates = []
    for orient, seq in [("forward", query), ("reverse_complement", str(Seq(query).reverse_complement()))]:
        aln = pairwise2.align.localms(ref, seq, 3, -5, -7, -2, one_alignment_only=True)[0]
        candidates.append((aln.score, orient, aln))
    candidates.sort(reverse=True, key=lambda x: x[0])
    return candidates[0][1], candidates[0][2]


def collect_diffs(ref: str, query: str):
    orient, aln = best_alignment(ref, query)
    ref_to_query = {}
    diffs = []
    ref_pos = 0
    for ca, cb in zip(aln.seqA, aln.seqB):
        if ca != "-":
            ref_pos += 1
            if cb != "-":
                ref_to_query[ref_pos] = cb.upper()
                if ca.upper() != cb.upper():
                    diffs.append((ref_pos, f"{ref_pos}:{ca}>{cb}", "sub"))
            else:
                diffs.append((ref_pos, f"{ref_pos}:del{ca}", "del"))
        else:
            if 0 < ref_pos < len(ref):
                diffs.append((ref_pos, f"{ref_pos}:ins{cb}", "ins"))
    return orient, aln.score, ref_to_query, diffs


def region_call(ref: str, ref_to_query: Dict[int, str], diffs: list, start0: int, end0: int, translate: bool):
    positions = range(start0 + 1, end0 + 1)
    covered = sum(1 for pos in positions if pos in ref_to_query)
    length = end0 - start0
    local_diffs = []
    for pos, diff, kind in diffs:
        if start0 < pos <= end0 or (kind == "ins" and start0 < pos < end0):
            local_diffs.append(f"{pos - start0}:{diff.split(':', 1)[1]}")
    dna_exact = covered == length and not local_diffs
    protein_exact = ""
    aa_diffs = ""
    if translate:
        expected = ref[start0:end0]
        observed = "".join(ref_to_query.get(pos, "N") for pos in positions)
        if covered == length and not any("ins" in d or "del" in d for d in local_diffs):
            exp_aa = str(Seq(expected).translate(to_stop=False))
            obs_aa = str(Seq(observed).translate(to_stop=False))
            protein_exact = exp_aa == obs_aa
            aa_diffs = ";".join(f"aa{i}:{a}>{b}" for i, (a, b) in enumerate(zip(exp_aa, obs_aa), 1) if a != b)
        else:
            protein_exact = False
            aa_diffs = "not fully covered or contains indel"
    return covered, length, len(local_diffs), dna_exact, protein_exact, ";".join(local_diffs), aa_diffs


def sanger_sample_id(path: Path) -> Optional[Tuple[int, int]]:
    m = re.match(r"^(\d+)-(\d+)\.", path.name)
    if not m:
        return None
    return int(m.group(1)), int(m.group(2))


def write_sanger_reports(constructs: List[Construct], second_seq: str, args: argparse.Namespace, out_dir: Path) -> None:
    if not args.sanger_success_dir:
        return
    success_dir = Path(args.sanger_success_dir)
    failed_dir = Path(args.sanger_failed_dir) if args.sanger_failed_dir else None
    if not success_dir.exists():
        raise FileNotFoundError(success_dir)

    left_len = len(clean_dna(args.left_overlap_seq)) if args.left_overlap_seq else 0
    right_len = len(clean_dna(args.right_overlap_seq)) if args.right_overlap_seq else 0
    by_index = {c.index: c for c in constructs}
    rows = []
    for path in sorted(list(success_dir.glob("*.ab1")) + list(success_dir.glob("*.seq"))):
        sid = sanger_sample_id(path)
        if not sid or path.suffix.lower() == ".seq" and (path.with_suffix(".ab1").exists()):
            continue
        ci, sample = sid
        c = by_index.get(ci)
        if c is None:
            continue
        ref = c.insert_dna + second_seq
        query = read_sanger(path)
        orient, score, ref_to_query, diffs = collect_diffs(ref, query)
        insert = region_call(ref, ref_to_query, diffs, 0, len(c.insert_dna), False)
        cds = region_call(ref, ref_to_query, diffs, left_len, len(c.insert_dna) - right_len, True)
        second = region_call(ref, ref_to_query, diffs, len(c.insert_dna), len(ref), bool(second_seq)) if second_seq else (0, 0, 0, "", "", "", "")
        if insert[3] and (not second_seq or second[3]):
            call = "INSERT_AND_SECOND_FRAGMENT_DNA_EXACT_PASS"
        elif cds[3] and (not second_seq or second[3]):
            call = "CDS_AND_SECOND_FRAGMENT_DNA_EXACT_PASS"
        elif cds[4] is True and (not second_seq or second[3]):
            call = "PROTEIN_AND_SECOND_FRAGMENT_DNA_EXACT_PASS"
        elif cds[4] is True and second_seq and second[4] is True:
            call = "PROTEIN_LEVEL_PASS_CONFIRM_DNA_IF_NEEDED"
        else:
            call = "FAIL_OR_INCOMPLETE"
        rows.append(
            {
                "construct_index": ci,
                "sample": sample,
                "file": path.name,
                "construct": c.display_name,
                "orientation": orient,
                "alignment_score": score,
                "insert_dna_exact": insert[3],
                "insert_differences": insert[5],
                "cds_dna_exact": cds[3],
                "cds_protein_exact": cds[4],
                "cds_differences": cds[5],
                "cds_aa_differences": cds[6],
                "second_fragment_dna_exact": second[3],
                "second_fragment_protein_exact": second[4],
                "second_fragment_differences": second[5],
                "second_fragment_aa_differences": second[6],
                "call": call,
            }
        )
    rows.sort(key=lambda r: (r["construct_index"], r["sample"]))
    per_read = out_dir / "sanger_results.csv"
    with per_read.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ["construct_index"])
        writer.writeheader()
        writer.writerows(rows)

    failed = defaultdict(list)
    if failed_dir and failed_dir.exists():
        for path in sorted(list(failed_dir.glob("*.ab1")) + list(failed_dir.glob("*.seq"))):
            sid = sanger_sample_id(path)
            if sid:
                failed[sid[0]].append(str(sid[1]))

    grouped = defaultdict(list)
    for row in rows:
        grouped[row["construct_index"]].append(row)
    summary = []
    rank = {
        "INSERT_AND_SECOND_FRAGMENT_DNA_EXACT_PASS": 4,
        "CDS_AND_SECOND_FRAGMENT_DNA_EXACT_PASS": 3,
        "PROTEIN_AND_SECOND_FRAGMENT_DNA_EXACT_PASS": 2,
        "PROTEIN_LEVEL_PASS_CONFIRM_DNA_IF_NEEDED": 1,
        "FAIL_OR_INCOMPLETE": 0,
    }
    for c in constructs:
        rs = sorted(grouped.get(c.index, []), key=lambda r: r["sample"])
        best = sorted(rs, key=lambda r: (-rank.get(r["call"], 0), r["sample"]))[0] if rs else None
        summary.append(
            {
                "construct_index": c.index,
                "construct": c.display_name,
                "success_samples": ",".join(str(r["sample"]) for r in rs),
                "insert_dna_exact_samples": ",".join(str(r["sample"]) for r in rs if r["insert_dna_exact"] is True),
                "second_fragment_dna_exact_samples": ",".join(str(r["sample"]) for r in rs if r["second_fragment_dna_exact"] is True),
                "both_dna_exact_samples": ",".join(str(r["sample"]) for r in rs if r["insert_dna_exact"] is True and (not second_seq or r["second_fragment_dna_exact"] is True)),
                "failed_folder_samples": ",".join(sorted(set(failed.get(c.index, [])), key=int)),
                "recommended_sample": best["sample"] if best and rank.get(best["call"], 0) > 0 else "",
                "final_call": best["call"] if best else "NO_SUCCESS_READ",
            }
        )
    summary_csv = out_dir / "sanger_summary.csv"
    with summary_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary[0].keys()) if summary else ["construct_index"])
        writer.writeheader()
        writer.writerows(summary)
    write_markdown_report(summary, out_dir / "sanger_report.md", args.second_feature_name)
    write_html_report(summary, out_dir / "sanger_report.html", args.second_feature_name)


def write_markdown_report(summary: List[dict], path: Path, second_name: str) -> None:
    with path.open("w") as f:
        f.write("# Sanger cloning validation report\n\n")
        if second_name:
            f.write(f"Second fragment reference: `{second_name}` from the input vector.\n\n")
        f.write("| # | Construct | Insert DNA exact | Second fragment DNA exact | Both DNA exact | Recommended | Final call |\n")
        f.write("|---:|---|---|---|---|---|---|\n")
        for r in summary:
            f.write(
                f"| {r['construct_index']} | {r['construct']} | {r['insert_dna_exact_samples']} | "
                f"{r['second_fragment_dna_exact_samples']} | {r['both_dna_exact_samples']} | "
                f"{r['recommended_sample']} | {r['final_call']} |\n"
            )


def badge(call: str) -> str:
    color = "#b91c1c" if "FAIL" in call or "NO_" in call else "#047857" if "DNA_EXACT" in call else "#b45309"
    return f'<span style="background:{color};color:white;border-radius:4px;padding:2px 6px;font-size:12px">{html.escape(call)}</span>'


def write_html_report(summary: List[dict], path: Path, second_name: str) -> None:
    rows = []
    for r in summary:
        rows.append(
            "<tr>"
            f"<td>{r['construct_index']}</td><td>{html.escape(r['construct'])}</td>"
            f"<td>{html.escape(str(r['insert_dna_exact_samples']))}</td>"
            f"<td>{html.escape(str(r['second_fragment_dna_exact_samples']))}</td>"
            f"<td>{html.escape(str(r['both_dna_exact_samples']))}</td>"
            f"<td>{html.escape(str(r['recommended_sample']))}</td><td>{badge(str(r['final_call']))}</td>"
            "</tr>"
        )
    path.write_text(
        "<!doctype html><meta charset='utf-8'><title>Sanger cloning validation</title>"
        "<style>body{font-family:Arial,sans-serif;margin:24px}table{border-collapse:collapse;width:100%}"
        "th,td{border:1px solid #ddd;padding:6px 8px;text-align:left}th{background:#f3f4f6}</style>"
        f"<h1>Sanger cloning validation</h1><p>Second fragment: {html.escape(second_name or 'none')}</p>"
        "<table><thead><tr><th>#</th><th>Construct</th><th>Insert DNA exact</th><th>Second fragment DNA exact</th>"
        "<th>Both DNA exact</th><th>Recommended</th><th>Final call</th></tr></thead><tbody>"
        + "".join(rows)
        + "</tbody></table>",
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate SnapGene maps and Sanger validation reports for Gibson assemblies.")
    parser.add_argument("--table", required=True, type=Path)
    parser.add_argument("--name-column", default="construct_id")
    parser.add_argument("--dna-column", default="DNA_final")
    parser.add_argument("--name-map", default="", help="Prefix map, e.g. A:targetA,B:targetB")
    parser.add_argument("--vector", required=True, type=Path)
    parser.add_argument("--left-anchor", default="")
    parser.add_argument("--right-anchor", default="")
    parser.add_argument("--replace-start", type=int, default=0, help="1-based inclusive")
    parser.add_argument("--replace-end", type=int, default=0, help="1-based inclusive")
    parser.add_argument("--left-overlap-seq", default="")
    parser.add_argument("--right-overlap-seq", default="")
    parser.add_argument("--second-feature-name", default="", help="Feature in vector to validate after the insert, e.g. 'fixed module'")
    parser.add_argument("--sanger-success-dir", default="")
    parser.add_argument("--sanger-failed-dir", default="")
    parser.add_argument("--snapgene-cli", default="")
    parser.add_argument("--skip-snapgene", action="store_true")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    out = args.out
    insert_gb = out / "01_synthesis_genbank"
    insert_dna = out / "02_synthesis_snapgene_dna"
    plasmid_gb = out / "03_cloned_plasmid_genbank"
    plasmid_dna = out / "04_cloned_plasmid_snapgene_dna"
    sanger_out = out / "05_sanger_alignment"
    for d in [insert_gb, insert_dna, plasmid_gb, plasmid_dna, sanger_out]:
        d.mkdir(parents=True, exist_ok=True)

    constructs = read_table(args.table, args.name_column, args.dna_column, parse_name_map(args.name_map))
    vector_seq, vector_features, topology = load_sequence_and_features(args.vector)
    replace_start, replace_end = locate_replacement(vector_seq, args)
    second_seq = ""
    if args.second_feature_name:
        feat = find_feature(vector_features, args.second_feature_name)
        if feat is None:
            raise ValueError(f"Second feature not found in vector: {args.second_feature_name}")
        second_seq = vector_seq[feat.start : feat.end]
        (sanger_out / "second_fragment_reference.fasta").write_text(
            f">{safe_name(args.second_feature_name)}_from_vector_1based_{feat.start + 1}_{feat.end}\n{second_seq}\n",
            encoding="utf-8",
        )

    snapgene_cli = None if args.skip_snapgene else find_snapgene_cli(args.snapgene_cli)
    validation_rows = []
    for c in constructs:
        prefix = f"{c.index:02d}_{safe_name(c.display_name)}"
        ins_gb_path = insert_gb / f"{prefix}.gb"
        write_genbank(c.insert_dna, make_insert_features(c, 0, len(c.insert_dna), args), ins_gb_path, c.display_name, "linear")
        if snapgene_cli:
            convert_to_snapgene(ins_gb_path, insert_dna / f"{prefix}.dna", snapgene_cli)
        cloned_seq = vector_seq[:replace_start] + c.insert_dna + vector_seq[replace_end:]
        delta = len(c.insert_dna) - (replace_end - replace_start)
        shifted = [to_seqfeature(f) for f in (shift_feature(feat, replace_start, replace_end, delta) for feat in vector_features) if f]
        shifted.extend(make_insert_features(c, replace_start, len(c.insert_dna), args))
        plasmid_name = f"{c.index:02d}_cloned_{safe_name(c.display_name)}"
        pl_gb_path = plasmid_gb / f"{plasmid_name}.gb"
        write_genbank(cloned_seq, shifted, pl_gb_path, plasmid_name, topology)
        if snapgene_cli:
            convert_to_snapgene(pl_gb_path, plasmid_dna / f"{plasmid_name}.dna", snapgene_cli)
        validation_rows.append({"construct_index": c.index, "construct": c.display_name, "insert_bp": len(c.insert_dna), "plasmid_bp": len(cloned_seq)})

    with (out / "generation_summary.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(validation_rows[0].keys()))
        writer.writeheader()
        writer.writerows(validation_rows)

    write_sanger_reports(constructs, second_seq, args, sanger_out)
    print(f"Generated outputs in {out}")
    print(f"Replacement region: 1-based {replace_start + 1}..{replace_end} ({replace_end - replace_start} bp)")
    if snapgene_cli:
        print(f"SnapGene CLI used: {snapgene_cli}")
    else:
        print("SnapGene .dna conversion skipped or SnapGene CLI not found; GenBank files were still generated.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
