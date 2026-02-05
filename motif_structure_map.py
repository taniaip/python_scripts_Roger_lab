#!/usr/bin/env python3

"""
Analyze the predicted RNA secondary-structure context of a motif across many sequences.

This script parses an RNAfold-style text file containing repeating 3-line blocks:
  >seq_id
  SEQUENCE
  DOTBRACKET (ENERGY)

For each sequence, it scans for all motif occurrences allowing up to a user-specified
number of mismatches (Hamming distance). For every motif hit, it summarizes the local
structural context over the motif window using the dot-bracket structure:

- Fraction of motif bases in stem (paired), loop (unpaired but enclosed), or external
  (unpaired and outside any pairing)
- Among paired motif positions, fraction that are “stacked” vs “loose” based on whether
  neighboring positions are also paired
- Start/end coordinates (0-based, end exclusive), matched motif sequence, mismatch count,
  and RNAfold free energy

Outputs:
  1) <outprefix>_motif_hits.csv
     All motif hits across all sequences with per-hit structural summaries.
  2) <outprefix>_best_hit_per_seq.csv
     One “best” hit per sequence (lowest mismatches; tie-break by higher loop_frac,
     then higher stem_frac, then more negative energy).

Also prints a brief summary of total hits and majority-context counts.

Dependencies:
  pandas

Author: Tania Iranpour
"""

import argparse
import re
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

import pandas as pd


@dataclass
class Record:
    seq_id: str
    seq: str
    dot: str
    energy: Optional[float]


def parse_rnafold_txt(path: str) -> List[Record]:
    """
    Parse a file with repeating blocks:
      >seq1
      UGA...
      (((...))) (-24.30)
    """
    records: List[Record] = []
    with open(path, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]

    i = 0
    while i < len(lines):
        if not lines[i].startswith(">"):
            raise ValueError(f"Expected header line starting with '>' at line {i+1}: {lines[i]}")
        seq_id = lines[i][1:].strip()
        if i + 2 >= len(lines):
            raise ValueError(f"Incomplete record for {seq_id}")

        seq = lines[i + 1].strip().upper()
        struct_line = lines[i + 2].strip()

        # dotbracket is first token, energy inside parentheses at end
        # Example: "(((...))) (-24.30)"
        m = re.match(r"^([().]+)\s+\(\s*([-+]?\d+(?:\.\d+)?)\s*\)\s*$", struct_line)
        if not m:
            raise ValueError(f"Could not parse structure+energy for {seq_id}: {struct_line}")

        dot = m.group(1)
        energy = float(m.group(2))

        if len(seq) != len(dot):
            raise ValueError(f"Length mismatch for {seq_id}: seq={len(seq)} dot={len(dot)}")

        records.append(Record(seq_id=seq_id, seq=seq, dot=dot, energy=energy))
        i += 3

    return records


def build_pair_map(dot: str) -> Dict[int, int]:
    """Return dict mapping i->j and j->i for paired positions (0-based)."""
    stack = []
    pairs: Dict[int, int] = {}
    for idx, ch in enumerate(dot):
        if ch == "(":
            stack.append(idx)
        elif ch == ")":
            if not stack:
                raise ValueError("Unbalanced dot-bracket: too many ')'")
            left = stack.pop()
            pairs[left] = idx
            pairs[idx] = left
    if stack:
        raise ValueError("Unbalanced dot-bracket: too many '('")
    return pairs


def structural_context_per_base(dot: str, pairs: Dict[int, int]) -> List[str]:
    """
    Label each position as:
      - 'stem' if paired
      - 'loop' if unpaired but enclosed inside any pair
      - 'external' if unpaired outside all pairs
    Enclosed check is done by tracking nesting depth.
    """
    ctx = []
    depth = 0
    for i, ch in enumerate(dot):
        if ch == "(":
            depth += 1
            ctx.append("stem")
        elif ch == ")":
            ctx.append("stem")
            depth -= 1
        else:  # '.'
            ctx.append("loop" if depth > 0 else "external")
    if depth != 0:
        raise ValueError("Dot-bracket depth not zero at end (unbalanced).")
    return ctx


def stacked_or_loose(dot: str, i: int) -> str:
    """
    Simple definition:
      - if position i is paired AND at least one neighbor is also paired -> 'stacked'
      - if paired but neighbors are unpaired/out of bounds -> 'loose'
      - if unpaired -> 'unpaired'
    """
    if dot[i] not in "()":
        return "unpaired"
    left_paired = (i - 1 >= 0 and dot[i - 1] in "()")
    right_paired = (i + 1 < len(dot) and dot[i + 1] in "()")
    return "stacked" if (left_paired or right_paired) else "loose"


def hamming_mismatches(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a, b) if x != y)


def find_motif_hits(seq: str, motif: str, max_mismatches: int) -> List[Tuple[int, int]]:
    """
    Return list of (start, mismatches) for all windows of length len(motif)
    with mismatches <= max_mismatches.
    """
    hits = []
    L = len(motif)
    for start in range(0, len(seq) - L + 1):
        window = seq[start:start + L]
        mm = hamming_mismatches(window, motif)
        if mm <= max_mismatches:
            hits.append((start, mm))
    return hits


def summarize_hit(seq: str, dot: str, ctx: List[str], start: int, motif: str) -> Dict:
    L = len(motif)
    end = start + L  # exclusive
    window_ctx = ctx[start:end]

    stem_ct = sum(1 for c in window_ctx if c == "stem")
    loop_ct = sum(1 for c in window_ctx if c == "loop")
    ext_ct = sum(1 for c in window_ctx if c == "external")

    # stacked vs loose among paired positions in motif window
    stacked_ct = 0
    loose_ct = 0
    for i in range(start, end):
        s = stacked_or_loose(dot, i)
        if s == "stacked":
            stacked_ct += 1
        elif s == "loose":
            loose_ct += 1

    paired_ct = stacked_ct + loose_ct
    stacked_frac = (stacked_ct / paired_ct) if paired_ct > 0 else 0.0
    loose_frac = (loose_ct / paired_ct) if paired_ct > 0 else 0.0

    return {
        "motif_seq": seq[start:end],
        "start_0based": start,
        "end_0based_excl": end,
        "stem_frac": stem_ct / L,
        "loop_frac": loop_ct / L,
        "external_frac": ext_ct / L,
        "paired_ct": paired_ct,
        "stacked_frac_of_paired": stacked_frac,
        "loose_frac_of_paired": loose_frac,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--infile", required=True, help="RNAfold-style txt with >id, seq, dotbracket (energy)")
    ap.add_argument("--motif", default="UGUUUGUU", help="Motif to scan (RNA alphabet)")
    ap.add_argument("--max_mismatches", type=int, default=1, help="Allow up to this many mismatches")
    ap.add_argument("--outprefix", default="motif_structure", help="Output prefix")
    args = ap.parse_args()

    motif = args.motif.upper()
    records = parse_rnafold_txt(args.infile)

    rows = []
    for r in records:
        pairs = build_pair_map(r.dot)
        ctx = structural_context_per_base(r.dot, pairs)

        hits = find_motif_hits(r.seq, motif, args.max_mismatches)
        for start, mm in hits:
            row = {
                "seq_id": r.seq_id,
                "energy": r.energy,
                "motif": motif,
                "mismatches": mm,
            }
            row.update(summarize_hit(r.seq, r.dot, ctx, start, motif))
            rows.append(row)

    df = pd.DataFrame(rows)
    out_all = f"{args.outprefix}_motif_hits.csv"
    df.to_csv(out_all, index=False)

    # Best hit per seq = lowest mismatches, then most stem/loop preference doesn’t matter,
    # so we tie-break by loop_frac then stem_frac then energy (more negative is stronger)
    if not df.empty:
        best = (
            df.sort_values(
                ["seq_id", "mismatches", "loop_frac", "stem_frac", "energy"],
                ascending=[True, True, False, False, True],
            )
            .groupby("seq_id", as_index=False)
            .head(1)
        )
    else:
        best = df

    out_best = f"{args.outprefix}_best_hit_per_seq.csv"
    best.to_csv(out_best, index=False)

    # Quick summary
    print(f"Parsed sequences: {len(records)}")
    print(f"Total motif hits (<= {args.max_mismatches} mismatches): {len(df)}")
    if len(df) > 0:
        # simple categorical “majority context”
        def majority_context(row):
            vals = {"stem": row["stem_frac"], "loop": row["loop_frac"], "external": row["external_frac"]}
            return max(vals, key=vals.get)

        df2 = df.copy()
        df2["majority_context"] = df2.apply(majority_context, axis=1)
        print("\nHit majority context counts:")
        print(df2["majority_context"].value_counts().to_string())

        print("\nBest-hit majority context counts (per seq):")
        if len(best) > 0:
            best2 = best.copy()
            best2["majority_context"] = best2.apply(majority_context, axis=1)
            print(best2["majority_context"].value_counts().to_string())

    print(f"\nWrote:\n  {out_all}\n  {out_best}")


if __name__ == "__main__":
    main()

