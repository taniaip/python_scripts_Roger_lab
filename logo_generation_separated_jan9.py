#!/usr/bin/env python3

#FOR THE NUMBERS AND COMMANDS FOR EACH OF THE CATEGORIES GO TO THE END OF THIS SCRIPT.

"""
Extract polyA-centered sequence windows, classify each window by the best-fitting
motif placement pattern, and generate sequence logos for multiple categories.

This script reads poly(A) sites from a GFF3 file (feature type 'polyA'), extracts a
fixed window around each site from a reference genome FASTA (default flank=24, so
window length = 2*flank+1), reverse-complements minus-strand sites, and then assigns
each window to motif-based categories using fuzzy matching against the motif
TGTTTGTT (up to 2 substitutions; no indels).

Two complementary classification schemes are produced:

1) Best-class (mutually exclusive):
   Each sequence is matched against a set of candidate placements that represent a
   leading base (e.g., 'T', 'TG', or 'TA'), followed by a gap of size k (1–8), followed
   by the fuzzy motif. If multiple candidates match a sequence, the single best label
   is chosen by:
     (i)  fewest substitutions vs motif,
     (ii) smallest gap size,
     (iii) deterministic label priority.

   Best-class labels include: gap1..gap8 and special 4-gap variants (4a, 4b, 4c).

2) pattern_k_all (first-match, not best-of-all):
   For k = 1..8, sequences are placed into the first k whose pattern matches:
     ^(.{24})([ACGT]{1})([ACGT]{k})(TGTTTGTT){s<=2}

Outputs:
- Sequence logo PNGs for:
    • all extracted sequences
    • the ELSE set (no candidate matched)
    • each best-class bucket (gap1..gap8, 4a/4b/4c)
    • each pattern_k_all bucket (k=1..8)
- A GFF3 file containing poly(A) entries whose extracted windows fell into ELSE.

Inputs:
  -g1 / --gff3_in   Input GFF3 of polyA sites
  -g2 / --gff3_else Output GFF3 for ELSE polyA sites
  -f  / --fasta     Reference genome FASTA
  --flank           Flank size around polyA site (default: 24)
  --outAll, --outElse, --out1..--out8, --out4a/--out4b/--out4c,
  --out1_all..--out8_all
                    Output logo basenames (PNG extension added)

Dependencies:
  Biopython, regex, matplotlib, logomaker

Author: Tania Iranpour
"""

import argparse
import regex
from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import logomaker


MOTIF = "TGTTTGTT"
MAX_SUBS = 2


def parse_args():
    p = argparse.ArgumentParser(
        description="Extract polyA-centered sequences, classify by best motif placement, and generate logos."
    )
    p.add_argument("-g1", "--gff3_in", dest="gff3_in", required=True,
                   help="Input GFF3 file with polyA site info (feature type must be 'polyA').")
    p.add_argument("-g2", "--gff3_else", dest="gff3_else", required=True,
                   help="Output GFF3 for sequences that do not match any candidate (else).")
    p.add_argument("-f", "--fasta", dest="fasta_file", required=True,
                   help="Input FASTA file of genomic sequences (seqid must match GFF3 column 1).")
    p.add_argument("--flank", type=int, default=24,
                   help="Number of bases upstream/downstream to extract around polyA site (default: 24).")

    # Logos for all sequences and else
    p.add_argument("--outAll",  dest="out_all",  required=True, help="Logo output prefix for ALL sequences.")
    p.add_argument("--outElse", dest="out_else", required=True, help="Logo output prefix for ELSE sequences.")

    # Best-class logos (chosen by least substitutions, then smallest gap)
    p.add_argument("--out1",  dest="out_1",  required=True, help="Logo output prefix for best-class gap1.")
    p.add_argument("--out2",  dest="out_2",  required=True, help="Logo output prefix for best-class gap2.")
    p.add_argument("--out3",  dest="out_3",  required=True, help="Logo output prefix for best-class gap3.")
    p.add_argument("--out4",  dest="out_4",  required=True, help="Logo output prefix for best-class gap4.")
    p.add_argument("--out5",  dest="out_5",  required=True, help="Logo output prefix for best-class gap5.")
    p.add_argument("--out6",  dest="out_6",  required=True, help="Logo output prefix for best-class gap6.")
    p.add_argument("--out7",  dest="out_7",  required=True, help="Logo output prefix for best-class gap7.")
    p.add_argument("--out8",  dest="out_8",  required=True, help="Logo output prefix for best-class gap8.")
    p.add_argument("--out4a", dest="out_4a", required=True, help="Logo output prefix for best-class 4a (just T).")
    p.add_argument("--out4b", dest="out_4b", required=True, help="Logo output prefix for best-class 4b (TG).")
    p.add_argument("--out4c", dest="out_4c", required=True, help="Logo output prefix for best-class 4c (TA).")

    # "k_all" logos: pattern_k_all for k = 1..8
    p.add_argument("--out1_all", dest="out_1_all", required=True, help="Logo output prefix for pattern_1_all.")
    p.add_argument("--out2_all", dest="out_2_all", required=True, help="Logo output prefix for pattern_2_all.")
    p.add_argument("--out3_all", dest="out_3_all", required=True, help="Logo output prefix for pattern_3_all.")
    p.add_argument("--out4_all", dest="out_4_all", required=True, help="Logo output prefix for pattern_4_all.")
    p.add_argument("--out5_all", dest="out_5_all", required=True, help="Logo output prefix for pattern_5_all.")
    p.add_argument("--out6_all", dest="out_6_all", required=True, help="Logo output prefix for pattern_6_all.")
    p.add_argument("--out7_all", dest="out_7_all", required=True, help="Logo output prefix for pattern_7_all.")
    p.add_argument("--out8_all", dest="out_8_all", required=True, help="Logo output prefix for pattern_8_all.")

    return p.parse_args()


def build_candidates():
    """
    Candidate placements that compete for the ONE best classification per sequence.
    Ranking rule:
      1) fewest substitutions vs motif
      2) if tie, smallest gap
      3) if tie, deterministic label order

    Returns list of tuples: (label, gap_size, compiled_regex)
    """
    motif_fuzzy = rf"({MOTIF}){{s<={MAX_SUBS}}}"
    candidates = []

    # gap1..8: ^(.{24})(T)([ACGT]{k})(TGTTTGTT){s<=2}
    for k in range(1, 9):
        candidates.append((
            f"gap{k}",
            k,
            regex.compile(rf"^(.{{24}})(T)([ACGT]{{{k}}}){motif_fuzzy}")
        ))

    # Special 4-gap variants
    candidates.extend([
        ("4a", 4, regex.compile(rf"^(.{{24}})(T)([ACGT]{{4}}){motif_fuzzy}")),
        ("4b", 4, regex.compile(rf"^(.{{23}})(TG)([ACGT]{{4}}){motif_fuzzy}")),
        ("4c", 4, regex.compile(rf"^(.{{23}})(TA)([ACGT]{{4}}){motif_fuzzy}")),
    ])

    return candidates


def build_pattern_k_all():
    """
    pattern_k_all[k] matches:
      ^(.{24})([ACGT]{1})([ACGT]{k})(TGTTTGTT){s<=2}
    for k = 1..8
    """
    patterns = {}
    motif_fuzzy = rf"({MOTIF}){{s<={MAX_SUBS}}}"
    for k in range(1, 9):
        patterns[k] = regex.compile(
            rf"^(.{{24}})([ACGT]{{1}})([ACGT]{{{k}}}){motif_fuzzy}"
        )
    return patterns


def best_motif_class(region_seq, candidates):
    """
    Return best label based on:
      - minimum substitutions (m.fuzzy_counts[0])
      - then smaller gap
      - then deterministic label ranking

    Returns: (label, match_obj, subs) or (None, None, None)
    """
    # deterministic tie-break among equal subs and gap
    label_rank = {
        "gap1": 1,
        "gap2": 2,
        "gap3": 3,
        "4a": 4,
        "4b": 5,
        "4c": 6,
        "gap4": 7,
        "gap5": 8,
        "gap6": 9,
        "gap7": 10,
        "gap8": 11,
    }

    hits = []
    for label, gap_size, pat in candidates:
        m = pat.match(region_seq)
        if not m:
            continue
        subs, ins, dels = m.fuzzy_counts  # substitutions, insertions, deletions
        hits.append((subs, gap_size, label_rank.get(label, 999), label, m))

    if not hits:
        return None, None, None

    hits.sort(key=lambda x: (x[0], x[1], x[2]))
    subs, gap_size, _, label, m = hits[0]
    return label, m, subs


def extract_polyA_sequences(gff3_in, fasta_file, flank=24):
    fasta_seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(fasta_file, "fasta")}

    candidates = build_candidates()
    pattern_k_all = build_pattern_k_all()

    sequences_all = []
    sequences_else = []
    else_gff_lines = []

    # Best-class buckets
    best = {
        "gap1": [], "gap2": [], "gap3": [], "gap4": [], "gap5": [], "gap6": [], "gap7": [], "gap8": [],
        "4a": [], "4b": [], "4c": []
    }

    # pattern_k_all buckets (k=1..8)
    sequences_k_all = {k: [] for k in range(1, 9)}

    with open(gff3_in, "r") as gf:
        for raw in gf:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 9:
                continue
            if parts[2] != "polyA":
                continue

            seqid = parts[0]
            if seqid not in fasta_seqs:
                continue

            start_0based = int(parts[3]) - 1
            strand = parts[6]

            full_seq = fasta_seqs[seqid]
            region_seq = full_seq[max(0, start_0based - flank): start_0based + flank + 1].upper()

            if strand == "-":
                region_seq = str(Seq(region_seq).reverse_complement())

            sequences_all.append(region_seq)

            # Fill pattern_k_all (choose smallest k if multiple match)
            for k in range(1, 9):
                if pattern_k_all[k].match(region_seq):
                    sequences_k_all[k].append(region_seq)
                    break

            # Best-class classification across all candidate placements
            label, m, subs = best_motif_class(region_seq, candidates)
            if label is None:
                sequences_else.append(region_seq)
                else_gff_lines.append(line)
            else:
                best[label].append(region_seq)

    # Print summary
    print(f"Total sequences: {len(sequences_all)}")
    for k in range(1, 9):
        print(f"pattern_{k}_all matches: {len(sequences_k_all[k])}")

    print("Best-class matches:")
    for k in range(1, 9):
        print(f"  gap{k}: {len(best[f'gap{k}'])}")
    print(f"  4a: {len(best['4a'])}")
    print(f"  4b: {len(best['4b'])}")
    print(f"  4c: {len(best['4c'])}")
    print(f"else (no match): {len(sequences_else)}")

    return sequences_all, sequences_else, else_gff_lines, best, sequences_k_all


def save_logo_logomaker(sequences, output_prefix):
    if not sequences:
        print(f"No sequences for {output_prefix}, skipping logo.")
        return

    info_df = logomaker.alignment_to_matrix(sequences, to_type="information")

    plt.figure(figsize=(max(2, len(info_df) / 3), 3))
    logo = logomaker.Logo(info_df, color_scheme="classic")
    logo.style_spines(visible=False)
    logo.style_spines(spines=["left", "bottom"], visible=True)
    logo.ax.set_ylabel("bits")
    logo.ax.set_ylim([0, 2.0])
    plt.title(f"Sequence Logo: {output_prefix}")
    plt.tight_layout()
    plt.savefig(output_prefix + ".png", dpi=200)
    plt.close()


def main():
    args = parse_args()

    seq_all, seq_else, else_gff_lines, best, sequences_k_all = extract_polyA_sequences(
        gff3_in=args.gff3_in,
        fasta_file=args.fasta_file,
        flank=args.flank
    )

    # Logos: overall + else
    save_logo_logomaker(seq_all, args.out_all)
    save_logo_logomaker(seq_else, args.out_else)

    # Logos: best-class
    save_logo_logomaker(best["gap1"], args.out_1)
    save_logo_logomaker(best["gap2"], args.out_2)
    save_logo_logomaker(best["gap3"], args.out_3)
    save_logo_logomaker(best["gap4"], args.out_4)
    save_logo_logomaker(best["gap5"], args.out_5)
    save_logo_logomaker(best["gap6"], args.out_6)
    save_logo_logomaker(best["gap7"], args.out_7)
    save_logo_logomaker(best["gap8"], args.out_8)
    save_logo_logomaker(best["4a"], args.out_4a)
    save_logo_logomaker(best["4b"], args.out_4b)
    save_logo_logomaker(best["4c"], args.out_4c)

    print("\nBest-class polyA counts:")
    for key in ["gap1","gap2","gap3","gap4","gap5","gap6","gap7","gap8","4a","4b","4c"]:
        print(f"  {key}: {len(best[key])}")

    # Logos: pattern_k_all (k=1..8)
    save_logo_logomaker(sequences_k_all[1], args.out_1_all)
    save_logo_logomaker(sequences_k_all[2], args.out_2_all)
    save_logo_logomaker(sequences_k_all[3], args.out_3_all)
    save_logo_logomaker(sequences_k_all[4], args.out_4_all)
    save_logo_logomaker(sequences_k_all[5], args.out_5_all)
    save_logo_logomaker(sequences_k_all[6], args.out_6_all)
    save_logo_logomaker(sequences_k_all[7], args.out_7_all)
    save_logo_logomaker(sequences_k_all[8], args.out_8_all)

    # Write ELSE GFF3
    with open(args.gff3_else, "w") as f:
        for gff_line in else_gff_lines:
            f.write(gff_line + "\n")


if __name__ == "__main__":
    main()


#FOR NON-STOP ST2:
# $ python ../../logo_generation_separated_jan9.py \
#   -g1 ../filtered_non_stop_polyA.gff3 \
#   -g2 else_polyAs_has_stop_ST2_logo.gff3 \
#   -f ../ST2_sorted_masked.fasta \
#   --outAll ST2_logo_ALL_NON \
#   --outElse ST2_logo_else \
#   --out1 ST2_logo_gap1 \
#   --out2 ST2_logo_gap2 \
#   --out3 ST2_logo_gap3 \
#   --out4 ST2_logo_gap4 \
#   --out5 ST2_logo_gap5 \
#   --out6 ST2_logo_gap6 \
#   --out7 ST2_logo_gap7 \
#   --out8 ST2_logo_gap8 \
#   --out4a ST2_logo_4a \
#   --out4b ST2_logo_4b \
#   --out4c ST2_logo_4c \
#   --out1_all ST2_logo_1_all \
#   --out2_all ST2_logo_2_all \
#   --out3_all ST2_logo_3_all \
#   --out4_all ST2_logo_4_all \
#   --out5_all ST2_logo_5_all \
#   --out6_all ST2_logo_6_all \
#   --out7_all ST2_logo_7_all \
#   --out8_all ST2_logo_8_all

# Total sequences: 938
# pattern_1_all matches: 8
# pattern_2_all matches: 4
# pattern_3_all matches: 35
# pattern_4_all matches: 843
# pattern_5_all matches: 11
# pattern_6_all matches: 4
# pattern_7_all matches: 0
# pattern_8_all matches: 1

# Best-class matches:
#   gap1: 0
#   gap2: 3
#   gap3: 1
#   gap4: 0
#   gap5: 4
#   gap6: 0
#   gap7: 0
#   gap8: 1
#   4a: 668
#   4b: 86
#   4c: 88
# else (no match): 87



#FOR HAS-STOP ST2:
# python ../../logo_generation_separated_jan9.py \
#   -g1 ../filtered_has_stop_polyA.gff3 \ 
#   -g2 else_polyAs_has_stop_ST2_logo.gff3 \
#   -f ../ST2_sorted_masked.fasta \
#   --outAll ST2_logo_ALL_HAS \
#   --outElse ST2_logo_else \
#   --out1 ST2_logo_gap1 \
#   --out2 ST2_logo_gap2 \
#   --out3 ST2_logo_gap3 \
#   --out4 ST2_logo_gap4 \
#   --out5 ST2_logo_gap5 \
#   --out6 ST2_logo_gap6 \
#   --out7 ST2_logo_gap7 \
#   --out8 ST2_logo_gap8 \
#   --out4a ST2_logo_4a \
#   --out4b ST2_logo_4b \
#   --out4c ST2_logo_4c \
#   --out1_all ST2_logo_1_all \
#   --out2_all ST2_logo_2_all \
#   --out3_all ST2_logo_3_all \
#   --out4_all ST2_logo_4_all \
#   --out5_all ST2_logo_5_all \
#   --out6_all ST2_logo_6_all \
#   --out7_all ST2_logo_7_all \
#   --out8_all ST2_logo_8_all

# Total sequences: 4277
# pattern_1_all matches: 126
# pattern_2_all matches: 77
# pattern_3_all matches: 340
# pattern_4_all matches: 3368
# pattern_5_all matches: 104
# pattern_6_all matches: 24
# pattern_7_all matches: 16
# pattern_8_all matches: 6

# Best-class matches:
#   gap1: 8
#   gap2: 14
#   gap3: 11
#   gap4: 0
#   gap5: 46
#   gap6: 1
#   gap7: 2
#   gap8: 3
#   4a: 1776
#   4b: 356
#   4c: 209
# else (no match): 1851



# FOR ALL POLYS IN ST2:
# python logo_generation_separated_jan9.py \
#   -g1 ST2_filtered_all_polyAs.gff3 \
#   -g2 else_polyAs_ALL_ST2_logo.gff3 \
#   -f ST2_sorted_masked.fasta \
#   --outAll ST2_logo_ALL \
#   --outElse ST2_logo_else \
#   --out1 ST2_logo_gap1 \
#   --out2 ST2_logo_gap2 \
#   --out3 ST2_logo_gap3 \
#   --out4 ST2_logo_gap4 \
#   --out5 ST2_logo_gap5 \
#   --out6 ST2_logo_gap6 \
#   --out7 ST2_logo_gap7 \
#   --out8 ST2_logo_gap8 \
#   --out4a ST2_logo_4a \
#   --out4b ST2_logo_4b \
#   --out4c ST2_logo_4c \
#   --out1_all ST2_logo_1_all \
#   --out2_all ST2_logo_2_all \
#   --out3_all ST2_logo_3_all \
#   --out4_all ST2_logo_4_all \
#   --out5_all ST2_logo_5_all \
#   --out6_all ST2_logo_6_all \
#   --out7_all ST2_logo_7_all \
#   --out8_all ST2_logo_8_all

# Total sequences: 5597
# pattern_1_all matches: 143
# pattern_2_all matches: 88
# pattern_3_all matches: 400
# pattern_4_all matches: 4435
# pattern_5_all matches: 131
# pattern_6_all matches: 31
# pattern_7_all matches: 19
# pattern_8_all matches: 10
#else: 340

# Best-class matches:
#   gap1: 9
#   gap2: 19
#   gap3: 13
#   gap4: 0
#   gap5: 60
#   gap6: 1
#   gap7: 2
#   gap8: 4
#   4a: 2589
#   4b: 461
#   4c: 310
# else (no match): 2129



#FOR ST7B:
#FOR ALL CATEGORIES:
# python ../../logo_generation_separated_jan9.py \
#   -g1 ../filtered_ST7B_all_polyA.gff3 \
#   -g2 else_polyAs_ALL_ST7B_logo.gff3 \
#   -f ../Blastocystis_ST7B_Complete_Assembly_pilon_revised.fasta \
#   --outAll ST7B_logo_ALL \
#   --outElse ST7B_logo_else \
#   --out1 ST2_logo_gap1 \
#   --out2 ST2_logo_gap2 \
#   --out3 ST2_logo_gap3 \
#   --out4 ST2_logo_gap4 \
#   --out5 ST2_logo_gap5 \
#   --out6 ST2_logo_gap6 \
#   --out7 ST2_logo_gap7 \
#   --out8 ST2_logo_gap8 \
#   --out4a ST2_logo_4a \
#   --out4b ST2_logo_4b \
#   --out4c ST2_logo_4c \
#   --out1_all ST7B_logo_1_all \
#   --out2_all ST7B_logo_2_all \
#   --out3_all ST7B_logo_3_all \
#   --out4_all ST7B_logo_4_all \
#   --out5_all ST7B_logo_5_all \
#   --out6_all ST7B_logo_6_all \
#   --out7_all ST7B_logo_7_all \
#   --out8_all ST7B_logo_8_all
# Matplotlib is building the font cache; this may take a moment.
# Total sequences: 5103
# pattern_1_all matches: 257
# pattern_2_all matches: 151
# pattern_3_all matches: 452
# pattern_4_all matches: 3640
# pattern_5_all matches: 111
# pattern_6_all matches: 28
# pattern_7_all matches: 16
# pattern_8_all matches: 15
# Best-class matches:
#   gap1: 26
#   gap2: 7
#   gap3: 14
#   gap4: 0
#   gap5: 62
#   gap6: 8
#   gap7: 8
#   gap8: 7
#   4a: 2630
#   4b: 263
#   4c: 196
# else (no match): 1882

# #FOR NON0STOP:
# python ../../logo_generation_separated_jan9.py \
#   -g1 ../filtered_ST7B_non_stop_polyA.gff3 \
#   -g2 else_polyAs_non_stop_ST7B_logo.gff3 \
#   -f ../Blastocystis_ST7B_Complete_Assembly_pilon_revised.fasta \
#   --outAll ST7B_logo_ALL_NON \
#   --outElse ST7B_logo_else \
#   --out1 ST2_logo_gap1 \
#   --out2 ST2_logo_gap2 \
#   --out3 ST2_logo_gap3 \
#   --out4 ST2_logo_gap4 \
#   --out5 ST2_logo_gap5 \
#   --out6 ST2_logo_gap6 \
#   --out7 ST2_logo_gap7 \
#   --out8 ST2_logo_gap8 \
#   --out4a ST2_logo_4a \
#   --out4b ST2_logo_4b \
#   --out4c ST2_logo_4c \
#   --out1_all ST7B_logo_1_all \
#   --out2_all ST7B_logo_2_all \
#   --out3_all ST7B_logo_3_all \
#   --out4_all ST7B_logo_4_all \
#   --out5_all ST7B_logo_5_all \
#   --out6_all ST7B_logo_6_all \
#   --out7_all ST7B_logo_7_all \
#   --out8_all ST7B_logo_8_all
# Matplotlib is building the font cache; this may take a moment.
# Total sequences: 669
# pattern_1_all matches: 17
# pattern_2_all matches: 0
# pattern_3_all matches: 29
# pattern_4_all matches: 598
# pattern_5_all matches: 8
# pattern_6_all matches: 0
# pattern_7_all matches: 0
# pattern_8_all matches: 0
# Best-class matches:
#   gap1: 1
#   gap2: 0
#   gap3: 0
#   gap4: 0
#   gap5: 4
#   gap6: 0
#   gap7: 0
#   gap8: 0
#   4a: 573
#   4b: 11
#   4c: 24
# else (no match): 56

# #FOR HAS STOP:
# $ python ../../logo_generation_separated_jan9.py \
#   -g1 ../filtered_ST7B_has_stop_polyA.gff3 \
#   -g2 else_polyAs_has_stop_ST7B_logo.gff3 \
#   -f ../Blastocystis_ST7B_Complete_Assembly_pilon_revised.fasta \
#   --outAll ST7B_logo_ALL_HAS \
#   --outElse ST7B_logo_else \
#   --out1 ST2_logo_gap1 \
#   --out2 ST2_logo_gap2 \
#   --out3 ST2_logo_gap3 \
#   --out4 ST2_logo_gap4 \
#   --out5 ST2_logo_gap5 \
#   --out6 ST2_logo_gap6 \
#   --out7 ST2_logo_gap7 \
#   --out8 ST2_logo_gap8 \
#   --out4a ST2_logo_4a \
#   --out4b ST2_logo_4b \
#   --out4c ST2_logo_4c \
#   --out1_all ST7B_logo_1_all \
#   --out2_all ST7B_logo_2_all \
#   --out3_all ST7B_logo_3_all \
#   --out4_all ST7B_logo_4_all \
#   --out5_all ST7B_logo_5_all \
#   --out6_all ST7B_logo_6_all \
#   --out7_all ST7B_logo_7_all \
#   --out8_all ST7B_logo_8_all
# Matplotlib is building the font cache; this may take a moment.
# Total sequences: 3610
# pattern_1_all matches: 201
# pattern_2_all matches: 125
# pattern_3_all matches: 367
# pattern_4_all matches: 2624
# pattern_5_all matches: 87
# pattern_6_all matches: 21
# pattern_7_all matches: 11
# pattern_8_all matches: 9
# Best-class matches:
#   gap1: 19
#   gap2: 6
#   gap3: 10
#   gap4: 0
#   gap5: 50
#   gap6: 8
#   gap7: 6
#   gap8: 7
#   4a: 1742
#   4b: 228
#   4c: 152
# else (no match): 1382




########################################################################################################

####THIS IS FULLY FUNCTIONAL, TO ACCOUNT FOR OTHER GAP VALUES(1,2,6,7,8), THE MODIFY VERSION IS USED.

# # python logo_generation_separated_jan9.py -g1 filtered_non_stop_polyA.gff3 -g2 else_polyAs_non_stop_logo.gff3 -f ST2_sorted_masked.fasta --outAll filtered_non_stop_all --out3 filtered_non_stop_gap3 --out4a filtered_non_stop_gap4a --out4b filtered_non_stop_gap4b --out4c filtered_non_stop_gap4c --out5 filtered_non_stop_gap5 --outElse filtered_non_stop_else --out3_all filtered_non_stop_gap3_all --out4_all filtered_non_stop_gap4_all --out5_all filtered_non_stop_gap5_all 

# # Matplotlib is building the font cache; this may take a moment.
# # Total sequences: 938
# # 3-gap all matches: 35
# # 3-gap matches: 1
# # 4-gap all matches: 853
# # 4-gap matches, just T: 671
# # 4-gap matches, TG: 86
# # 4-gap matches, TA: 88
# # 5-gap all matches: 13
# # 5-gap matches: 4
# # else (no match): 88

# # python logo_generation_separated_jan9.py -g1 march4_has_stop_polyA.gff3 -g2 else_polyAs_has_stop_logo.gff3 -f ST2_sorted_masked.fasta --outAll has_stop_all --out3 has_stop_gap3 --out4a has_stop_gap4a --out4b has_stop_gap4b --out4c has_stop_gap4c --out5 has_stop_gap5 --outElse has_stop_else --out3_all has_stop_gap3_all --out4_all has_stop_gap4_all --out5_all has_stop_gap5_all

# # Total sequences: 5163
# # 3-gap all matches: 516
# # 3-gap matches: 52
# # 4-gap all matches: 3582
# # 4-gap matches, just T: 1828
# # 4-gap matches, TG: 372
# # 4-gap matches, TA: 215
# # 5-gap all matches: 501
# # 5-gap matches: 110
# # else (no match): 2586

# # Total sequences: 5597
# # 3-gap all matches: 403
# # 3-gap matches: 13
# # 4-gap all matches: 4564
# # 4-gap matches, just T: 2606
# # 4-gap matches, TG: 461
# # 4-gap matches, TA: 310
# # 5-gap all matches: 155
# # 5-gap matches: 63
# # else (no match): 2144

# import argparse
# from Bio import SeqIO
# from Bio.Seq import Seq
# import regex
# import pandas as pd
# import logomaker
# import matplotlib.pyplot as plt
# import io

# def parse_args():
#     parser = argparse.ArgumentParser(description="Generate sequence logos (PNG) and write 'else' GFF3 lines.")
#     parser.add_argument("-g1", "--gff3_in", dest="gff3_in", type=str, required=True,
#                         help="Input GFF3 file with polyA site info.")
#     parser.add_argument("-g2", "--gff3_else", dest="gff3_else", type=str, required=True,
#                         help="Output GFF3 for 'else' category.")
#     parser.add_argument("-f", "--fasta", dest="fasta_file", type=str, required=True,
#                         help="Input FASTA file of genomic sequences.")
#     parser.add_argument("--outAll",  dest="out_all",  required=True, help="Logo for all sequences.")
#     parser.add_argument("--out3",    dest="out_3",    required=True, help="Logo for 3-gap sequences.")
#     parser.add_argument("--out4a",   dest="out_4a",   required=True, help="Logo for 4-gap (just T).")
#     parser.add_argument("--out4b",   dest="out_4b",   required=True, help="Logo for 4-gap (TG).")
#     parser.add_argument("--out4c",   dest="out_4c",   required=True, help="Logo for 4-gap (TA).")
#     parser.add_argument("--out5",    dest="out_5",    required=True, help="Logo for 5-gap sequences.")
#     parser.add_argument("--outElse", dest="out_else", required=True, help="Logo for 'else' sequences.")
#     parser.add_argument("--out3_all",    dest="out_3_all",    required=True, help="Logo for 3-gap sequences.")
#     parser.add_argument("--out4_all",    dest="out_4_all",    required=True, help="Logo for 4-gap sequences.")
#     parser.add_argument("--out5_all",    dest="out_5_all",    required=True, help="Logo for 5-gap sequences.")

#     return parser.parse_args()

# def extract_polyA_sequences(gff3_in, fasta_file, flank=24):
#     sequences_all  = []
#     sequences_3    = []
#     sequences_4a   = []
#     sequences_4b   = []
#     sequences_4c   = []
#     sequences_5    = []
#     sequences_else = []
#     else_gff_lines = []
#     sequences_3_all= []
#     sequences_4_all= []
#     sequences_5_all= []


#     pattern_3  = regex.compile(r"^(.{24})(T)([ACGT]{3})(TGTTTGTT){s<=2}")

#     pattern_4a = regex.compile(r"^(.{24})(T)([ACGT]{4})(TGTTTGTT){s<=2}")
#     pattern_4b = regex.compile(r"^(.{23})(TG)([ACGT]{4})(TGTTTGTT){s<=2}")
#     pattern_4c = regex.compile(r"^(.{23})(TA)([ACGT]{4})(TGTTTGTT){s<=2}")

#     pattern_5  = regex.compile(r"^(.{24})(T)([[ACGT]{5})(TGTTTGTT){s<=2}")

#     pattern_3_all = regex.compile(r"^(.{24})([ACGT]{1})([ACGT]{3})(TGTTTGTT){s<=2}")
#     pattern_4_all = regex.compile(r"^(.{24})([ACGT]{1})([ACGT]{4})(TGTTTGTT){s<=2}")
#     pattern_5_all = regex.compile(r"^(.{24})([ACGT]{1})([ACGT]{5})(TGTTTGTT){s<=2}")

#     fasta_seqs = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}

#     with open(gff3_in, "r") as gf:
#         for line in gf:
#             line = line.strip()
#             if line.startswith("#"):
#                 continue
#             parts = line.split("\t")
#             if len(parts) < 9 or parts[2] != "polyA":
#                 continue

#             seqid = parts[0]
#             start = int(parts[3]) - 1
#             strand = parts[6]
#             if seqid not in fasta_seqs:
#                 continue

#             full_seq = fasta_seqs[seqid]
#             region_seq = full_seq[max(0, start - flank): start + flank + 1].upper()
#             if strand == "-":
#                 region_seq = str(Seq(region_seq).reverse_complement())

#             sequences_all.append(region_seq)
#             matched_3  = pattern_3.match(region_seq)
#             matched_4a = pattern_4a.match(region_seq)
#             matched_4b = pattern_4b.match(region_seq)
#             matched_4c = pattern_4c.match(region_seq)
#             matched_5  = pattern_5.match(region_seq)
#             matched_3_all  = pattern_3_all.match(region_seq)
#             matched_4_all  = pattern_4_all.match(region_seq)
#             matched_5_all  = pattern_5_all.match(region_seq)


#             if matched_3_all:
#                 sequences_3_all.append(region_seq)
#             elif matched_4_all:
#                 sequences_4_all.append(region_seq)                
#             elif matched_5_all:
#                 sequences_5_all.append(region_seq)  
            
#             if matched_4a:
#                 sequences_4a.append(region_seq)
#             elif matched_4b:
#                 sequences_4b.append(region_seq)
#             elif matched_4c:
#                 sequences_4c.append(region_seq)
#             elif matched_3:
#                 sequences_3.append(region_seq)
#             elif matched_5:
#                 sequences_5.append(region_seq)   
#             else:
#                 sequences_else.append(region_seq)
#                 else_gff_lines.append(line)

#     print(f"Total sequences: {len(sequences_all)}")

#     print(f"3-gap all matches: {len(sequences_3_all)}")
#     print(f"3-gap matches: {len(sequences_3)}")

#     print(f"4-gap all matches: {len(sequences_4_all)}")
#     print(f"4-gap matches, just T: {len(sequences_4a)}")
#     print(f"4-gap matches, TG: {len(sequences_4b)}")
#     print(f"4-gap matches, TA: {len(sequences_4c)}")

#     print(f"5-gap all matches: {len(sequences_5_all)}")
#     print(f"5-gap matches: {len(sequences_5)}")

#     print(f"else (no match): {len(sequences_else)}")
#     return sequences_all, sequences_3, sequences_4a, sequences_4b, sequences_4c, sequences_5, sequences_else, else_gff_lines, sequences_3_all, sequences_4_all, sequences_5_all

# def save_logo_logomaker(sequences, output_file):
#     if not sequences:
#         print(f"No sequences found for {output_file}, skipping logo.")
#         return
#     # Convert the alignment to an information content matrix (bits)
#     info_df = logomaker.alignment_to_matrix(sequences, to_type='information')

#     if "4_all" in str(output_file):
#         print(sequences)


#     # Set up the logo figure
#     plt.figure(figsize=(len(info_df)/3, 3))
#     logo = logomaker.Logo(info_df, color_scheme='classic')
#     logo.style_spines(visible=False)
#     logo.style_spines(spines=['left', 'bottom'], visible=True)
#     logo.ax.set_ylabel('bits')
#     logo.ax.set_ylim([0, 2.0])  # Y-axis from 0 to 2 bits
#     plt.title(f"Sequence Logo: {output_file}")
#     plt.tight_layout()
#     # Save as PNG
#     plt.savefig(output_file + ".png", dpi=200)
#     plt.close()

# def main():
#     args = parse_args()
#     (seq_all, seq_3,seq_4a,
#      seq_4b, seq_4c, seq_5, seq_else, else_gff_lines, seq_3_all, seq_4_all, seq_5_all) = extract_polyA_sequences(
#         gff3_in=args.gff3_in,
#         fasta_file=args.fasta_file,
#         flank=24
#     )

#     # Save logos as PNG (will be named e.g., outAll.png, out3.png, ...)
#     save_logo_logomaker(seq_all, args.out_all)
#     save_logo_logomaker(seq_3, args.out_3)
#     save_logo_logomaker(seq_4a, args.out_4a)
#     save_logo_logomaker(seq_4b, args.out_4b)
#     save_logo_logomaker(seq_4c, args.out_4c)
#     save_logo_logomaker(seq_5, args.out_5)
#     save_logo_logomaker(seq_else, args.out_else)
#     save_logo_logomaker(seq_3_all, args.out_3_all)
#     save_logo_logomaker(seq_4_all, args.out_4_all)
#     save_logo_logomaker(seq_5_all, args.out_5_all)

#     # Write else GFF3
#     with open(args.gff3_else, "w") as f:
#         for gff_line in else_gff_lines:
#             f.write(gff_line + "\n")

# if __name__ == "__main__":
#     main()






#!/usr/bin/env python

###For has_stop category:
# python logo_generation_separated_jan9.py \
#     -g1 march4_has_stop_polyA.gff3 \
#     -g2 else_polyAs_has_stop_logo.gff3 \
#     -f ST2_sorted_masked.fasta \
#     --outAll has_stop_all \
#     --out3 has_stop_gap3 \
#     --out4a has_stop_gap4a \
#     --out4b has_stop_gap4b \
#     --out4c has_stop_gap4c \
#     --out5 has_stop_gap5 \
#     --outElse has_stop_else

###Results for has_stop category:
# Total sequences: 5163
# 3-gap matches: 52
# 4-gap matches, just T: 1828
# 4-gap matches, TG: 372
# 4-gap matches, TA: 215
# 5-gap matches: 110
# else (no match): 2586

###Results for non_stop category:
# Total sequences: 968
# 3-gap matches: 3
# 4-gap matches, just T: 673
# 4-gap matches, TG: 86
# 4-gap matches, TA: 88
# 5-gap matches: 5
# else (no match): 113

# import argparse
# from Bio import SeqIO
# from Bio.Seq import Seq
# import regex
# import pandas as pd
# import logomaker
# import matplotlib.pyplot as plt
# import io

# def parse_args():
#     parser = argparse.ArgumentParser(description="Generate sequence logos (PNG) and write 'else' GFF3 lines.")
#     parser.add_argument("-g1", "--gff3_in", dest="gff3_in", type=str, required=True,
#                         help="Input GFF3 file with polyA site info.")
#     parser.add_argument("-g2", "--gff3_else", dest="gff3_else", type=str, required=True,
#                         help="Output GFF3 for 'else' category.")
#     parser.add_argument("-f", "--fasta", dest="fasta_file", type=str, required=True,
#                         help="Input FASTA file of genomic sequences.")
#     parser.add_argument("--outAll",  dest="out_all",  required=True, help="Logo for all sequences.")
#     parser.add_argument("--out3",    dest="out_3",    required=True, help="Logo for 3-gap sequences.")
#     parser.add_argument("--out4a",   dest="out_4a",   required=True, help="Logo for 4-gap (just T).")
#     parser.add_argument("--out4b",   dest="out_4b",   required=True, help="Logo for 4-gap (TG).")
#     parser.add_argument("--out4c",   dest="out_4c",   required=True, help="Logo for 4-gap (TA).")
#     parser.add_argument("--out5",    dest="out_5",    required=True, help="Logo for 5-gap sequences.")
#     parser.add_argument("--outElse", dest="out_else", required=True, help="Logo for 'else' sequences.")
#     return parser.parse_args()

# def extract_polyA_sequences(gff3_in, fasta_file, flank=24):
#     sequences_all  = []
#     sequences_3    = []
#     sequences_4a   = []
#     sequences_4b   = []
#     sequences_4c   = []
#     sequences_5    = []
#     sequences_else = []
#     else_gff_lines = []

#     pattern_3  = regex.compile(r"^(.{24})(T)([ACGT]{3})(TGTTTGTT){s<=2}")
#     pattern_4a = regex.compile(r"^(.{24})(T)([ACGT]{4})(TGTTTGTT){s<=2}")
#     pattern_4b = regex.compile(r"^(.{23})(TG)([ACGT]{4})(TGTTTGTT){s<=2}")
#     pattern_4c = regex.compile(r"^(.{23})(TA)([ACGT]{4})(TGTTTGTT){s<=2}")
#     pattern_5  = regex.compile(r"^(.{24})(T)([CT]{1}[ACGT]{4})(TGTTTGTT){s<=2}")

#     fasta_seqs = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}

#     with open(gff3_in, "r") as gf:
#         for line in gf:
#             line = line.strip()
#             if line.startswith("#"):
#                 continue
#             parts = line.split("\t")
#             if len(parts) < 9 or parts[2] != "polyA":
#                 continue

#             seqid = parts[0]
#             start = int(parts[3]) - 1
#             strand = parts[6]
#             if seqid not in fasta_seqs:
#                 continue

#             full_seq = fasta_seqs[seqid]
#             region_seq = full_seq[max(0, start - flank): start + flank + 1].upper()
#             if strand == "-":
#                 region_seq = str(Seq(region_seq).reverse_complement())

#             sequences_all.append(region_seq)
#             matched_3  = pattern_3.match(region_seq)
#             matched_4a = pattern_4a.match(region_seq)
#             matched_4b = pattern_4b.match(region_seq)
#             matched_4c = pattern_4c.match(region_seq)
#             matched_5  = pattern_5.match(region_seq)

#             if matched_4a:
#                 sequences_4a.append(region_seq)
#             elif matched_4b:
#                 sequences_4b.append(region_seq)
#             elif matched_4c:
#                 sequences_4c.append(region_seq)
#             elif matched_3:
#                 sequences_3.append(region_seq)
#             elif matched_5:
#                 sequences_5.append(region_seq)
#             else:
#                 sequences_else.append(region_seq)
#                 else_gff_lines.append(line)

#     print(f"Total sequences: {len(sequences_all)}")
#     print(f"3-gap matches: {len(sequences_3)}")
#     print(f"4-gap matches, just T: {len(sequences_4a)}")
#     print(f"4-gap matches, TG: {len(sequences_4b)}")
#     print(f"4-gap matches, TA: {len(sequences_4c)}")
#     print(f"5-gap matches: {len(sequences_5)}")
#     print(f"else (no match): {len(sequences_else)}")
#     return sequences_all, sequences_3, sequences_4a, sequences_4b, sequences_4c, sequences_5, sequences_else, else_gff_lines

# def save_logo_logomaker(sequences, output_file):
#     if not sequences:
#         print(f"No sequences found for {output_file}, skipping logo.")
#         return
#     # Convert the alignment to an information content matrix (bits)
#     info_df = logomaker.alignment_to_matrix(sequences, to_type='information')

#     # Set up the logo figure
#     plt.figure(figsize=(len(info_df)/3, 3))
#     logo = logomaker.Logo(info_df, color_scheme='classic')
#     logo.style_spines(visible=False)
#     logo.style_spines(spines=['left', 'bottom'], visible=True)
#     logo.ax.set_ylabel('bits')
#     logo.ax.set_ylim([0, 2.0])  # Y-axis from 0 to 2 bits
#     plt.title(f"Sequence Logo: {output_file}")
#     plt.tight_layout()
#     # Save as PNG
#     plt.savefig(output_file + ".png", dpi=200)
#     plt.close()

# def main():
#     args = parse_args()
#     (seq_all, seq_3,seq_4a,
#      seq_4b, seq_4c, seq_5, seq_else, else_gff_lines) = extract_polyA_sequences(
#         gff3_in=args.gff3_in,
#         fasta_file=args.fasta_file,
#         flank=24
#     )

#     # Save logos as PNG (will be named e.g., outAll.png, out3.png, ...)
#     save_logo_logomaker(seq_all, args.out_all)
#     save_logo_logomaker(seq_3, args.out_3)
#     save_logo_logomaker(seq_4a, args.out_4a)
#     save_logo_logomaker(seq_4b, args.out_4b)
#     save_logo_logomaker(seq_4c, args.out_4c)
#     save_logo_logomaker(seq_5, args.out_5)
#     save_logo_logomaker(seq_else, args.out_else)

#     # Write else GFF3
#     with open(args.gff3_else, "w") as f:
#         for gff_line in else_gff_lines:
#             f.write(gff_line + "\n")

# if __name__ == "__main__":
#     main()




#### This script worked on the lab's computer but not on my laptop, but the new version works on my laptop:
# # $ python logo_generation_separated_jan9.py -g1 march4_non_stop_polyA.gff3 -f ST2_sorted_masked.fasta --outAll non_stop_all --out3 non_stop_gap3 --out4a non_stop_gap4a --out4b non_stop_gap4b --out4c non_stop_gap4c --out5 non_stop_gap5 --outElse non_stop_else -g2 else_polyAs_logo.gff3
# # Total sequences: 968
# # 3-gap matches: 3
# # 4-gap matches, just T: 673
# # 4-gap matches, TG: 86
# # 4-gap matches, TA: 88
# # 5-gap matches: 10
# # else (no match): 108

# #!/usr/bin/env python
# import argparse
# from Bio import SeqIO
# from Bio.Seq import Seq
# import weblogo as wl
# import regex
# import io

# def parse_args():
#     parser = argparse.ArgumentParser(
#         description="Generate separate sequence logos and save 'else' category lines to a separate GFF3."
#     )
#     # GFF3 input (main polyA sites)
#     parser.add_argument(
#         "-g1", "--gff3_in",
#         dest="gff3_in",
#         type=str,
#         required=True,
#         help="Input GFF3 file with polyA site information."
#     )
#     # GFF3 output (for 'else' category lines)
#     parser.add_argument(
#         "-g2", "--gff3_else",
#         dest="gff3_else",
#         type=str,
#         required=True,
#         help="Output GFF3 file to store lines that don't match any pattern."
#     )
#     # FASTA input
#     parser.add_argument(
#         "-f", "--fasta",
#         dest="fasta_file",
#         type=str,
#         required=True,
#         help="Input FASTA file containing genomic sequences."
#     )
#     # Output logos
#     parser.add_argument("--outAll",  dest="out_all",  required=True, help="Logo for all sequences.")
#     parser.add_argument("--out3",    dest="out_3",    required=True, help="Logo for 3-gap sequences.")
#     parser.add_argument("--out4a",   dest="out_4a",   required=True, help="Logo for 4-gap (just T).")
#     parser.add_argument("--out4b",   dest="out_4b",   required=True, help="Logo for 4-gap (TG).")
#     parser.add_argument("--out4c",   dest="out_4c",   required=True, help="Logo for 4-gap (TA).")
#     parser.add_argument("--out5",    dest="out_5",    required=True, help="Logo for 5-gap sequences.")
#     parser.add_argument("--outElse", dest="out_else", required=True, help="Logo for 'else' category sequences.")

#     return parser.parse_args()

# def extract_polyA_sequences(gff3_in, fasta_file, flank=24):
#     """
#     Extract sequences around polyA sites, classify them by various fuzzy patterns,
#     and store lines that match "else" category in a separate list for writing later.
#     """
#     sequences_all  = []
#     sequences_3    = []
#     sequences_4a   = []
#     sequences_4b   = []
#     sequences_4c   = []
#     sequences_5    = []
#     sequences_else = []
#     else_gff_lines = []

#     # Patterns for 3-gap, 4-gap, 5-gap
#     pattern_3  = regex.compile(r"^(.{24})(T)([ACGT]{3})(TGTTTGTT){s<=2}")
#     pattern_4a = regex.compile(r"^(.{24})(T)([ACGT]{4})(TGTTTGTT){s<=2}")
#     pattern_4b = regex.compile(r"^(.{23})(TG)([ACGT]{4})(TGTTTGTT){s<=2}")
#     pattern_4c = regex.compile(r"^(.{23})(TA)([ACGT]{4})(TGTTTGTT){s<=2}")
#     pattern_5  = regex.compile(r"^(.{24})(T)([ACGT]{5})(TGTTTGTT){s<=2}")

#     # Load FASTA
#     fasta_seqs = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}

#     with open(gff3_in, "r") as gf:
#         for line in gf:
#             line = line.strip()
#             if line.startswith("#"):
#                 continue
#             parts = line.split("\t")
#             if len(parts) < 9 or parts[2] != "polyA":
#                 continue

#             seqid = parts[0]
#             start = int(parts[3]) - 1
#             strand = parts[6]

#             if seqid not in fasta_seqs:
#                 continue

#             full_seq = fasta_seqs[seqid]
#             region_seq = full_seq[max(0, start - flank) : start + flank + 1].upper()

#             # Reverse complement if negative
#             if strand == "-":
#                 region_seq = str(Seq(region_seq).reverse_complement())

#             sequences_all.append(region_seq)

#             # Check patterns
#             matched_3  = pattern_3.match(region_seq)
#             matched_4a = pattern_4a.match(region_seq)
#             matched_4b = pattern_4b.match(region_seq)
#             matched_4c = pattern_4c.match(region_seq)
#             matched_5  = pattern_5.match(region_seq)

#             if matched_4a:
#                 sequences_4a.append(region_seq)
#             elif matched_4b:
#                 sequences_4b.append(region_seq)
#             elif matched_4c:
#                 sequences_4c.append(region_seq)
#             elif matched_3:
#                 sequences_3.append(region_seq)
#             elif matched_5:
#                 sequences_5.append(region_seq)
#             else:
#                 sequences_else.append(region_seq)
#                 else_gff_lines.append(line)

#     # Print summary
#     print(f"Total sequences: {len(sequences_all)}")
#     print(f"3-gap matches: {len(sequences_3)}")
#     print(f"4-gap matches, just T: {len(sequences_4a)}")
#     print(f"4-gap matches, TG: {len(sequences_4b)}")
#     print(f"4-gap matches, TA: {len(sequences_4c)}")
#     print(f"5-gap matches: {len(sequences_5)}")
#     print(f"else (no match): {len(sequences_else)}")

#     return sequences_all, sequences_3, sequences_4a, sequences_4b, sequences_4c, sequences_5, sequences_else, else_gff_lines

# def generate_logo(sequences, output_file):
#     """
#     Generate a sequence logo from a list of sequences and save it as an image.
#     """
#     if not sequences:
#         print(f"No sequences found for {output_file}, skipping logo.")
#         return

#     seq_data = "\n".join(sequences)
#     seqs = wl.read_seq_data(io.StringIO(seq_data))

#     logo_data = wl.LogoData.from_seqs(seqs)
#     logo_options = wl.LogoOptions()
#     logo_options.title = f"Sequence Logo: {output_file}"
#     logo_options.color_scheme = wl.std_color_schemes['classic']

#     logo_format = wl.LogoFormat(logo_data, logo_options)
#     png = wl.png_print_formatter(logo_data, logo_format)

#     with open(output_file, "wb") as out:
#         out.write(png)

# def main():
#     args = parse_args()

#     # 1. Extract & classify sequences
#     (seq_all, seq_3,seq_4a,
#      seq_4b, seq_4c, seq_5, seq_else, else_gff_lines) = extract_polyA_sequences(
#         gff3_in=args.gff3_in,
#         fasta_file=args.fasta_file,
#         flank=24
#     )

#     # 2. Generate logos for each group
#     generate_logo(seq_all, f"{args.out_all}.png")
#     generate_logo(seq_3, f"{args.out_3}.png")
#     generate_logo(seq_4a, f"{args.out_4a}.png")
#     generate_logo(seq_4b, f"{args.out_4b}.png")
#     generate_logo(seq_4c, f"{args.out_4c}.png")
#     generate_logo(seq_5, f"{args.out_5}.png")
#     generate_logo(seq_else, f"{args.out_else}.png")

#     # 3. Write "else" category lines to an output GFF3
#     with open(args.gff3_else, "w") as f:
#         for gff_line in else_gff_lines:
#             f.write(gff_line + "\n")

# if __name__ == "__main__":
#     main()
