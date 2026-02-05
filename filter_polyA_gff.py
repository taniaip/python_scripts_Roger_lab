#!/usr/bin/env python3

"""
Collapse nearby poly(A) sites by retaining only the feature with the highest
read count within each local cluster.

This script reads a sorted polyA GFF3 file and compares adjacent polyA features
on the same contig and strand. If two neighboring sites are within 4 nucleotides
of each other, they are treated as belonging to the same cluster (“bucket”), and
only the site with the higher read count is kept (ties resolved by later position).

Read counts are extracted from the polyA feature ID, assumed to follow the format:
  polyA_site_<index>_<strand>_readcount_<N>

After collapsing nearby sites, duplicate entries are removed and the filtered
polyA features are written to a new GFF3 file.

Inputs:
  --polyA           Sorted GFF3 of polyA sites
  --filtered_polyA Output GFF3 path

Output:
  GFF3 containing de-duplicated polyA sites with local clusters collapsed.

Dependencies:
  gffutils

Author: Tania Iranpour
"""

import argparse
import gffutils

# extract readcount from ID like: polyA_site_00421_-_readcount_55  -> 55
def find_readcount(id_str: str) -> int:
    parts = id_str.split('_')
    return int(parts[5])


def same_bucket(a, b) -> bool:
    return (a.seqid == b.seqid) and (a.strand == b.strand)

def main():
    ap = argparse.ArgumentParser(description="Replace nearby polyA sites with the higher-readcount line; dedupe at end.")
    ap.add_argument("-g", "--polyA", required=True, help="Input GFF3 of polyA sites (already sorted).")
    ap.add_argument("--filtered_polyA", required=True, help="Output GFF3 path.")
    args = ap.parse_args()

    db = gffutils.create_db(args.polyA, dbfn=":memory:", force=True,
                            keep_order=True, merge_strategy="create_unique")

    # put features into a list that we can index
    feats = list(db.features_of_type("polyA"))

    # compare neighbors, if they're in the same bucket, replace both with the "winner"(has the higher readcount)
    length_db = len(feats)
    for i in range(0, length_db - 1):
        curr = feats[i]
        next = feats[i + 1]

        if same_bucket(curr, next) and (next.start - curr.start) <= 4:
            rc_curr = find_readcount(curr.id)
            rc_next = find_readcount(next.id)

            # kept_feat: higher readcount; tie -> smaller start
            if (rc_next > rc_curr) or (rc_next == rc_curr):
                kept_feat = next
            else:
                kept_feat = curr

            # replace both lines with the chosen winner
            feats[i] = kept_feat
            feats[i + 1] = kept_feat

    # To keep one copy of the duplicated lines
    seen = set()
    kept_lines = []
    for f in feats:
        line = str(f)
        if line not in seen:
            seen.add(line)
            kept_lines.append(line)

    # Write output
    with open(args.filtered_polyA, "w") as out:
        for line in kept_lines:
            out.write(line + "\n")

if __name__ == "__main__":
    main()

