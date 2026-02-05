#!/usr/bin/env python

# python categorize_polyA.py -g1 ST2_orfs.gff3 -g2 polyAsites_updated.gff3 -f ST2_sorted_masked.fasta --output_has_stop march4_has_stop_polyA.gff3
#  --output_non_stop march4_non_stop_polyA.gff3 --output_within_gene march4_within_genes_polyA.gff3 --output_not_matched march4_not_matched_polyA.gff3

"""
Categorize poly(A) sites relative to annotated genes and write separate GFF3 outputs.

This script builds in-memory gffutils databases from:
  1) a gene-annotation GFF3, and
  2) a polyA-site GFF3,
then uses the reference genome FASTA to classify each poly(A) feature into one of
four categories:

1) has_stop:
   Poly(A) site is assigned to the nearest non-overlapped (“big”) gene on the same
   strand, and the gene’s terminal codon is a canonical stop (TAA/TGA/TAG on + strand;
   reverse-complement equivalents on - strand).

2) non_stop:
   Poly(A) site is assigned to the nearest non-overlapped (“big”) gene on the same
   strand, but the gene does not end with a canonical stop codon at the expected end.

3) within_gene:
   Poly(A) site is not assignable to a flanking gene interval and remains within an
   annotated gene span.

4) not_matched:
   Poly(A) site falls in an intergenic interval but does not match the expected strand
   of the neighboring gene(s), so it cannot be assigned confidently.

Key steps:
- Genes are grouped per contig, sorted by start (ties broken by longer end), and
  overlapping smaller genes are skipped to define clean intergenic intervals.
- Each poly(A) site is tested against these intervals with a small ±3 nt tolerance.
- For assigned sites, the output GFF3 attributes are augmented with the associated
  gene ID, gene coordinates/strand, and the assigned category.

Inputs:
  --genes            Gene annotation GFF3
  --polyA            PolyA site GFF3
  --fasta            Reference genome FASTA
  --output_has_stop  GFF3 of assigned sites whose genes end in a stop codon
  --output_non_stop  GFF3 of assigned sites whose genes lack a terminal stop codon
  --output_within_gene   GFF3 of polyA sites within genes
  --output_not_matched   GFF3 of unassigned intergenic polyA sites

Dependencies:
  gffutils, Biopython

Author: Tania Iranpour
"""

import sys
import argparse

import gffutils
from Bio import SeqIO

from typing import List, Dict, Tuple, Optional
from gffutils.feature import Feature
from gffutils.interface import FeatureDB
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict #for adding a column to the output gff3 file for the gene ID of the associated gene.

# Argument Parser
parser = argparse.ArgumentParser(description="Categorize polyA sites into genomic contexts.")
parser.add_argument("-g1", "--genes", required=True, help="Input GFF3 file with gene information.")
parser.add_argument("-g2", "--polyA", required=True, help="Input GFF3 file with polyA site information.")
parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file containing genomic sequences.")
parser.add_argument("--output_has_stop", required=True, help="Output file for polyA sites associated with genes having stop codons.")
parser.add_argument("--output_non_stop", required=True, help="Output file for polyA sites associated with genes not having stop codons.")
parser.add_argument("--output_within_gene", required=True, help="Output file for polyA sites within genes.")
parser.add_argument("--output_not_matched", required=True, help="Output file for polyA sites not associated with any genes.")
args = parser.parse_args()


class Interval:
    """Class representing an interval between two genes."""
    def __init__(self, contig: str, gene_left: Optional[Feature], gene_right: Optional[Feature]):
        self.contig = contig
        self.start = gene_left.end if gene_left else 0
        self.end = gene_right.start if gene_right else sys.maxsize
        self.gene_left = gene_left
        self.gene_right = gene_right
        self.strand_left = gene_left.strand if gene_left else None
        self.strand_right = gene_right.strand if gene_right else None


def main() -> None:
    db_genes: FeatureDB = gffutils.create_db(args.genes, dbfn=":memory:", merge_strategy="create_unique")
    db_polyA: FeatureDB = gffutils.create_db(args.polyA, dbfn=":memory:", merge_strategy="create_unique")
    sequences = SeqIO.index(args.fasta, 'fasta')

    # Generate intervals
    gene_intervals = process_gene_intervals(db_genes)

    # Assign polyA to categories and write directly to files
    assign_polyA_to_genes(
        gene_intervals,
        db_polyA,
        sequences,
        output_has_stop=args.output_has_stop,
        output_non_stop=args.output_non_stop,
        output_within_gene=args.output_within_gene,
        output_not_matched=args.output_not_matched
    )

def process_gene_intervals(db: FeatureDB) -> List[Interval]:
    """
    Generate intervals between 'big' genes in the GFF3 database,
    skipping smaller genes overlapped by a bigger gene.
    """
    intervals: List[Interval] = []
    genes_by_contig: Dict[str, List[Feature]] = {}

    # Collect genes by contig
    for gene in db.features_of_type("gene"):
        if gene.seqid not in genes_by_contig:
            genes_by_contig[gene.seqid] = []
        genes_by_contig[gene.seqid].append(gene)

    # Process each contig individually
    for contig, genes in genes_by_contig.items():
        # Sort genes by their start positions and in case of starting at the same point prioritise the ones with further ends
        genes.sort(key=lambda g: (g.start, -g.end))

        # Create interval from "start of genome" to the first gene
        first_gene = genes[0]
        intervals.append(Interval(contig=contig, gene_left=None, gene_right=first_gene))
        correct_start = first_gene

        # Iterate over the remaining genes
        for i in range(1, len(genes)):
            # If next gene is overshadowed by correct_start (start < correct_start.end), skip it
            if genes[i].start < correct_start.end:
                continue

            # Otherwise, create an interval from the last big gene to this new gene
            intervals.append(Interval(contig=contig, gene_left=correct_start, gene_right=genes[i]))

            correct_start = genes[i]

        # Create interval from the last big gene to "end of genome"
        intervals.append(Interval(contig=contig, gene_left=correct_start, gene_right=None))

    return intervals


def has_stop_codon(gene: Feature, sequences: Dict[str, SeqRecord]) -> bool:
    """Check if a gene has a stop codon based on its strand."""
    sequence = str(sequences[gene.seqid].seq[gene.start-1 : gene.end]).upper()
    if gene.strand == "+":
        return sequence[-3:] in ["TAA", "TGA", "TAG"]
    if gene.strand == "-":
        return sequence[:3] in ["TTA", "TCA", "CTA"]


# To add a column to the output gff3 file for the gene ID that each polyA in has_stop and non_stop categories are associated with.
#written by ChatGPT
def _format_attributes(attr_dict: Dict[str, List[str]]) -> str:
    """Turn a gffutils attributes dict into a GFF3 'key=val;key2=val2' string."""
    parts = []
    for k, vals in attr_dict.items():
        for v in vals:
            parts.append(f"{k}={v}")
    return ";".join(parts)

def _write_polyA_with_attrs(polyA: Feature, out_fh, extra_attrs: Dict[str, str]) -> None:
    """Write a polyA feature line, adding/updating attributes."""
    # Start from existing attributes, append ours
    merged = OrderedDict((k, list(v)) for k, v in polyA.attributes.items())
    for k, v in extra_attrs.items():
        merged[k] = [str(v)]
    attr_str = _format_attributes(merged)

    # Keep columns 1–8 intact (score/phase stay spec-compliant)
    line = "\t".join([
        polyA.seqid,
        polyA.source,
        polyA.featuretype,
        str(polyA.start),
        str(polyA.end),
        polyA.score if polyA.score is not None else ".",
        polyA.strand if polyA.strand else ".",
        polyA.frame if polyA.frame not in (None, "") else ".",
        attr_str if attr_str else ".",
    ]) + "\n"
    out_fh.write(line)
###end


def assign_polyA_to_genes(
    gene_intervals: List[Interval],
    db_polyA: FeatureDB,
    sequences: Dict[str, SeqRecord],
    output_has_stop: str,
    output_non_stop: str,
    output_within_gene: str,
    output_not_matched: str,
) -> None:
    """Assign polyA sites to genes, categorize them, and write directly to output files."""
    with open(output_has_stop, "w") as has_stop_file, \
         open(output_non_stop, "w") as non_stop_file, \
         open(output_within_gene, "w") as within_file, \
         open(output_not_matched, "w") as unmatched_file:

        for polyA in db_polyA.features_of_type('polyA'):
            matched = False
            within_gene = True

            for interval in gene_intervals:
                if polyA.seqid != interval.contig:
                    continue

                # Check if polyA falls within an interval and the strands match
                if interval.start <= polyA.start <= interval.end + 3:
                    within_gene = False
                    if polyA.strand == "-" and interval.strand_right == "-":
                        matched = True
                        gene = interval.gene_right
                        if has_stop_codon(gene, sequences):
                            _write_polyA_with_attrs(polyA, has_stop_file, {
                                "associated_gene": gene.id,
                                "gene_start": gene.start,
                                "gene_end": gene.end,
                                "gene_strand": gene.strand,
                                "category": "has_stop"
                            })
                        else:
                            _write_polyA_with_attrs(polyA, non_stop_file, {
                                "associated_gene": gene.id,
                                "gene_start": gene.start,
                                "gene_end": gene.end,
                                "gene_strand": gene.strand,
                                "category": "non_stop"
                            })

                if interval.start - 3 <= polyA.start <= interval.end:
                    within_gene = False
                    if polyA.strand == "+" and interval.strand_left == "+":
                        matched = True
                        gene = interval.gene_left
                        if has_stop_codon(gene, sequences):
                            _write_polyA_with_attrs(polyA, has_stop_file, {
                                "associated_gene": gene.id,
                                "gene_start": gene.start,
                                "gene_end": gene.end,
                                "gene_strand": gene.strand,
                                "category": "has_stop"
                            })
                        else:
                            _write_polyA_with_attrs(polyA, non_stop_file, {
                                "associated_gene": gene.id,
                                "gene_start": gene.start,
                                "gene_end": gene.end,
                                "gene_strand": gene.strand,
                                "category": "non_stop"
                            })    

            if not matched and within_gene:
                within_file.write(str(polyA) + "\n")
            elif not matched:
                unmatched_file.write(str(polyA) + "\n")


if __name__ == "__main__":
    main()




# #!/usr/bin/env python
# import sys
# import argparse

# import gffutils
# from Bio import SeqIO

# from typing import List, Dict, Tuple, Optional
# from gffutils.feature import Feature
# from gffutils.interface import FeatureDB
# from Bio.SeqRecord import SeqRecord

# # Argument Parser
# parser = argparse.ArgumentParser(description="Categorize polyA sites into genomic contexts.")
# parser.add_argument("-g1", "--genes", required=True, help="Input GFF3 file with gene information.")
# parser.add_argument("-g2", "--polyA", required=True, help="Input GFF3 file with polyA site information.")
# parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file containing genomic sequences.")
# parser.add_argument("--output_has_stop", required=True, help="Output file for polyA sites associated with genes having stop codons.")
# parser.add_argument("--output_non_stop", required=True, help="Output file for polyA sites associated with genes not having stop codons.")
# parser.add_argument("--output_within_gene", required=True, help="Output file for polyA sites within genes.")
# parser.add_argument("--output_not_matched", required=True, help="Output file for polyA sites not associated with any genes.")
# args = parser.parse_args()


# class Interval:
#     """Class representing an interval between two genes."""
#     def __init__(self, contig: str, gene_left: Optional[Feature], gene_right: Optional[Feature]):
#         self.contig = contig
#         self.start = gene_left.end if gene_left else 0
#         self.end = gene_right.start if gene_right else sys.maxsize
#         self.gene_left = gene_left
#         self.gene_right = gene_right
#         self.strand_left = gene_left.strand if gene_left else None
#         self.strand_right = gene_right.strand if gene_right else None


# def main() -> None:
#     db_genes: FeatureDB = gffutils.create_db(args.genes, dbfn=":memory:", merge_strategy="create_unique")
#     db_polyA: FeatureDB = gffutils.create_db(args.polyA, dbfn=":memory:", merge_strategy="create_unique")
#     sequences = SeqIO.index(args.fasta, 'fasta')

#     # Generate intervals
#     gene_intervals = process_gene_intervals(db_genes)

#     # Assign polyA to categories and write directly to files
#     assign_polyA_to_genes(
#         gene_intervals,
#         db_polyA,
#         sequences,
#         output_has_stop=args.output_has_stop,
#         output_non_stop=args.output_non_stop,
#         output_within_gene=args.output_within_gene,
#         output_not_matched=args.output_not_matched
#     )


# def process_gene_intervals(db: FeatureDB) -> List[Interval]:
#     """Generate intervals between genes in the GFF3 database."""
#     intervals = []
#     genes_by_contig: Dict[str, List[Feature]] = {}

#     for gene in db.features_of_type("gene"):
#         if gene.seqid not in genes_by_contig:
#             genes_by_contig[gene.seqid] = []
#         genes_by_contig[gene.seqid].append(gene)

#     for contig, genes in genes_by_contig.items():
#         genes.sort(key=lambda g: g.start)

#         # First gene interval
#         first_gene = genes[0]
#         intervals.append(Interval(contig = contig, gene_left = None, gene_right = first_gene))

#         # Gene-to-gene intervals
#         for i in range(len(genes) - 1):
#             intervals.append(Interval(contig = contig, gene_left = genes[i], gene_right = genes[i + 1]))

#         # Last gene interval
#         last_gene = genes[-1]
#         intervals.append(Interval(contig = contig, gene_left = last_gene, gene_right = None))

#     return intervals


# def has_stop_codon(gene: Feature, sequences: Dict[str, SeqRecord]) -> bool:
#     """Check if a gene has a stop codon based on its strand."""
#     sequence = str(sequences[gene.seqid].seq[gene.start-1 : gene.end]).upper()
#     if gene.strand == "+":
#         return sequence[-3:] in ["TAA", "TGA", "TAG"]
#     if gene.strand == "-":
#         return sequence[:3] in ["TTA", "TCA", "CTA"]


# def assign_polyA_to_genes(
#     gene_intervals: List[Interval],
#     db_polyA: FeatureDB,
#     sequences: Dict[str, SeqRecord],
#     output_has_stop: str,
#     output_non_stop: str,
#     output_within_gene: str,
#     output_not_matched: str,
# ) -> None:
#     """Assign polyA sites to genes, categorize them, and write directly to output files."""
#     with open(output_has_stop, "w") as has_stop_file, \
#          open(output_non_stop, "w") as non_stop_file, \
#          open(output_within_gene, "w") as within_file, \
#          open(output_not_matched, "w") as unmatched_file:

#         for polyA in db_polyA.features_of_type('polyA'):
#             matched = False
#             within_gene = True

#             for interval in gene_intervals:
#                 if polyA.seqid != interval.contig:
#                     continue

#                 # Check if polyA falls within an interval and the strands match
#                 if interval.start <= polyA.start <= interval.end + 3:
#                     within_gene = False
#                     if polyA.strand == "-" and interval.strand_right == "-":
#                         matched = True
#                         gene = interval.gene_right
#                         if has_stop_codon(gene, sequences):
#                             has_stop_file.write(str(polyA) + "\n")
#                         else:
#                             non_stop_file.write(str(polyA) + "\n")

#                 if interval.start - 3 <= polyA.start <= interval.end:
#                     within_gene = False
#                     if polyA.strand == "+" and interval.strand_left == "+":
#                         matched = True
#                         gene = interval.gene_left
#                         if has_stop_codon(gene, sequences):
#                             has_stop_file.write(str(polyA) + "\n")
#                         else:
#                             non_stop_file.write(str(polyA) + "\n")

#             if not matched and within_gene:
#                 within_file.write(str(polyA) + "\n")
#             elif not matched:
#                 unmatched_file.write(str(polyA) + "\n")


# if __name__ == "__main__":
#     main()








# ### Everything is perefctly edited based on the 4th feedback on Jan 30;
# #!/usr/bin/env python
# import sys
# import argparse

# import gffutils
# from Bio import SeqIO

# from typing import List, Dict, Tuple, Optional
# from gffutils.feature import Feature
# from gffutils.interface import FeatureDB
# from Bio.SeqRecord import SeqRecord

# # Argument Parser
# parser = argparse.ArgumentParser(description="Categorize polyA sites into genomic contexts.")
# parser.add_argument("-g1", "--genes", required=True, help="Input GFF3 file with gene information.")
# parser.add_argument("-g2", "--polyA", required=True, help="Input GFF3 file with polyA site information.")
# parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file containing genomic sequences.")
# parser.add_argument("--output_true", required=True, help="Output file for polyA sites associated with genes having stop codons.")
# parser.add_argument("--output_false", required=True, help="Output file for polyA sites associated with genes not having stop codons.")
# parser.add_argument("--output_within_gene", required=True, help="Output file for polyA sites within genes.")
# parser.add_argument("--output_not_matched", required=True, help="Output file for polyA sites not associated with any genes.")
# args = parser.parse_args()


# class Interval:
#     """Class representing an interval between two genes."""
#     def __init__(self, contig: str, gene_left: Optional[Feature], gene_right: Optional[Feature]):
#         self.contig = contig
#         self.start = gene_left.end if gene_left else 0
#         self.end = gene_right.start if gene_right else sys.maxsize
#         self.gene_left = gene_left
#         self.gene_right = gene_right
#         self.strand_left = gene_left.strand if gene_left else None
#         self.strand_right = gene_right.strand if gene_right else None


# def main() -> None:
#     db_genes: FeatureDB = gffutils.create_db(args.genes, dbfn=":memory:", merge_strategy="create_unique")
#     db_polyA: FeatureDB = gffutils.create_db(args.polyA, dbfn=":memory:", merge_strategy="create_unique")
#     sequences = SeqIO.index(args.fasta, 'fasta')

#     # Generate intervals and extract polyA features
#     gene_intervals = process_gene_intervals(db_genes)

#     # Assign polyA to categories
#     has_stop_matched_ids, non_stop_matched_ids, within_gene_polyA_ids, unmatched_polyA_ids = assign_polyA_to_genes(
#         gene_intervals, db_polyA, sequences
#     )

#     # Write output files
#     with open(args.output_true, "w") as true_file, open(args.output_false, "w") as false_file, open(args.output_within_gene, "w") as within_file, open(args.output_not_matched, "w") as unmatched_file:
#         for feature in db_polyA.features_of_type("polyA"):
#             if feature.id in has_stop_matched_ids:
#                 true_file.write(str(feature) + "\n")
#             elif feature.id in non_stop_matched_ids:
#                 false_file.write(str(feature) + "\n")
#             elif feature.id in within_gene_polyA_ids:
#                 within_file.write(str(feature) + "\n")
#             elif feature.id in unmatched_polyA_ids:
#                 unmatched_file.write(str(feature) + "\n")


# def process_gene_intervals(db: FeatureDB) -> List[Interval]:
#     """Generate intervals between genes in the GFF3 database."""
#     intervals = []
#     genes_by_contig: Dict[str, List[Feature]] = {}

#     for gene in db.features_of_type("gene"):
#         if gene.seqid not in genes_by_contig:
#             genes_by_contig[gene.seqid] = []
#         genes_by_contig[gene.seqid].append(gene)

#     for contig, genes in genes_by_contig.items():
#         genes.sort(key=lambda g: g.start)

#         # First gene interval
#         first_gene = genes[0]
#         intervals.append(Interval(contig = contig, gene_left = None, gene_right = first_gene))

#         # Gene-to-gene intervals
#         for i in range(len(genes) - 1):
#             intervals.append(Interval(contig = contig, gene_left = genes[i], gene_right = genes[i + 1]))

#         # Last gene interval
#         last_gene = genes[-1]
#         intervals.append(Interval(contig = contig, gene_left = last_gene, gene_right = None))

#     return intervals


# def has_stop_codon(gene: Feature, sequences: Dict[str, SeqRecord]) -> bool:
#     """Check if a gene has a stop codon based on its strand."""
#     sequence = str(sequences[gene.seqid].seq[gene.start-1 : gene.end]).upper()
#     if gene.strand == "+":
#         return sequence[-3:] in ["TAA", "TGA", "TAG"]
#     if gene.strand == "-":
#         return sequence[:3] in ["TTA", "TCA", "CTA"]

# def assign_polyA_to_genes(
#     gene_intervals: List[Interval],
#     db_polyA: FeatureDB,
#     sequences: Dict[str, SeqRecord],
# ) -> Tuple[List[str], List[str], List[str], List[str]]:
#     """Assign polyA sites to genes and categorize them."""
#     has_stop_matched = []
#     non_stop_matched = []
#     within_gene_polyA = []
#     unmatched_polyA = []

#     for polyA in db_polyA.features_of_type('polyA'):
#         matched = False
#         within_gene = True

#         for interval in gene_intervals:
#             if polyA.seqid != interval.contig:
#                 continue
                  
#             # Check if polyA falls within an interval and the strands match 
#             if interval.start <= polyA.start <= interval.end+3:
#                 within_gene = False
#                 if polyA.strand == "-" and interval.strand_right == "-":
#                     matched = True
#                     gene = interval.gene_right
#                     if has_stop_codon(gene, sequences):
#                         has_stop_matched.append(polyA.id)
#                     else:
#                         non_stop_matched.append(polyA.id)

#             if interval.start-3 <= polyA.start <= interval.end:
#                 within_gene = False
#                 if polyA.strand == "+" and interval.strand_left == "+":
#                     matched = True
#                     gene = interval.gene_left
#                     if has_stop_codon(gene, sequences):
#                         has_stop_matched.append(polyA.id)
#                     else:
#                         non_stop_matched.append(polyA.id)

#         if not matched and within_gene:
#             within_gene_polyA.append(polyA.id)
#         else:
#             unmatched_polyA.append(polyA.id)

#     return has_stop_matched, non_stop_matched, within_gene_polyA, unmatched_polyA


# if __name__ == "__main__":
#     main()



# #The following script works perfectly fine. The new script is written for the 3rd feedback from Joran. 
# 
# !/usr/bin/env python

# from typing import List, Dict, Union, Tuple
# import argparse
# import gffutils
# from gffutils.feature import Feature
# from gffutils.interface import FeatureDB
# from Bio import SeqIO

# parser = argparse.ArgumentParser(description="Categorize polyA sites into 'has stop codon', 'not having stop codon', and 'not associated with any genes'.")
# parser.add_argument("-g1", "--genes", required=True, help="Input GFF3 file with gene information.")
# parser.add_argument("-g2", "--polyA", required=True, help="Input GFF3 file with polyA site information.")
# parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file containing genomic sequences.")
# parser.add_argument("--output_true", required=True, help="Output GFF3 file of polyA sites associated with genes having stop codons.")
# parser.add_argument("--output_false", required=True, help="Output GFF3 file of polyA sites associated with genes not having stop codons.")
# parser.add_argument("--output_within_gene", required=True, help="Output GFF3 file of polyA sites which are within the genes")
# parser.add_argument("--output_not_matched", required=True, help="Output GFF3 file of polyA sites not associated with any genes.")
# args = parser.parse_args()


# def main() -> None:
#     db_genes: FeatureDB = gffutils.create_db(args.genes, dbfn=":memory:", merge_strategy="create_unique")
#     db_polyA: FeatureDB = gffutils.create_db(args.polyA, dbfn=":memory:", merge_strategy="create_unique")
#     sequences: Dict[str, str] = load_fasta_sequences(args.fasta)

#     gene_intervals = process_gene_intervals(db_genes)
#     polyA_sites = extract_polyA_sites(db_polyA)

#     matched_polyA, within_gene_polyA_ids, unmatched_polyA_ids = assign_polyA_to_genes(gene_intervals, polyA_sites, db_genes, sequences)

#     # Write output files directly while iterating through polyA features
#     with open(args.output_true, "w") as true_file, open(args.output_false, "w") as false_file, open(args.output_within_gene, "w") as within_file, open(args.output_not_matched, "w") as unmatched_file:
#         for feature in db_polyA.features_of_type("polyA"):
#             if feature.id in {entry[1] for entry in matched_polyA if entry[-1] is True}:
#                 true_file.write(str(feature) + "\n")
#             elif feature.id in {entry[1] for entry in matched_polyA if entry[-1] is False}:
#                 false_file.write(str(feature) + "\n")
#             elif feature.id in within_gene_polyA_ids:
#                 within_file.write(str(feature) + "\n")
#             elif feature.id in unmatched_polyA_ids:
#                 unmatched_file.write(str(feature) + "\n")


# def load_fasta_sequences(fasta_file: str) -> Dict[str, str]:
#     """
#     Load sequences from a FASTA file into a dictionary.
#     """
#     sequences = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}
#     return sequences


# def process_gene_intervals(db) -> List[List[Union[str, int]]]:
#     """
#     Generate intervals between genes in the GFF3 database.
#     """
#     intervals = []
#     genes_by_contig = {}

#     for gene in db.features_of_type("gene"):
#         if gene.seqid not in genes_by_contig:
#             genes_by_contig[gene.seqid] = []
#         genes_by_contig[gene.seqid].append(gene)

#     for contig, genes in genes_by_contig.items():
#         genes.sort(key=lambda g: g.start)

#         # First gene
#         first_gene = genes[0]
#         intervals.append([contig, 0, first_gene.start, None, first_gene.id, None, first_gene.strand])

#         # Gene intervals
#         for i in range(len(genes) - 1):
#             gene_left = genes[i]
#             gene_right = genes[i + 1]
#             intervals.append([contig, gene_left.end, gene_right.start, gene_left.id, gene_right.id, gene_left.strand, gene_right.strand])

#         # Last gene
#         last_gene = genes[-1]
#         intervals.append([contig, last_gene.end, 0, last_gene.id, "", last_gene.strand, ""])

#     return intervals


# def extract_polyA_sites(db) -> List[List[Union[str, int]]]:
#     """
#     Extract polyA site information, including ID, from the GFF3 database.
#     """
#     return [[feature.seqid, feature.start, feature.strand, feature.id] for feature in db.features_of_type("polyA")]


# def has_stop_codon(gene: Feature, sequences: Dict[str, str]) -> bool:
#     """
#     Check if a gene has a stop codon based on its strand.
#     """
#     sequence = sequences[gene.seqid][gene.start - 1:gene.end].upper()

#     if gene.strand == "+":
#         return sequence[-3:] in ["TAA", "TGA", "TAG"]
#     elif gene.strand == "-":
#         return sequence[:3] in ["TTA", "TCA", "CTA"]
#     return False


# def assign_polyA_to_genes(
#     gene_intervals: List[List[Union[str, int]]],
#     polyA_sites: List[List[Union[str, int]]],
#     db_genes: FeatureDB,
#     sequences: Dict[str, str]
# ) -> (List[List[Union[str, str]]], List[str]):
#     """
#     Assign polyA sites to genes based on intervals and categorize them as True or False.
#     Return a list of matched polyA sites and unmatched polyA IDs.
#     """
#     matched_polyA = []
#     within_gene_polyA_ids= []
#     unmatched_polyA_ids = []


#     for polyA in polyA_sites:
#         contig, start, polyA_strand, polyA_id = polyA
#         matching_intervals = [interval for interval in gene_intervals if interval[0] == contig]
#         within_gene= True
#         matched = False

#         for interval in matching_intervals:
#             _, end_left, start_right, gene_left, gene_right, strand_left, strand_right = interval

#             # First half-interval
#             if end_left == 0 and start < start_right + 3:
#                 within_gene= False 
#                 if polyA_strand == "-" and strand_right == "-":
#                     gene = db_genes[gene_right]
#                     matched_polyA.append([contig, polyA_id, gene_right, "-", has_stop_codon(gene, sequences)])
#                     matched = True

#             # Middle intervals
#             if start_right != 0:
#                 if (end_left - 3) < start < start_right:
#                     within_gene= False
#                     if polyA_strand == "+" and strand_left == "+":
#                         gene = db_genes[gene_left]
#                         matched_polyA.append([contig, polyA_id, gene_left, "+", has_stop_codon(gene, sequences)])
#                         matched = True
#                 if end_left < start < (start_right + 3):    
#                     within_gene= False
#                     if polyA_strand == "-" and strand_right == "-":
#                         gene = db_genes[gene_right]
#                         matched_polyA.append([contig, polyA_id, gene_right, "-", has_stop_codon(gene, sequences)])
#                         matched = True

#             # Last half-interval
#             if start_right == 0 and end_left - 3 < start:
#                 within_gene= False
#                 if polyA_strand == "+" and strand_left == "+":
#                     gene = db_genes[gene_left]
#                     matched_polyA.append([contig, polyA_id, gene_left, "+", has_stop_codon(gene, sequences)])
#                     matched = True


#         if not matched and within_gene:
#             within_gene_polyA_ids.append(polyA_id)

#         elif not matched and not within_gene:
#             unmatched_polyA_ids.append(polyA_id)

#     return matched_polyA, within_gene_polyA_ids, unmatched_polyA_ids

# if __name__ == "__main__":
#     main()









# ##python categorize_polyA.py -g1 ST2_orfs.gff3 -g2 polyAsites_updated.gff3 -f ST2_sorted_masked.fasta --output_true new_has_stop_polyA.gff3 --output_false new_non_stop_polyA.gff3 --output_not_matched new_not_matched_polyA.gff3

# #!/usr/bin/env python

# from typing import List, Dict, Union, Tuple
# import argparse
# import gffutils
# from gffutils.feature import Feature
# from Bio import SeqIO

# parser = argparse.ArgumentParser(description="Categorize polyA sites into 'has stop codon', 'not having stop codon', and 'not associated with any genes'.")
# parser.add_argument("-g1", "--genes", required=True, help="Input GFF3 file with gene information.")
# parser.add_argument("-g2", "--polyA", required=True, help="Input GFF3 file with polyA site information.")
# parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file containing genomic sequences.")
# parser.add_argument("--output_true", required=True, help="Output GFF3 file of polyA sites associated with gens having stop codons.")
# parser.add_argument("--output_false", required=True, help="Output GFF3 file of polyA sites associated with gens not having stop codons.")
# parser.add_argument("--output_not_matched", required=True, help="Output GFF3 file of polyA sites not associated with any gens")
# args = parser.parse_args()


# def main():
#     db_genes: FeatureDB = gffutils.create_db(args.genes, dbfn=":memory:", merge_strategy="create_unique")
#     db_polyA: FeatureDB = gffutils.create_db(args.polyA, dbfn=":memory:", merge_strategy="create_unique")
#     sequences = load_fasta_sequences(args.fasta)

#     gene_intervals = process_gene_intervals(db_genes)
#     polyA_sites = extract_polyA_sites(db_polyA)

#     matched_polyA, unmatched_polyA_ids = assign_polyA_to_genes(gene_intervals, polyA_sites, db_genes, sequences)

#     true_ids = {entry[1] for entry in matched_polyA if entry[-1] is True}
#     false_ids = {entry[1] for entry in matched_polyA if entry[-1] is False}

#     true_lines = filter_polyAs_by_ids(args.polyA, true_ids)
#     false_lines = filter_polyAs_by_ids(args.polyA, false_ids)
#     unmatched_lines = filter_polyAs_by_ids(args.polyA, unmatched_polyA_ids)

#     with open(args.output_true, "w") as true_file, open(args.output_false, "w") as false_file, open(args.output_not_matched, "w") as unmatched_file:
#         true_file.writelines(true_lines)
#         false_file.writelines(false_lines)
#         unmatched_file.writelines(unmatched_lines)


# def load_fasta_sequences(fasta_file: str) -> Dict[str, str]:
#     """
#     Load sequences from a FASTA file into a dictionary.
#     """
#     sequences = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}
#     return sequences


# def process_gene_intervals(db: FeatureDB) -> List[List[Union[str, int]]]:
#     """
#     Generate intervals between genes in the GFF3 database.
#     """
#     intervals = []
#     genes_by_contig = {}

#     for gene in db.features_of_type("gene"):
#         if gene.seqid not in genes_by_contig:
#             genes_by_contig[gene.seqid] = []
#         genes_by_contig[gene.seqid].append(gene)

#     for contig, genes in genes_by_contig.items():
#         genes.sort(key=lambda g: g.start)

#         # First gene
#         first_gene = genes[0]
#         intervals.append([contig, 0, first_gene.start, None, first_gene.id, None, first_gene.strand])

#         # Gene intervals
#         for i in range(len(genes) - 1):
#             gene_left = genes[i]
#             gene_right = genes[i + 1]
#             intervals.append([contig, gene_left.end, gene_right.start, gene_left.id, gene_right.id, gene_left.strand, gene_right.strand])

#         # Last gene
#         last_gene = genes[-1]
#         intervals.append([contig, last_gene.end, 0, last_gene.id, "", last_gene.strand, ""])

#     return intervals


# def extract_polyA_sites(db: FeatureDB) -> List[List[Union[str, int]]]:
#     """
#     Extract polyA site information, including ID, from the GFF3 database.
#     """
#     return [[feature.seqid, feature.start, feature.strand, feature.id] for feature in db.features_of_type("polyA")]


# def has_stop_codon(gene: gffutils.Feature, sequences: Dict[str, str]) -> bool:
#     """
#     Check if a gene has a stop codon based on its strand.
#     """
#     sequence = sequences[gene.seqid][gene.start - 1:gene.end].upper()

#     if gene.strand == "+":
#         return sequence[-3:] in ["TAA", "TGA", "TAG"]
#     elif gene.strand == "-":
#         return sequence[:3] in ["TTA", "TCA", "CTA"]
#     return False


# def assign_polyA_to_genes(
#     gene_intervals: List[List[Union[str, int]]],
#     polyA_sites: List[List[Union[str, int]]],
#     db_genes: FeatureDB,
#     sequences: Dict[str, str]
# ) -> (List[List[Union[str, str]]], List[str]):
#     """
#     Assign polyA sites to genes based on intervals and categorize them as True or False.
#     Return a list of matched polyA sites and unmatched polyA IDs.
#     """
#     matched_polyA = []
#     unmatched_polyA_ids = []

#     for polyA in polyA_sites:
#         contig, start, polyA_strand, polyA_id = polyA
#         matching_intervals = [interval for interval in gene_intervals if interval[0] == contig]
#         matched = False

#         for interval in matching_intervals:
#             _, end_left, start_right, gene_left, gene_right, strand_left, strand_right = interval

#             # First half-interval
#             if end_left == 0 and start < start_right + 3 and polyA_strand == "-" and strand_right == "-":
#                 gene = db_genes[gene_right]
#                 matched_polyA.append([contig, polyA_id, gene_right, "-", has_stop_codon(gene, sequences)])
#                 matched = True

#             # Middle intervals
#             elif start_right != 0:
#                 if polyA_strand == "+" and strand_left == "+" and (end_left - 3) < start < start_right:
#                     gene = db_genes[gene_left]
#                     matched_polyA.append([contig, polyA_id, gene_left, "+", has_stop_codon(gene, sequences)])
#                     matched = True
#                 elif polyA_strand == "-" and strand_right == "-" and end_left < start < (start_right + 3):
#                     gene = db_genes[gene_right]
#                     matched_polyA.append([contig, polyA_id, gene_right, "-", has_stop_codon(gene, sequences)])
#                     matched = True

#             # Last half-interval
#             elif start_right == 0 and end_left - 3 < start and polyA_strand == "+" and strand_left == "+":
#                 gene = db_genes[gene_left]
#                 matched_polyA.append([contig, polyA_id, gene_left, "+", has_stop_codon(gene, sequences)])
#                 matched = True

#         if not matched:
#             unmatched_polyA_ids.append(polyA_id)

#     return matched_polyA, unmatched_polyA_ids


# def filter_polyAs_by_ids(gff3_file: str, polyA_ids: List[str]) -> List[str]:
#     """
#     Filter lines from the GFF3 file matching the given polyA IDs.
#     """
#     filtered_lines = []
#     with open(gff3_file, "r") as gff3:
#         for line in gff3:
#             if line.startswith("#"):
#                 continue
#             columns = line.strip().split("\t")
#             attributes = columns[8] if len(columns) > 8 else ""
#             polyA_id = next((attr.split("=")[1] for attr in attributes.split(";") if attr.startswith("ID=")), None)
#             if polyA_id in polyA_ids:
#                 filtered_lines.append(line)
#     return filtered_lines


# if __name__ == "__main__":
#     main()
  







### The script works perefctly!
## $ python categorize_polyA.py -g1 ST2_orfs.gff3 -g2 polyAsites_updated.gff3 -f ST2_sorted_masked.fasta --output_yes has_stop_polyA.gff3 --output_no non_stop_polyA.gff3 --output_not_matched not_matched_polyA.gff3

# import argparse
# import gffutils
# from Bio import SeqIO

# def assign_polyA_to_genes(gene_pairs, polyA_sites: list[list[str, int, str, str]], db_genes, sequences):
#     """
#     Assign polyA sites to genes based on the criteria described, including the polyA ID in the results.
#     Also adds a 'yes' or 'no' field indicating the presence of a stop codon.
#     Tracks unmatched polyA sites.
#     """
#     results2 = []
#     unmatched_polyA_ids = []  # Store unmatched polyA IDs

#     for polyA in polyA_sites:
#         contig, start, polyA_strand, polyA_id = polyA  # Include polyA ID here

#         # Find the matching lists for this contig
#         matching_gene_pairs = [gp for gp in gene_pairs if str(gp[0]) == str(contig)]
#         matched = False  # Flag to track if the polyA was matched

#         for gene_pair in matching_gene_pairs:
#             contig, end_left, start_right, gene_left, gene_right, strand_left, strand_right = gene_pair

#             # Check for the first half-interval
#             if end_left == 0 and start < start_right + 3:
#                 if polyA_strand == "-" and strand_right == "-":
#                     gene = db_genes[gene_right]
#                     stop_codon = has_stop_codon(gene, sequences)
#                     results2.append([contig, polyA_id, gene_right, "-", stop_codon])
#                     matched = True

#             elif start_right != 0:
#                 # For polyA and gene both on "+" strand: shift interval 3 to the left
#                 if polyA_strand == "+" and strand_left == "+":
#                     if (end_left - 3) < start < start_right:  # Shift interval 3 to the left
#                         gene = db_genes[gene_left]
#                         stop_codon = has_stop_codon(gene, sequences)
#                         results2.append([contig, polyA_id, gene_left, "+", stop_codon])
#                         matched = True

#                 # For polyA and gene both on "-" strand: shift interval 3 to the right
#                 elif polyA_strand == "-" and strand_right == "-":
#                     if end_left < start < (start_right + 3):  # Shift interval 3 to the right
#                         gene = db_genes[gene_right]
#                         stop_codon = has_stop_codon(gene, sequences)
#                         results2.append([contig, polyA_id, gene_right, "-", stop_codon])
#                         matched = True

#             # Check for the last half-interval
#             elif start_right == 0 and end_left - 3 < start:
#                 if polyA_strand == "+" and strand_left == "+":
#                     gene = db_genes[gene_left]
#                     stop_codon = has_stop_codon(gene, sequences)
#                     results2.append([contig, polyA_id, gene_left, "+", stop_codon])
#                     matched = True

#         if not matched:
#             unmatched_polyA_ids.append(polyA_id)  # Record unmatched polyA ID

#     return results2, unmatched_polyA_ids

# ##ibelieve the error is rising from ere that it is not considering the reverse comliment of the negative strands
# # def has_stop_codon(gene, sequences):
# #     """
# #     Check if a given gene has a stop codon at the end.
# #     """
# #     sequence = sequences[gene.seqid][gene.start - 1:gene.end].upper()
# #     return "yes" if sequence[-3:] in ["TAA", "TGA", "TAG"] else "no"

# def has_stop_codon(gene, sequences):
#     """
#     Check if a given gene has a stop codon at the end, considering the strand direction.
#     """
#     # Extract the gene sequence
#     sequence = sequences[gene.seqid][gene.start - 1:gene.end].upper()

#     if gene.strand == "+":
#         # Check for stop codon on the positive strand
#         return "yes" if sequence[-3:] in ["TAA", "TGA", "TAG"] else "no"
#     elif gene.strand == "-":
#         return "yes" if sequence[:3] in ["TTA", "TCA", "CTA"] else "no"
        
#     else:
#         # If strand information is unavailable, return "unknown"
#         return "unknown"

# def process_gff3(db):
#     """
#     Process the genes GFF3 database and generate gene pair information.
#     """
#     results = []
#     genes_by_contig = {}

#     for gene in db.features_of_type("gene"):
#         if gene.seqid not in genes_by_contig:
#             genes_by_contig[gene.seqid] = []
#         genes_by_contig[gene.seqid].append(gene)

#     for contig, genes in genes_by_contig.items():
#         genes = sorted(genes, key=lambda g: g.start)

#         first_gene = genes[0]
#         results.append([contig, 0, first_gene.start, None, first_gene.id, None, first_gene.strand])

#         for i in range(len(genes) - 2):
#             gene_left = genes[i]
#             gene_right = genes[i + 1]
#             results.append([contig, gene_left.end, gene_right.start, gene_left.id, gene_right.id, gene_left.strand, gene_right.strand])

#         last_gene = genes[-1]
#         results.append([contig, last_gene.end, 0, last_gene.id, "", last_gene.strand, ""])

#     return results


# def process_polyA_sites(db):
#     """
#     Process the polyA GFF3 database and extract polyA site information, including the polyA ID.
#     """
#     polyA_sites = []
#     for feature in db.features_of_type("polyA"):
#         contig = feature.seqid
#         start = feature.start
#         strand = feature.strand
#         polyA_id = feature.id
#         polyA_sites.append([contig, start, strand, polyA_id])
#     return polyA_sites


# def filter_polyAs_by_ids(gff3_file, polyA_ids):
#     """
#     Extract GFF3 lines matching the given polyA IDs.
#     """
#     filtered_lines = []
#     with open(gff3_file, "r") as gff3:
#         for line in gff3:
#             if line.startswith("#"):
#                 continue
#             columns = line.strip().split("\t")
#             if len(columns) > 8:
#                 attributes = columns[8]
#                 polyA_id = None
#                 for attr in attributes.split(";"):
#                     if attr.startswith("ID="):
#                         polyA_id = attr.split("=")[1]
#                         break
#                 if polyA_id and polyA_id in polyA_ids:
#                     filtered_lines.append(line)
#     return filtered_lines


# def load_fasta_sequences(fasta_file):
#     """
#     Load sequences from a FASTA file into a dictionary.
#     """
#     sequences = {}
#     for record in SeqIO.parse(fasta_file, "fasta"):
#         sequences[record.id] = str(record.seq)
#     return sequences


# def main():
#     parser = argparse.ArgumentParser(description="Categorize polyA sites into 'yes', 'no', and 'not matched'.")
#     parser.add_argument("-g1", "--genes", required=True, help="Input GFF3 file with gene information.")
#     parser.add_argument("-g2", "--polyA", required=True, help="Input GFF3 file with polyA site information.")
#     parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file containing genomic sequences.")
#     parser.add_argument("--output_yes", required=True, help="Output GFF3 file of polyA sites associated with gens having stop codons.")
#     parser.add_argument("--output_no", required=True, help="Output GFF3 file of polyA sites associated with gens not having stop codons.")
#     parser.add_argument("--output_not_matched", required=True, help="Output GFF3 file of polyA sites not associated with any gens")
#     args = parser.parse_args()

#     db_genes = gffutils.create_db(args.genes, dbfn=":memory:", merge_strategy="create_unique", keep_order=True)
#     db_polyA = gffutils.create_db(args.polyA, dbfn=":memory:", merge_strategy="create_unique", keep_order=True)
#     sequences = load_fasta_sequences(args.fasta)

#     gene_pairs = process_gff3(db_genes)
#     polyA_sites = process_polyA_sites(db_polyA)

#     categorized_polyA, unmatched_polyA_ids = assign_polyA_to_genes(gene_pairs, polyA_sites, db_genes, sequences)

#     yes_ids = {entry[1] for entry in categorized_polyA if entry[-1] == "yes"}
#     no_ids = {entry[1] for entry in categorized_polyA if entry[-1] == "no"}

#     yes_lines = filter_polyAs_by_ids(args.polyA, yes_ids)
#     no_lines = filter_polyAs_by_ids(args.polyA, no_ids)
#     unmatched_lines = filter_polyAs_by_ids(args.polyA, unmatched_polyA_ids)

#     with open(args.output_yes, "w") as yes_file, open(args.output_no, "w") as no_file, open(args.output_not_matched, "w") as unmatched_file:
#         yes_file.writelines(yes_lines)
#         no_file.writelines(no_lines)
#         unmatched_file.writelines(unmatched_lines)


# if __name__ == "__main__":
#     main()




##ATTENTION!!!! It is working perfectly for yes and no groups(not considering non-matched!):

# import argparse
# import gffutils
# import pprint
# from Bio import SeqIO

# # Custom PrettyPrinter with increased width
# pprint.pprint = lambda obj, **kwargs: pprint.PrettyPrinter(width=120).pprint(obj)

# parser = argparse.ArgumentParser()
# parser.add_argument(
#     "-g1", "--genes",
#     dest="genes_gff3_file",
#     type=str,
#     required=True,
#     help="Input GFF3 file with gene information."
# )
# parser.add_argument(
#     "-g2", "--polyA",
#     dest="polyA_gff3_file",
#     type=str,
#     required=True,
#     help="Input GFF3 file with polyA site information."
# )
# parser.add_argument(
#     "-f", "--fasta",
#     dest="fasta_file",
#     type=str,
#     required=True,
#     help="Input FASTA file containing genomic sequences."
# )
# parser.add_argument(
#     "--output_yes",
#     type=str,
#     required=True,
#     help="Output GFF3 file for 'yes' polyA entries."
# )
# parser.add_argument(
#     "--output_no",
#     type=str,
#     required=True,
#     help="Output GFF3 file for 'no' polyA entries."
# )


# def load_fasta_sequences(fasta_file):
#     """
#     Load sequences from a FASTA file into a dictionary.
#     """
#     sequences = {}
#     for record in SeqIO.parse(fasta_file, "fasta"):
#         sequences[record.id] = str(record.seq)
#     return sequences


# def process_gff3(db):
#     """
#     Process the genes GFF3 database and generate gene pair information.
#     """
#     results = []
#     genes_by_contig = {}

#     # Group genes by contig
#     for gene in db.features_of_type("gene"):
#         if gene.seqid not in genes_by_contig:
#             genes_by_contig[gene.seqid] = []  # Initialize list for this contig
#         genes_by_contig[gene.seqid].append(gene)

#     # Process each contig
#     for contig, genes in genes_by_contig.items():
#         genes = sorted(genes, key=lambda g: g.start)  # Sort genes by start position

#         # First gene
#         first_gene = genes[0]
#         results.append([
#             contig,
#             0,
#             first_gene.start,
#             None,
#             first_gene.id,
#             None,
#             first_gene.strand,
#         ])

#         # All genes in between first and last gene
#         for i in range(len(genes) - 2):
#             gene_left = genes[i]
#             gene_right = genes[i + 1]

#             results.append([
#                 contig,
#                 gene_left.end,
#                 gene_right.start,
#                 gene_left.id,
#                 gene_right.id,
#                 gene_left.strand,
#                 gene_right.strand
#             ])

#         # Last gene
#         last_gene = genes[-1]
#         results.append([
#             contig,
#             last_gene.end,
#             0,
#             last_gene.id,
#             "",
#             last_gene.strand,
#             ""
#         ])

#     return results


# def process_polyA_sites(db) -> list[list[str, int, str, str]]:
#     """
#     Process the polyA GFF3 database and extract polyA site information, including the polyA ID.
#     """
#     polyA_sites = []
#     for feature in db.features_of_type("polyA"):
#         contig = feature.seqid
#         start = feature.start
#         strand = feature.strand
#         polyA_id = feature.id  # Extract the polyA ID
#         polyA_sites.append([contig, start, strand, polyA_id])
#     return polyA_sites


# def has_stop_codon(gene, sequences):
#     """
#     Check if a given gene has a stop codon at the end.
#     """
#     # Extract the gene sequence
#     sequence = sequences[gene.seqid][gene.start - 1:gene.end] 

#     # Convert to uppercase for consistency
#     sequence = sequence.upper()

#     # Check the last three bases for stop codons
#     return "yes" if sequence[-3:] in ["TAA", "TGA", "TAG"] else "no"


# def assign_polyA_to_genes(gene_pairs, polyA_sites: list[list[str, int, str, str]], db_genes, sequences):
#     """
#     Assign polyA sites to genes based on the criteria described, including the polyA ID in the results.
#     Also adds a 'yes' or 'no' field indicating the presence of a stop codon.
#     """
#     results2 = []

#     for polyA in polyA_sites:
#         contig, start, polyA_strand, polyA_id = polyA  # Include polyA ID here

#         # Find the matching lists for this contig
#         matching_gene_pairs = [gp for gp in gene_pairs if str(gp[0]) == str(contig)]

#         for gene_pair in matching_gene_pairs:
#             contig, end_left, start_right, gene_left, gene_right, strand_left, strand_right = gene_pair

#             # Check for the first half-interval
#             if end_left == 0 and start < start_right + 3:
#                 if polyA_strand == "-" and strand_right == "-":
#                     gene = db_genes[gene_right]
#                     stop_codon = has_stop_codon(gene, sequences)
#                     results2.append([contig, polyA_id, gene_right, "-", stop_codon])

#             elif start_right != 0:
#                 # For polyA and gene both on "+" strand: shift interval 3 to the left
#                 if polyA_strand == "+" and strand_left == "+":
#                     if (end_left - 3) < start < start_right:  # Shift interval 3 to the left
#                         gene = db_genes[gene_left]
#                         stop_codon = has_stop_codon(gene, sequences)
#                         results2.append([contig, polyA_id, gene_left, "+", stop_codon])

#                 # For polyA and gene both on "-" strand: shift interval 3 to the right
#                 elif polyA_strand == "-" and strand_right == "-":
#                     if end_left < start < (start_right + 3):  # Shift interval 3 to the right
#                         gene = db_genes[gene_right]
#                         stop_codon = has_stop_codon(gene, sequences)
#                         results2.append([contig, polyA_id, gene_right, "-", stop_codon])

#             # Check for the last half-interval
#             elif start_right == 0 and end_left - 3 < start:
#                 if polyA_strand == "+" and strand_left == "+":
#                     gene = db_genes[gene_left]
#                     stop_codon = has_stop_codon(gene, sequences)
#                     results2.append([contig, polyA_id, gene_left, "+", stop_codon])

#     return results2


# def filter_polyAs_by_ids(gff3_file, gene_polyA_list, output_yes, output_no):
#     """
#     Write polyA entries from the GFF3 file into separate 'yes' and 'no' output GFF3 files
#     based on their IDs in the gene_polyA_list.
#     """
#     yes_ids = {entry[1] for entry in gene_polyA_list if entry[-1] == "yes"}
#     no_ids = {entry[1] for entry in gene_polyA_list if entry[-1] == "no"}

#     with open(gff3_file, "r") as gff3, open(output_yes, "w") as yes_file, open(output_no, "w") as no_file:
#         for line in gff3:
#             if line.startswith("#"):  # Write header lines
#                 yes_file.write(line)
#                 no_file.write(line)
#                 continue

#             columns = line.strip().split("\t")
#             if len(columns) > 8:  # Ensure the line has attributes
#                 attributes = columns[8]
#                 polyA_id = None
#                 for attr in attributes.split(";"):
#                     if attr.startswith("ID="):  # Extract the ID attribute
#                         polyA_id = attr.split("=")[1]
#                         break
#                 if polyA_id in yes_ids:
#                     yes_file.write(line)
#                 elif polyA_id in no_ids:
#                     no_file.write(line)


# def main(args):
#     # Create a database for the genes GFF3 file
#     db_genes = gffutils.create_db(
#         args.genes_gff3_file,
#         dbfn=":memory:",
#         merge_strategy="create_unique",
#         keep_order=True,
#     )
#     # Create a database for the polyA sites GFF3 file
#     db_polyA = gffutils.create_db(
#         args.polyA_gff3_file,
#         dbfn=":memory:",
#         merge_strategy="create_unique",
#         keep_order=True,
#     )

#     # Load sequences from FASTA
#     sequences = load_fasta_sequences(args.fasta_file)

#     # Generate gene pairs from the genes GFF3 file
#     result_list_gene_pairs = process_gff3(db_genes)

#     # Extract polyA sites from the polyA sites GFF3 file
#     polyA_sites = process_polyA_sites(db_polyA)

#     # Assign polyA sites to genes
#     gene_polyA_list = assign_polyA_to_genes(result_list_gene_pairs, polyA_sites, db_genes, sequences)

#     # Filter and write 'yes' and 'no' polyA entries to separate GFF3 files
#     filter_polyAs_by_ids(args.polyA_gff3_file, gene_polyA_list, args.output_yes, args.output_no)


# if __name__ == "__main__":
#     args = parser.parse_args()
#     main(args)







