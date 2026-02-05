#!/usr/bin/env python

# $ python UTR_length_finder.py -g1 ST2_orfs.gff3 -g2 march4_has_stop_polyA.gff3 -f ST2_sorted_masked.fasta > UTR_length_info_last_version.txt

"""
Estimate 3′ UTR lengths by assigning poly(A) sites (pre-filtered to genes with stop codons)
to their flanking genes and measuring the distance from the gene boundary to the poly(A) site.

This script:
1) Builds in-memory gffutils databases for a gene annotation GFF3 and a polyA-site GFF3
   (where polyA sites are already associated with genes that have stop codons).
2) Loads the reference genome FASTA with SeqIO.index for coordinate-based sequence access.
3) Constructs non-overlapping “big gene” intervals per contig (overlapping smaller genes are
   skipped in favor of the longest spanning gene) and defines intergenic intervals between them.
4) Assigns each poly(A) site to the appropriate flanking gene when the site falls inside an
   interval (with a ±3 bp tolerance) and the strands match.
5) Computes the inferred 3′ UTR length as the distance between the poly(A) site and the gene-end
   boundary implied by the interval:
     - '+' strand:  UTR = polyA.start - interval.start + 1
     - '-' strand:  UTR = interval.end - polyA.start + 1
6) Saves a histogram of UTR lengths (filtered to < 1000 bp) as 'UTR_length_distribution.jpg'.

QC / diagnostics:
- Prints (gene_id, polyA_id, UTR_length) for each assigned site.
- Warns if a UTR length is unusually long (≥ 1000 bp), suggesting a potential missing annotation.
- Collects polyA entries with UTR length in {-3, -2, -1, 0} into 'stop_polyA_align.gff3' to help
  quantify polyA sites that align on/near the stop-codon boundary.
- Warns if the same gene is encountered more than once with UTR lengths differing by > 50 bp.

Example:
  python UTR_length_finder.py -g1 genes.gff3 -g2 has_stop_polyA.gff3 -f genome.fasta > UTR_length_info.txt

Outputs:
  UTR_length_distribution.jpg
  stop_polyA_align.gff3

Dependencies:
  gffutils, Biopython, pandas, matplotlib, seaborn

Author: Tania Iranpour
"""

import sys
import argparse

import gffutils
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


from typing import List, Dict, Optional
from gffutils.feature import Feature
from gffutils.interface import FeatureDB
from Bio.SeqRecord import SeqRecord


# Argument Parser
parser = argparse.ArgumentParser(description="Find UTR length of genes containing stop codon.")
parser.add_argument("-g1", "--genes", required=True, help="Input GFF3 file with gene annotation.")
parser.add_argument("-g2", "--polyA", required=True, help="Input GFF3 file with polyA sites which are associated with genes having stop codons.")
parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file containing genomic sequences.")
args = parser.parse_args()


class Interval:
    """Class representing an interval between two genes."""
    def __init__(self, contig: str, gene_left: Optional[Feature], gene_right: Optional[Feature]):
        self.contig = contig
        self.start = gene_left.end+1 if gene_left else 0
        self.end = gene_right.start-1 if gene_right else sys.maxsize
        self.gene_left = gene_left
        self.gene_right = gene_right
        self.strand_left = gene_left.strand if gene_left else None
        self.strand_right = gene_right.strand if gene_right else None


def main() -> None:
    # FIX: Use in-memory DB cautiously, large files might cause issues.
    db_genes: FeatureDB = gffutils.create_db(args.genes, dbfn=':memory:', merge_strategy="create_unique")
    db_polyA: FeatureDB = gffutils.create_db(args.polyA, dbfn=':memory:', merge_strategy="create_unique")
    sequences = SeqIO.index(args.fasta, 'fasta')

    # Generate intervals
    # Explicitly declare the types of the returned variables
    gene_intervals: List[Interval] = process_gene_intervals(db_genes)
    UTR_lengths_list: List[int] = find_UTR_lengths(gene_intervals, db_polyA, sequences)

    UTR_lengths_series = pd.Series(UTR_lengths_list, name="UTR Length")

    plot = sns.histplot(
        UTR_lengths_series,
        bins=1000,
        kde=False
    )
    plot.set_xlim(-50, 400)
    plot.get_figure().savefig("UTR_length_distribution.jpg", dpi=300)



def process_gene_intervals(db: FeatureDB) -> List[Interval]:
    """
    Build intervals by:
      1) Iterating through each seqid in db.seqids().
      2) Retrieving genes for that seqid using db.region(...).
      3) Sorting them by start coordinate; if the starts are the same, prioritize descending ends.
      4) Skipping overshadowed genes by keeping the largest(big gene) (non-overlapping).
      5) Creating intervals between these final big genes, plus intervals
    """
    intervals: List[Interval] = []

    for contig in db.seqids():
        # Get all 'gene' features on this contig
        genes = list(db.region(seqid=contig, featuretype="gene", completely_within=False))

        if not genes:
            continue

        # Sort by ascending start and descending ends in the second priority
        genes.sort(key=lambda g: (g.start, -g.end))

        # Keep only non-overlapping 'big' genes
        big_genes: List[Feature] = []
        current = genes[0]
        for i in range(1, len(genes)):
            g = genes[i]
            # Overlap check
            if g.start <= current.end:
                # They overlap => keep the one that extends 
                if g.end > current.end:
                    current = g
            else:
                # No overlap => finalize 'current'
                big_genes.append(current)
                current = g

        # Add the last 'current'
        big_genes.append(current)

        # half-interval from "start-of-seqid" to the first big gene
        first_gene = big_genes[0]
        intervals.append(
            Interval(contig=contig, gene_left=None, gene_right=first_gene)
        )

        # Intervals between big genes
        for i in range(len(big_genes) - 1):
            left_gene = big_genes[i]
            right_gene = big_genes[i + 1]
            intervals.append(
                Interval(contig=contig, gene_left=left_gene, gene_right=right_gene)
            )

        # half-interval from the last gene to "end-of-seqid"
        last_gene = big_genes[-1]
        intervals.append(
            Interval(contig=contig, gene_left=last_gene, gene_right=None)
        )

    return intervals


def find_UTR_lengths(
    gene_intervals: List[Interval],
    db_polyA: FeatureDB,
    sequences: Dict[str, SeqRecord],
) -> List[int]:
    """
    Assign polyA sites to genes, categorize them, and return UTR lengths.
    Additionally, if the UTR length is in [-3, -2, -1, 0], write those polyA lines
    to 'stop_polyA_align.gff3'. Also print a warning if the same gene is encountered
    twice with a difference in UTR lengths > 50 bp.
    """
    UTR_lengths = []
    align_polyAs = []  # will collect lines for polyAs with UTR length in [-3, -2, -1, 0]

    # These variables track the last gene and last UTR length we encountered
    last_gene = None
    last_utr_length = None
    last_polyA_id = None

    with open("stop_polyA_align.gff3", "w") as align_file:
        for polyA in db_polyA.features_of_type("polyA"):
            for interval in gene_intervals:
                if polyA.seqid != interval.contig:
                    continue

                # Check if polyA falls within an interval (±3 bp) and if strand matches
                if interval.start - 3 <= polyA.start <= interval.end + 3:
                    if polyA.strand == "-" and interval.strand_right == "-":
                        current_gene = interval.gene_right.id
                        current_UTR_length = interval.end - polyA.start + 1
                    elif polyA.strand == "+" and interval.strand_left == "+":
                        current_gene = interval.gene_left.id
                        current_UTR_length = polyA.start - interval.start + 1
                    else:
                        # If the strand doesn't match, skip
                        continue

                    # Save the UTR length if it's < 1000
                    if current_UTR_length < 1000:
                        UTR_lengths.append(current_UTR_length)
                    else:
                        print(f"WARNING: UTR ({current_UTR_length} bp) for gene {current_gene} (against {polyA.id}) is unusually long; a missing gene may be possible.")

                    print(current_gene, polyA.id, current_UTR_length)

                    # If the length is in [-3, -2, -1, 0], record the line in stop_polyA_align.gff3 
                    #(want to check what ratio of the has stop codon genes also have a polyA which is aligining on the stop codon)
                    if current_UTR_length in [-3, -2, -1, 0]:
                        align_polyAs.append(str(polyA))

                    # Check if we are hitting the same gene again
                    if last_gene == current_gene:
                        # If so, see if the difference is > 50
                        if abs(last_utr_length - current_UTR_length) > 50:
                            print(f"WARNING: For gene {current_gene}, polyA {last_polyA_id} and {polyA.id} differ in UTR length by more than 50 bp. ")

                    # Update the tracking variables for the next iteration
                    last_gene = current_gene
                    last_utr_length = current_UTR_length
                    last_polyA_id = polyA.id

                    # Once assigned, we break out of the intervals loop
                    break

        # After collecting them, write them all at once
        for line in align_polyAs:
            align_file.write(line + "\n")

    return UTR_lengths


if __name__ == "__main__":
    main()



# # $ python UTR_length_finder.py -g1 ST2_orfs.gff3 -g2 march4_has_stop_polyA.gff3 -f ST2_sorted_masked.fasta > UTR_length_info_new.txt

# # $ wc -l stop_polyA_align.gff3
# # 410 stop_polyA_align.gff3

# #!/usr/bin/env python
# import sys
# import argparse

# import gffutils
# from Bio import SeqIO
# import matplotlib.pyplot as plt
# import seaborn as sns
# import pandas as pd


# from typing import List, Dict, Optional
# from gffutils.feature import Feature
# from gffutils.interface import FeatureDB
# from Bio.SeqRecord import SeqRecord


# # Argument Parser
# parser = argparse.ArgumentParser(description="Find UTR length of genes containing stop codon.")
# parser.add_argument("-g1", "--genes", required=True, help="Input GFF3 file with gene information.")
# parser.add_argument("-g2", "--polyA", required=True, help="Input GFF3 file with polyA sites which are associated with genes having stop codons.")
# parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file containing genomic sequences.")
# args = parser.parse_args()


# class Interval:
#     """Class representing an interval between two genes."""
#     def __init__(self, contig: str, gene_left: Optional[Feature], gene_right: Optional[Feature]):
#         self.contig = contig
#         self.start = gene_left.end+1 if gene_left else 0
#         self.end = gene_right.start-1 if gene_right else sys.maxsize
#         self.gene_left = gene_left
#         self.gene_right = gene_right
#         self.strand_left = gene_left.strand if gene_left else None
#         self.strand_right = gene_right.strand if gene_right else None


# def main() -> None:
#     # FIX: Use in-memory DB cautiously, large files might cause issues.
#     db_genes: FeatureDB = gffutils.create_db(args.genes, dbfn='memory', merge_strategy="create_unique")
#     db_polyA: FeatureDB = gffutils.create_db(args.polyA, dbfn='memory', merge_strategy="create_unique")
#     sequences = SeqIO.index(args.fasta, 'fasta')

#     # Generate intervals
#     # Explicitly declare the types of the returned variables
#     gene_intervals: List[Interval] = process_gene_intervals(db_genes)
#     UTR_lengths_list: List[int] = find_UTR_lengths(gene_intervals, db_polyA, sequences)


#     UTR_lengths_series = pd.Series(UTR_lengths_list, name="UTR Length")
#     sns.histplot(
#         UTR_lengths_series,
#         bins=1000,
#         kde=True
#     ).get_figure().savefig("UTR_length_distribution.jpg", dpi=300)

      

# def process_gene_intervals(db: FeatureDB) -> List[Interval]:
#     """
#     Generate intervals between 'big' genes in the GFF3 database,
#     skipping smaller genes overlapped by a bigger gene.
#     """
#     intervals: List[Interval] = []
#     genes_by_contig: Dict[str, List[Feature]] = {}

#     # Collect genes by contig
#     for gene in db.features_of_type("gene"):
#         if gene.seqid not in genes_by_contig:
#             genes_by_contig[gene.seqid] = []
#         genes_by_contig[gene.seqid].append(gene)

#     # Process each contig individually
#     for contig, genes in genes_by_contig.items():
        
#         # Sort genes by their start positions and in case of starting at the same point prioritise the ones with further ends
#         genes.sort(key=lambda g: (g.start, -g.end))

#         # Create interval from "start of genome" to the first gene
#         first_gene = genes[0]
#         intervals.append(Interval(contig=contig, gene_left=None, gene_right=first_gene))
#         correct_start = first_gene

#         # Iterate over the remaining genes
#         for i in range(1, len(genes)):
            
#             # If next gene is overshadowed by correct_start (start < correct_start.end), skip it
#             if genes[i].start < correct_start.end:
#                 continue

#             # Otherwise, create an interval from the last big gene to this new gene
#             intervals.append(Interval(contig=contig, gene_left=correct_start, gene_right=genes[i]))

#             correct_start = genes[i]

#         # Create interval from the last big gene to "end of genome"
#         intervals.append(Interval(contig=contig, gene_left=correct_start, gene_right=None))

#     return intervals

# def find_UTR_lengths(
#     gene_intervals: List[Interval],
#     db_polyA: FeatureDB,
#     sequences: Dict[str, SeqRecord],
# ) -> List[int]:
#     """
#     Assign polyA sites to genes, categorize them, and return UTR lengths.
#     Additionally, if the UTR length is in [-3, -2, -1, 0], write those polyA lines
#     to 'stop_polyA_align.gff3'. Also print a warning if the same gene is encountered
#     twice with a difference in UTR lengths > 50 bp.
#     """
#     UTR_lengths = []
#     align_polyAs = []  # will collect lines for polyAs with UTR length in [-3, -2, -1, 0]

#     # These variables track the last gene and last UTR length we encountered
#     last_gene = None
#     last_utr_length = None
#     last_polyA_id = None

#     with open("stop_polyA_align.gff3", "w") as align_file:
#         for polyA in db_polyA.features_of_type("polyA"):
#             for interval in gene_intervals:
#                 if polyA.seqid != interval.contig:
#                     continue

#                 # Check if polyA falls within an interval (±3 bp) and if strand matches
#                 if interval.start - 3 <= polyA.start <= interval.end + 3:
#                     if polyA.strand == "-" and interval.strand_right == "-":
#                         current_gene = interval.gene_right.id
#                         current_UTR_length = interval.end - polyA.start + 1
#                     elif polyA.strand == "+" and interval.strand_left == "+":
#                         current_gene = interval.gene_left.id
#                         current_UTR_length = polyA.start - interval.start + 1
#                     else:
#                         # If the strand doesn't match, skip
#                         continue

#                     # Save the UTR length if it's < 1000
#                     if current_UTR_length < 1000:
#                         UTR_lengths.append(current_UTR_length)
#                     else:
#                         print(f"WARNING: UTR ({current_UTR_length} bp) for gene {current_gene} is unusually long; a missing gene may be possible.")

#                     print(current_gene, polyA.id, current_UTR_length)

#                     # If the length is in [-3, -2, -1, 0], record the line in stop_polyA_align.gff3 
#                     #(want to check what ratio of the has stop codon genes also have a polyA which is aligining on the stop codon)
#                     if current_UTR_length in [-3, -2, -1, 0]:
#                         align_polyAs.append(str(polyA))

#                     # Check if we are hitting the same gene again
#                     if last_gene == current_gene:
#                         # If so, see if the difference is > 50
#                         if abs(last_utr_length - current_UTR_length) > 50:
#                             print(f"WARNING: For gene {current_gene}, polyA {last_polyA_id} and {polyA.id} differ in UTR length by more than 50 bp. ")

#                     # Update the tracking variables for the next iteration
#                     last_gene = current_gene
#                     last_utr_length = current_UTR_length
#                     last_polyA_id = polyA.id

#                     # Once assigned, we break out of the intervals loop
#                     break

#         # After collecting them, write them all at once
#         for line in align_polyAs:
#             align_file.write(line + "\n")

#     return UTR_lengths


# if __name__ == "__main__":
#     main()




# # $ python UTR_length_finder.py -g1 ST2_orfs.gff3 -g2 march4_has_stop_polyA.gff3 -f ST2_sorted_masked.fasta > UTR_length_info_new.txt

# # $ grep -c "there must be" UTR_length_info_new.txt
# # 52 


# # $ python UTR_length_finder.py -g1 ST2_orfs.gff3 -g2 march4_has_stop_polyA.gff3 -f ST2_sorted_masked.fasta > UTR_length_info_new.txt
# #!/usr/bin/env python
# import sys
# import argparse

# import gffutils
# from Bio import SeqIO
# import matplotlib.pyplot as plt
# import seaborn as sns
# import pandas as pd


# from typing import List, Dict, Optional
# from gffutils.feature import Feature
# from gffutils.interface import FeatureDB
# from Bio.SeqRecord import SeqRecord


# # Argument Parser
# parser = argparse.ArgumentParser(description="Categorize polyA sites into genomic contexts.")
# parser.add_argument("-g1", "--genes", required=True, help="Input GFF3 file with gene information.")
# parser.add_argument("-g2", "--polyA", required=True, help="Input GFF3 file with polyA sites which are associated with genes having stop codons.")
# parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file containing genomic sequences.")
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
#     # FIX: Use in-memory DB cautiously, large files might cause issues.
#     db_genes: FeatureDB = gffutils.create_db(args.genes, dbfn="genes.db", force=True, merge_strategy="create_unique")
#     db_polyA: FeatureDB = gffutils.create_db(args.polyA, dbfn="polyA.db", force=True, merge_strategy="create_unique")
#     sequences = SeqIO.index(args.fasta, 'fasta')

#     # Generate intervals
#     gene_intervals = process_gene_intervals(db_genes)

#     UTR_lengths_list = find_UTR_lengths(gene_intervals, db_polyA, sequences)

#     df = pd.DataFrame(UTR_lengths_list, columns=["UTR Length"])
#     plt.xlim(-10, 200)  
#     sns.histplot(df, x="UTR Length", bins=1000, kde=True).get_figure().savefig("UTR_length_distribution.jpg", dpi=300)
      

# def process_gene_intervals(db: FeatureDB) -> List[Interval]:
#     """
#     Generate intervals between 'big' genes in the GFF3 database,
#     skipping smaller genes overlapped by a bigger gene.
#     """
#     intervals: List[Interval] = []
#     genes_by_contig: Dict[str, List[Feature]] = {}

#     # Collect genes by contig
#     for gene in db.features_of_type("gene"):
#         if gene.seqid not in genes_by_contig:
#             genes_by_contig[gene.seqid] = []
#         genes_by_contig[gene.seqid].append(gene)

#     # Process each contig individually
#     for contig, genes in genes_by_contig.items():
#         # Sort genes by their start positions and in case of starting at the same point prioritise the ones with further ends
#         genes.sort(key=lambda g: (g.start, -g.end))

#         # Create interval from "start of genome" to the first gene
#         first_gene = genes[0]
#         intervals.append(Interval(contig=contig, gene_left=None, gene_right=first_gene))
#         correct_start = first_gene

#         # Iterate over the remaining genes
#         for i in range(1, len(genes)):
#             # If next gene is overshadowed by correct_start (start < correct_start.end), skip it
#             if genes[i].start < correct_start.end:
#                 continue

#             # Otherwise, create an interval from the last big gene to this new gene
#             intervals.append(Interval(contig=contig, gene_left=correct_start, gene_right=genes[i]))

#             correct_start = genes[i]

#         # Create interval from the last big gene to "end of genome"
#         intervals.append(Interval(contig=contig, gene_left=correct_start, gene_right=None))

#     return intervals


# def has_stop_codon(gene: Feature, sequences: Dict[str, SeqRecord]) -> bool:
#     """Check if a gene has a stop codon based on its strand."""
#     sequence = str(sequences[gene.seqid].seq[gene.start - 1: gene.end]).upper()
#     if gene.strand == "+":
#         return sequence[-3:] in ["TAA", "TGA", "TAG"]
#     if gene.strand == "-":
#         return sequence[:3] in ["TTA", "TCA", "CTA"]
#     return False


# def find_UTR_lengths(
#     gene_intervals: List[Interval],
#     db_polyA: FeatureDB,
#     sequences: Dict[str, SeqRecord],
# ) -> List[int]: 
#     """Assign polyA sites to genes, categorize them, and return UTR lengths."""
#     UTR_lengths = []
#     geneIdToCheck = None
#     UTRlengthToCheck = None
#     for polyA in db_polyA.features_of_type('polyA'):
#         for interval in gene_intervals:
#             if polyA.seqid != interval.contig:
#                 continue

#             # Check if polyA falls within an interval and the strands match
#             if interval.start - 3 <= polyA.start <= interval.end + 3:
#                 if polyA.strand == "-" and interval.strand_right == "-":
#                     current_gene = interval.gene_right.id
#                     current_UTR_length = interval.end - polyA.start


#                 elif polyA.strand == "+" and interval.strand_left == "+":
#                     current_gene = interval.gene_left.id
#                     current_UTR_length =polyA.start - interval.start


#                 if current_UTR_length < 1000:
#                     UTR_lengths.append(current_UTR_length)
#                 else:
#                     print(f"WARNING: UTR ({current_UTR_length} bp) for gene {current_gene} is unusually long; a missing gene may be possible.")

#                 print(current_gene, polyA.id, current_UTR_length)

#                 # Only do the comparison if both geneIdToCheck and UTRlengthToCheck are already defined
#                 if geneIdToCheck is not None and UTRlengthToCheck is not None:
#                 # Check if we have the same gene, and if the difference is bigger than 50
#                     if geneIdToCheck == current_gene and abs(UTRlengthToCheck - current_UTR_length) > 50:
#                         print(f"WARNING: there might be a gene missing between {polyAtoCheck} and {polyA.id}")

#                 # Update these variables for the next comparison
#                 geneIdToCheck = current_gene
#                 UTRlengthToCheck = current_UTR_length
#                 polyAtoCheck = polyA.id

#     return UTR_lengths


# if __name__ == "__main__":
#     main()
