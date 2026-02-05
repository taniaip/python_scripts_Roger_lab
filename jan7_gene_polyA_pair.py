#!/usr/bin/env python3

"""
Assign poly(A) sites to flanking genes using a simple interval model, and report
the matched gene ID and stop-codon status for each assigned site.

This script builds in-memory gffutils databases from:
  - a gene annotation GFF3 (features_of_type("gene"))
  - a polyA site GFF3 (features_of_type("polyA"))
and loads the reference genome FASTA into a dictionary for quick slicing.

Workflow:
1) Group genes by contig and create “gene-pair” intervals (left gene → right gene)
   spanning intergenic regions, plus two boundary intervals (start→first gene,
   last gene→end).
2) Extract polyA sites (contig, start, strand, polyA ID).
3) For each polyA site, find the corresponding interval on the same contig and
   assign it to the upstream (+ strand) or downstream (- strand) flanking gene
   when the gene and polyA site share the same strand. A ±3 nt tolerance is
   applied at interval boundaries to account for polyA sites near gene ends.
4) For each assigned gene, check whether the gene ends in a canonical stop codon
   (TAA/TGA/TAG) based on the reference FASTA sequence.
5) Print a list of assignments as:
   [contig, polyA_id, associated_gene_id, strand, stop_codon_yes/no]

Note:
- This script is primarily for inspection/debugging: it prints assignments via
  pretty-print rather than writing output files.

Inputs:
  --genes  Gene annotation GFF3
  --polyA  PolyA site GFF3
  --fasta  Reference genome FASTA

Dependencies:
  gffutils, Biopython

Author: Tania Iranpour
"""

import argparse
import gffutils
import pprint
from Bio import SeqIO

# Custom PrettyPrinter with increased width
pprint.pprint = lambda obj, **kwargs: pprint.PrettyPrinter(width=120).pprint(obj)

parser = argparse.ArgumentParser()
parser.add_argument(
    "-g1", "--genes",
    dest="genes_gff3_file",
    type=str,
    required=True,
    help="Input GFF3 file with gene information."
)
parser.add_argument(
    "-g2", "--polyA",
    dest="polyA_gff3_file",
    type=str,
    required=True,
    help="Input GFF3 file with polyA site information."
)
parser.add_argument(
    "-f", "--fasta",
    dest="fasta_file",
    type=str,
    required=True,
    help="Input FASTA file containing genomic sequences."
)


def load_fasta_sequences(fasta_file):
    """
    Load sequences from a FASTA file into a dictionary.
    """
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences


def process_gff3(db):
    """
    Process the genes GFF3 database and generate gene pair information.
    """
    results = []
    genes_by_contig = {}

    # Group genes by contig
    for gene in db.features_of_type("gene"):
        if gene.seqid not in genes_by_contig:
            genes_by_contig[gene.seqid] = []  # Initialize list for this contig
        genes_by_contig[gene.seqid].append(gene)

    # Process each contig
    for contig, genes in genes_by_contig.items():
        genes = sorted(genes, key=lambda g: g.start)  # Sort genes by start position

        # First gene
        first_gene = genes[0]
        results.append([
            contig,
            0,
            first_gene.start,
            None,
            first_gene.id,
            None,
            first_gene.strand,
        ])

        # All genes in between first and last gene
        for i in range(len(genes) - 2):
            gene_left = genes[i]
            gene_right = genes[i + 1]

            results.append([
                contig,
                gene_left.end,
                gene_right.start,
                gene_left.id,
                gene_right.id,
                gene_left.strand,
                gene_right.strand
            ])

        # Last gene
        last_gene = genes[-1]
        results.append([
            contig,
            last_gene.end,
            0,
            last_gene.id,
            "",
            last_gene.strand,
            ""
        ])

    return results


def process_polyA_sites(db) -> list[list[str, int, str, str]]:
    """
    Process the polyA GFF3 database and extract polyA site information, including the polyA ID.
    """
    polyA_sites = []
    for feature in db.features_of_type("polyA"):
        contig = feature.seqid
        start = feature.start
        strand = feature.strand
        polyA_id = feature.id  # Extract the polyA ID
        polyA_sites.append([contig, start, strand, polyA_id])
    return polyA_sites


def has_stop_codon(gene, sequences):
    """
    Check if a given gene has a stop codon at the end.
    """
    # Extract the gene sequence
    sequence = sequences[gene.seqid][gene.start - 1:gene.end] 
    
    # Convert to uppercase for consistency
    sequence = sequence.upper()
    
    # Check the last three bases for stop codons
    return "yes" if sequence[-3:] in ["TAA", "TGA", "TAG"] else "no"



def assign_polyA_to_genes(gene_pairs, polyA_sites: list[list[str, int, str, str]], db_genes, sequences):
    """
    Assign polyA sites to genes based on the criteria described, including the polyA ID in the results.
    Also adds a 'yes' or 'no' field indicating the presence of a stop codon.
    """
    results2 = []

    for polyA in polyA_sites:
        contig, start, polyA_strand, polyA_id = polyA  # Include polyA ID here

        # Find the matching lists for this contig
        matching_gene_pairs = [gp for gp in gene_pairs if str(gp[0]) == str(contig)]

        for gene_pair in matching_gene_pairs:
            contig, end_left, start_right, gene_left, gene_right, strand_left, strand_right = gene_pair

            # Check for the first half-interval
            if end_left == 0 and start < start_right+3:
                if polyA_strand == "-" and strand_right == "-":
                    gene = db_genes[gene_right]
                    stop_codon = has_stop_codon(gene, sequences)
                    results2.append([contig, polyA_id, gene_right, "-", stop_codon])

            # # Check if polyA start site falls into the interval
            # # Not accounting for the ones who are located in the last three bases of the gene!(the new one does!)
            # if start_right != 0 and end_left < start < start_right:
            #     if polyA_strand == "+" and strand_left == "+":
            #         gene = db_genes[gene_left]
            #         stop_codon = has_stop_codon(gene, sequences)
            #         results2.append([contig, polyA_id, gene_left, "+", stop_codon])
            #     elif polyA_strand == "-" and strand_right == "-":
            #         gene = db_genes[gene_right]
            #         stop_codon = has_stop_codon(gene, sequences)
            #         results2.append([contig, polyA_id, gene_right, "-", stop_codon])


            # Check if polyA start site falls into the interval
            elif start_right != 0:
            # For polyA and gene both on "+" strand: shift interval 3 to the left
                if polyA_strand == "+" and strand_left == "+":
                    if (end_left - 3) < start < start_right:  # Shift interval 3 to the left
                        gene = db_genes[gene_left]
                        stop_codon = has_stop_codon(gene, sequences)
                        results2.append([contig, polyA_id, gene_left, "+", stop_codon])

                    # For polyA and gene both on "-" strand: shift interval 3 to the right
                elif polyA_strand == "-" and strand_right == "-":
                    if end_left < start < (start_right + 3):  # Shift interval 3 to the right
                        gene = db_genes[gene_right]
                        stop_codon = has_stop_codon(gene, sequences)
                        results2.append([contig, polyA_id, gene_right, "-", stop_codon])

            # Check for the last half-interval
            elif start_right == 0 and end_left-3 < start:
                if polyA_strand == "+" and strand_left == "+":
                    gene = db_genes[gene_left]
                    stop_codon = has_stop_codon(gene, sequences)
                    results2.append([contig, polyA_id, gene_left, "+", stop_codon])

    return results2


def main(args):
    # Create a database for the genes GFF3 file
    db_genes = gffutils.create_db(
        args.genes_gff3_file,
        dbfn=":memory:",
        merge_strategy="create_unique",
        keep_order=True,
    )
    # Create a database for the polyA sites GFF3 file
    db_polyA = gffutils.create_db(
        args.polyA_gff3_file,
        dbfn=":memory:",
        merge_strategy="create_unique",
        keep_order=True,
    )

    # Load sequences from FASTA
    sequences = load_fasta_sequences(args.fasta_file)

    # Generate gene pairs from the genes GFF3 file
    result_list_gene_pairs = process_gff3(db_genes)

    # Extract polyA sites from the polyA sites GFF3 file
    polyA_sites = process_polyA_sites(db_polyA)

    # Assign polyA sites to genes
    matched_polyA_genes = assign_polyA_to_genes(result_list_gene_pairs, polyA_sites, db_genes, sequences)

    # Pretty-print the results
    pprint.pprint(matched_polyA_genes)


if __name__ == "__main__":
    args = parser.parse_args()
    main(args)




# import argparse
# import gffutils
# import pprint

# # Custom PrettyPrinter with increased width (default is 80) 
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

#         # first gene
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

#         # all genes in between first and last gene
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

#         # last gene
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
#         polyA_id = feature.id  
#         polyA_sites.append([contig, start, strand, polyA_id])
#     return polyA_sites


# def assign_polyA_to_genes(gene_pairs, polyA_sites: list[list[str, int, str, str]]):
#     """
#     Assign polyA sites to genes based on the criteria described, including the polyA ID in the results.
#     """
#     results2 = []

#     for polyA in polyA_sites:
#         contig, start, polyA_strand, polyA_id = polyA  

#         # Find the matching lists for this contig
#         # Find gene pairs of the contig of the current polyA site
#         matching_gene_pairs = [gp for gp in gene_pairs if str(gp[0]) == str(contig)]

#         for gene_pair in matching_gene_pairs:
#             contig, end_left, start_right, gene_left, gene_right, strand_left, strand_right = gene_pair
            
#             # Check for the first half-interval    
#             if end_left == 0 and start < start_right:
#                 if polyA_strand == "-" and strand_right == "-":
#                     results2.append([contig, polyA_id, gene_right, "-"])  
               
#             # Check if polyA start site falls into the interval
#             elif start_right != 0 and end_left < start < start_right:
#                 if polyA_strand == "+" and strand_left == "+":
#                     results2.append([contig, polyA_id, gene_left, "+"]) 
#                 elif polyA_strand == "-" and strand_right == "-":
#                     results2.append([contig, polyA_id, gene_right, "-"]) 

#             # Check for the last half-interval         
#             elif start_right == 0 and end_left < start:
#                 if polyA_strand == "+" and strand_left == "+":
#                     results2.append([contig, polyA_id, gene_left, "+"])        

#     return results2


# def main(args):
#     # Create a database for the genes GFF3 file
#     db_genes = gffutils.create_db(
#         args.genes_gff3_file,
#         dbfn=":memory:",
#         merge_strategy="create_unique"
#     )
#     # Create a database for the polyA sites GFF3 file
#     db_polyA = gffutils.create_db(
#         args.polyA_gff3_file,
#         dbfn=":memory:",
#         merge_strategy="create_unique"
#     )

#     # Generate gene pairs from the genes GFF3 file
#     result_list_gene_pairs: list[list[str, int, int, str, str, str, str]] = process_gff3(db_genes)

#     # Extract polyA sites from the polyA sites GFF3 file
#     polyA_sites: list[list[str, int, str, str]] = process_polyA_sites(db_polyA)

#     # Assign polyA sites to genes
#     matched_polyA_genes: list[list[str,str,str,str]] = assign_polyA_to_genes(result_list_gene_pairs, polyA_sites)

#     # Pretty-print the results
#     pprint.pprint(matched_polyA_genes)


# if __name__ == "__main__":
#     args = parser.parse_args()
#     main(args)



