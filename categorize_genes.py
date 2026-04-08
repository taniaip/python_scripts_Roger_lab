import argparse
import gffutils
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--genes", required=True, help="Input gene GFF3")
parser.add_argument("-f", "--fasta", required=True, help="Genome FASTA")
parser.add_argument("--output_has_stop", required=True, help="Output GFF3 for genes with stop codon")
parser.add_argument("--output_non_stop", required=True, help="Output GFF3 for genes without stop codon")
args = parser.parse_args()


def has_stop_codon(gene, sequences):
    seq = str(sequences[gene.seqid].seq[gene.start - 1:gene.end]).upper()

    if gene.strand == "+":
        return seq[-3:] in {"TAA", "TAG", "TGA"}
    elif gene.strand == "-":
        return seq[:3] in {"TTA", "CTA", "TCA"}

    return False


def write_gene_block(gene, db, out_file):
    out_file.write(str(gene) + "\n")
    for child in db.children(gene, level=None, order_by="start"):
        out_file.write(str(child) + "\n")


def main():
    db = gffutils.create_db(args.genes, dbfn=":memory:", merge_strategy="create_unique")
    sequences = SeqIO.index(args.fasta, "fasta")

    with open(args.output_has_stop, "w") as has_stop_file, \
         open(args.output_non_stop, "w") as non_stop_file:

        for gene in db.features_of_type("gene", order_by=("seqid", "start")):
            if has_stop_codon(gene, sequences):
                write_gene_block(gene, db, has_stop_file)
            else:
                write_gene_block(gene, db, non_stop_file)


if __name__ == "__main__":
    main()