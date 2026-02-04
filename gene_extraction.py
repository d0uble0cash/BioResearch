from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
from os import path
# Input files (your filenames)
genbank_file = "NZ_CP097882.1[1..506252].flat"
fasta_file = "NZ_CP097882.1[1..4675188].fa"
output_file1 = f"{path.splitext(fasta_file)[-2]}_protein_info.txt"
output_file2 = f"{path.splitext(fasta_file)[-2]}_upstream.fa"

# Load genome sequence from FASTA file
genome_record = SeqIO.read(fasta_file, "fasta")
genome_seq = genome_record.seq

upstream_count1 = int(input("Enter desired count1: "))
upstream_count2 = int(input("Enter desired count2: "))

def get_cds(feature):
    return int(feature.location.start), int(feature.location.end), feature.location.strand

def get_qualifiers(feature):
    return (feature.qualifiers.get("gene", ["unknown"])[0],
        feature.qualifiers.get("product", ["unknown product"])[0],
        feature.qualifiers.get("protein_id", ["no protein_id"])[0])

def write_files(files, content):
    for file in files:
        file.write(content)

# Open output file for writing
with open(output_file1, "w") as protein_info:
    with open(output_file2, "w") as upstream_fasta:
        # Parse GenBank file for features
        for record in SeqIO.parse(genbank_file, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    gene, product, protein_id = get_qualifiers(feature)

                    # Get start, end, and strand of CDS
                    start, end, strand = get_cds(feature)

                    # Extract CDS sequence, reverse complement if on minus strand
                    cds_seq = genome_seq[start:end]
                    if strand == -1:
                        cds_seq = cds_seq.reverse_complement()

                    # Get start and stop codons from CDS sequence
                    start_codon = cds_seq[:3]
                    stop_codon = cds_seq[-3:]

                    # Extract upstream 80 bp sequence (strand-aware)
                    if strand == 1:
                        upstream_start1 = max(0, start - upstream_count1)
                        upstream_seq1 = genome_seq[upstream_start1:start]
                        upstream_start2 = max(0, start - upstream_count2)
                        upstream_seq2 = genome_seq[upstream_start2:start]
                    else:
                        upstream_end1 = min(len(genome_seq), end + upstream_count1)
                        upstream_seq1 = genome_seq[end:upstream_end1].reverse_complement()
                        upstream_end2 = min(len(genome_seq), end + upstream_count2)
                        upstream_seq2 = genome_seq[end:upstream_end2].reverse_complement()

                    # Count bases in upstream sequence
                    base_counts1 = Counter(str(upstream_seq1).upper())
                    base_counts_complete1 = {base: base_counts1.get(base, 0) for base in ['A', 'T', 'G', 'C']}
                    base_counts2 = Counter(str(upstream_seq2).upper())
                    base_counts_complete2 = {base: base_counts2.get(base, 0) for base in ['A', 'T', 'G', 'C']}

                    # Write info to output file
                    write_files([protein_info, upstream_fasta], f"Gene: {gene}\n")
                    protein_info.write(f"Protein ID: {protein_id}\n")
                    protein_info.write(f"Product: {product}\n")
                    protein_info.write(f"Strand: {'+' if strand == 1 else '-'}\n")
                    protein_info.write(f"CDS Range: {start} to {end}\n")
                    protein_info.write(f"Start codon: {start_codon}, Stop codon: {stop_codon}\n")
                    protein_info.write(f"CDS DNA Sequence:\n{cds_seq}\n")
                    write_files([protein_info, upstream_fasta], f"Upstream {upstream_count1} bp:\n{upstream_seq1}\n")
                    write_files([protein_info, upstream_fasta], f"Base counts: {base_counts_complete1}\n")
                    write_files([protein_info, upstream_fasta], f"Upstream {upstream_count2} bp:\n{upstream_seq2}\n")
                    write_files([protein_info, upstream_fasta], f"Base counts: {base_counts_complete2}\n")
                    write_files([protein_info, upstream_fasta], "-" * 60 + "\n")


print(f"✅ Protein info saved to: {output_file1}")
print(f"✅ Protein info saved to: {output_file2}")