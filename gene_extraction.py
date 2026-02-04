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

# Open output file for writing
with open(output_file1, "w") as out_handle:
    # Parse GenBank file for features
    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                qualifiers = feature.qualifiers

                # Extract metadata with fallback defaults
                gene = qualifiers.get("gene", ["unknown"])[0]
                product = qualifiers.get("product", ["unknown product"])[0]
                protein_id = qualifiers.get("protein_id", ["no protein_id"])[0]

                # Get start, end, and strand of CDS
                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = feature.location.strand  # 1 or -1

                # Extract CDS sequence, reverse complement if on minus strand
                cds_seq = genome_seq[start:end]
                if strand == -1:
                    cds_seq = cds_seq.reverse_complement()

                # Get start and stop codons from CDS sequence
                start_codon = cds_seq[:3]
                stop_codon = cds_seq[-3:]

                # Extract upstream 80 bp sequence (strand-aware)
                if strand == 1:
                    upstream_start80 = max(0, start - 80)
                    upstream_seq80 = genome_seq[upstream_start80:start]
                    upstream_start40 = max(0, start - 40)
                    upstream_seq40 = genome_seq[upstream_start40:start]
                else:
                    upstream_end80 = min(len(genome_seq), end + 80)
                    upstream_seq80 = genome_seq[end:upstream_end80].reverse_complement()
                    upstream_end40 = min(len(genome_seq), end + 40)
                    upstream_seq40 = genome_seq[end:upstream_end40].reverse_complement()

                # Count bases in upstream sequence
                base_counts80 = Counter(str(upstream_seq80).upper())
                base_counts_complete80 = {base: base_counts80.get(base, 0) for base in ['A', 'T', 'G', 'C']}
                base_counts40 = Counter(str(upstream_seq40).upper())
                base_counts_complete40 = {base: base_counts40.get(base, 0) for base in ['A', 'T', 'G', 'C']}

                # Write info to output file
                out_handle.write(f"Gene: {gene}\n")
                out_handle.write(f"Protein ID: {protein_id}\n")
                out_handle.write(f"Product: {product}\n")
                out_handle.write(f"Strand: {'+' if strand == 1 else '-'}\n")
                out_handle.write(f"CDS Range: {start} to {end}\n")
                out_handle.write(f"Start codon: {start_codon}, Stop codon: {stop_codon}\n")
                out_handle.write(f"CDS DNA Sequence:\n{cds_seq}\n")
                out_handle.write(f"Upstream 80 bp:\n{upstream_seq80}\n")
                out_handle.write(f"Base counts: {base_counts_complete80}\n")
                out_handle.write(f"Upstream 40 bp:\n{upstream_seq40}\n")
                out_handle.write(f"Base counts: {base_counts_complete40}\n")
                out_handle.write("-" * 60 + "\n")
with open(output_file2, "w") as f:
    # Parse GenBank file for features
    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                qualifiers = feature.qualifiers

                # Extract metadata with fallback defaults
                gene = qualifiers.get("gene", ["unknown"])[0]
                product = qualifiers.get("product", ["unknown product"])[0]
                protein_id = qualifiers.get("protein_id", ["no protein_id"])[0]

                # Get start, end, and strand of CDS
                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = feature.location.strand  # 1 or -1

                # Extract upstream 80 bp sequence (strand-aware)
                if strand == 1:
                    upstream_start80 = max(0, start - 80)
                    upstream_seq80 = genome_seq[upstream_start80:start]
                    upstream_start40 = max(0, start - 40)
                    upstream_seq40 = genome_seq[upstream_start40:start]
                else:
                    upstream_end80 = min(len(genome_seq), end + 80)
                    upstream_seq80 = genome_seq[end:upstream_end80].reverse_complement()
                    upstream_end40 = min(len(genome_seq), end + 40)
                    upstream_seq40 = genome_seq[end:upstream_end40].reverse_complement()

                # Count bases in upstream sequence
                base_counts80 = Counter(str(upstream_seq80).upper())
                base_counts_complete80 = {base: base_counts80.get(base, 0) for base in ['A', 'T', 'G', 'C']}
                base_counts40 = Counter(str(upstream_seq40).upper())
                base_counts_complete40 = {base: base_counts40.get(base, 0) for base in ['A', 'T', 'G', 'C']}
    
                f.write(f"Gene: {gene}\n")
                f.write(f"Upstream 80 bp:\n{upstream_seq80}\n")
                f.write(f"Base counts: {base_counts_complete80}\n")
                f.write(f"Upstream 40 bp:\n{upstream_seq40}\n")
                f.write(f"Base counts: {base_counts_complete40}\n")
                f.write(f"-" * 60 + "\n")


print(f"✅ Protein info saved to: {output_file1}")
print(f"✅ Protein info saved to: {output_file2}")