from time import perf_counter

t1_start = perf_counter()

with open("NZ_CP097882.1[1..499457].fa", "r") as gfg_file:
    genome = gfg_file.read()

start_codons = "ATG", "GTG", "TTG", "CTG"
stop_codons = "TAA", "TAG", "TGA"
protein_found = 0
i = 0
while i < len(genome) - 3:
    codon = genome[i:i+3]
    if codon in start_codons:
        j = i + 3
        while j < len(genome) - 3:
            stop_codon = genome[j:j+3]
            if stop_codon in stop_codons:
                # Found stop codon; extract ORF
                orf = genome[i:j+3]
                
                # Upstream region for G-rich 5-mer check
                upstream_12_start = max(0, i - 12)
                upstream_region = genome[upstream_12_start:i]

                # Check for 5-mer with at least 3 G's
                g_rich_found = False
                for k in range(len(upstream_region) - 4):  # 5-mer sliding window
                    five_mer = upstream_region[k:k+5]
                    if five_mer.count("G") >= 3:
                        g_rich_found = True
                        break

                if not g_rich_found:
                    break  # Skip this ORF, go to next start codon

                # Passed G-rich test, continue
                upstream_40_start = max(0, i - 40)
                upstream_seq = genome[upstream_40_start:i]
                protein_found += 1

                base_counts = {
                    "A": upstream_seq.count("A"),
                    "T": upstream_seq.count("T"),
                    "G": upstream_seq.count("G"),
                    "C": upstream_seq.count("C")
                }

                print(f"\nORF found from {i} to {j+3}")
                print(f"Start codon: {codon}, Stop codon: {stop_codon}")
                print(f"Protein Sequence found: {orf}")
                print(f"Upstream (40 nt): {upstream_seq}")
                print(f"Base counts: {base_counts}")
                break  # Stop at first stop codon
            j += 3
        i = j  # Move past this ORF
    else:
        i += 3

print(f"\nNumber of Proteins found: {protein_found}")

t1_stop = perf_counter()

print("Elapsed time during the whole program in seconds:",t1_stop - t1_start)
