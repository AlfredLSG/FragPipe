#!/usr/bin/env python3
import argparse
import sys
import logging
from collections import Counter, defaultdict
from pathlib import Path
import math

# --- Setup Logging ---
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

# --- Core Bio-functions ---

def read_fasta(fasta_file: Path) -> dict[str, str]:
    """
    Reads a FASTA file efficiently and returns a dictionary of sequences.
    Handles multi-line sequences.
    """
    logging.info(f"Reading reference genome from {fasta_file}...")
    sequences = {}
    current_seq_id = None
    current_seq = []

    with open(fasta_file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_seq_id:
                    sequences[current_seq_id] = "".join(current_seq)
                # Extract sequence ID up to the first space
                current_seq_id = line.lstrip(">").split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())
    
    if current_seq_id:
        sequences[current_seq_id] = "".join(current_seq)
    
    logging.info(f"Finished reading {len(sequences)} sequences from FASTA file.")
    return sequences

def reverse_complement(dna_seq: str) -> str:
    """Returns the reverse complement of a DNA sequence."""
    complement = str.maketrans("ACGTN", "TGCAN")
    return dna_seq.translate(complement)[::-1]

# --- Main Pipeline Functions (mimicking Perl scripts) ---

def step2_extract_end_motifs(
    fragment_bed: Path,
    motif_bed: Path,
    fasta_data: dict,
    extend_size: int = 10
):
    """
    Corresponds to 'bed_to_end_motif.pl'.
    Extracts motifs from the 5' and 3' ends of each fragment.
    """
    logging.info(f"Step 2: Extracting {extend_size}-mer end motifs...")
    records_processed = 0
    records_written = 0
    with open(fragment_bed, "r") as f_in, open(motif_bed, "w") as f_out:
        for line in f_in:
            records_processed += 1
            parts = line.strip().split("\t")
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])

            if chrom not in fasta_data:
                logging.warning(f"Chromosome '{chrom}' not found in reference FASTA. Skipping line: {line.strip()}")
                continue

            # Skip chrM if start is out of bounds (as in Perl script)
            if chrom == "chrM" and start < extend_size:
                continue

            ref_seq = fasta_data[chrom]
            
            # Extract 5' end motif (from start of fragment)
            e5_motif = ref_seq[start : start + extend_size]
            
            # Extract 3' end motif (from start of the end position)
            e3_motif_start = end - extend_size
            e3_motif = ref_seq[e3_motif_start : end]

            combined_motif = e5_motif + e3_motif
            if "N" not in combined_motif and len(e5_motif) == extend_size and len(e3_motif) == extend_size:
                f_out.write(f"{line.strip()}\t{e5_motif}\t{e3_motif}\n")
                records_written += 1

    logging.info(f"Step 2 finished. Processed {records_processed} records, wrote {records_written} to {motif_bed}.")

def step3_calculate_size_and_4mer(
    motif_bed: Path,
    output_prefix: str
):
    """
    Corresponds to 'bed_to_size_and_4mer_motif.pl'.
    Calculates fragment size distribution and 4-mer motif frequencies.
    """
    logging.info("Step 3: Calculating fragment size distribution and 4-mer frequencies...")
    
    size_counts = defaultdict(lambda: defaultdict(int)) # {chr: {len: count}}
    chr_totals = defaultdict(int)
    motif_counts = Counter()
    total_motif_sites = 0

    # The original perl script used a predefined list of 136 4-mers.
    # We generate them dynamically for robustness.
    canonical_4mers = set()
    from itertools import product
    for p in product('ACGT', repeat=4):
        mer = "".join(p)
        rev_mer = reverse_complement(mer)
        canonical_4mers.add(min(mer, rev_mer))
    MER_LIST = sorted(list(canonical_4mers))

    with open(motif_bed, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            # bed format: chr, start, end, ..., e5_motif, e3_motif
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            e5_motif, e3_motif = parts[-2], parts[-1]

            # 1. Calculate 4-mer frequencies
            motif_up = e5_motif[:4]
            motif_down = reverse_complement(e3_motif)[:4]
            motif_counts[motif_up] += 1
            motif_counts[motif_down] += 1
            total_motif_sites += 2

            # 2. Calculate size distribution
            length = end - start
            if length > 600:
                continue
            
            size_counts[chrom][length] += 1
            chr_totals[chrom] += 1
            if chrom.startswith("chr") and chrom[3:].isdigit():
                size_counts["chrA"][length] += 1
                chr_totals["chrA"] += 1
    
    # Write size file
    size_file = Path(f"{output_prefix}.size")
    with open(size_file, "w") as f_out:
        for length in range(601):
            line_data = [str(length)]
            for chr_type in ["chrA", "chrM"]:
                count = size_counts[chr_type].get(length, 0)
                total = chr_totals.get(chr_type, 1) # Avoid division by zero
                pct = (count / total) * 100
                line_data.extend([str(count), f"{pct:.6f}"])
            f_out.write("\t".join(line_data) + "\n")
    logging.info(f"Step 3a: Size distribution written to {size_file}")

    # Write 4-mer frequency file
    freq_file = Path(f"{output_prefix}.4mer.motif.freq")
    sample_name = Path(output_prefix).name
    with open(freq_file, "w") as f_out:
        f_out.write(f"Type\t{sample_name}\n")
        for mer in MER_LIST:
            count = motif_counts.get(mer, 0)
            ratio = count / total_motif_sites if total_motif_sites > 0 else 0
            f_out.write(f"{mer}\t{ratio:.10f}\n")
    logging.info(f"Step 3b: 4-mer frequencies written to {freq_file}")

def step4_calculate_size_ratio(
    size_file: Path,
    size_ratio_file: Path
):
    """
    Corresponds to 'STEP4'. Since the original script is missing, this function
    implements a common cfDNA metric: the ratio of short to long fragments.
    """
    logging.info("Step 4: Calculating short vs. long fragment size ratios...")
    # Definitions can vary, a common one is ~100-150bp vs ~151-220bp
    short_range = (100, 150)
    long_range = (151, 220)
    
    counts = {"chrA_short": 0, "chrA_long": 0, "chrM_short": 0, "chrM_long": 0}

    with open(size_file, "r") as f_in:
        for line in f_in:
            parts = line.strip().split("\t")
            length, chrA_count, _, chrM_count, _ = int(parts[0]), int(parts[1]), parts[2], int(parts[3]), parts[4]

            if short_range[0] <= length <= short_range[1]:
                counts["chrA_short"] += chrA_count
                counts["chrM_short"] += chrM_count
            elif long_range[0] <= length <= long_range[1]:
                counts["chrA_long"] += chrA_count
                counts["chrM_long"] += chrM_count

    ratio_chrA = counts["chrA_short"] / counts["chrA_long"] if counts["chrA_long"] > 0 else 0
    ratio_chrM = counts["chrM_short"] / counts["chrM_long"] if counts["chrM_long"] > 0 else 0

    with open(size_ratio_file, "w") as f_out:
        f_out.write("Metric\tValue\n")
        f_out.write(f"Autosomal_Short_Fragments\t{counts['chrA_short']}\n")
        f_out.write(f"Autosomal_Long_Fragments\t{counts['chrA_long']}\n")
        f_out.write(f"Autosomal_Short_Long_Ratio\t{ratio_chrA:.4f}\n")
        f_out.write(f"Mitochondrial_Short_Fragments\t{counts['chrM_short']}\n")
        f_out.write(f"Mitochondrial_Long_Fragments\t{counts['chrM_long']}\n")
        f_out.write(f"Mitochondrial_Short_Long_Ratio\t{ratio_chrM:.4f}\n")
    
    logging.info(f"Step 4 finished. Size ratios written to {size_ratio_file}")

def step5_calculate_mds(
    freq_file: Path,
    mds_file: Path,
    sample_label: str,
    sample_name: str,
):
    """
    Corresponds to 'calcu_motif_entropy_MDS.pl'.
    Calculates the Motif Diversity Score (MDS) using Shannon entropy.
    """
    logging.info("Step 5: Calculating Motif Diversity Score (MDS)...")
    diversity_score = 0.0
    
    with open(freq_file, "r") as f_in:
        next(f_in)  # Skip header
        for line in f_in:
            parts = line.strip().split("\t")
            freq = float(parts[1])
            
            # The formula is sum(-p * log256(p)) which is sum(-p * log(p)/log(256))
            # This avoids log(0) and correctly handles the sum.
            if freq > 0:
                motif_score = -freq * (math.log(freq) / math.log(256))
                diversity_score += motif_score

    with open(mds_file, "w") as f_out:
        f_out.write(f"{sample_label}\t{sample_name}\t{diversity_score:.10f}\n")
    
    logging.info(f"Step 5 finished. MDS score written to {mds_file}")

# --- Main Orchestrator ---

def main():
    parser = argparse.ArgumentParser(
        description="A comprehensive Python script to process a fragment BED file for motif and size analysis.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--fragment-bed",
        type=Path,
        required=True,
        help="Input fragment BED file generated from the upstream pipeline.",
    )
    parser.add_argument(
        "--ref-fasta",
        type=Path,
        required=True,
        help="Path to the reference genome FASTA file (e.g., hg38.fa).",
    )
    parser.add_argument(
        "--output-prefix",
        type=str,
        required=True,
        help="Prefix for all output files (e.g., /path/to/results/Sample1). \n"
             "This will generate Sample1.motif.bed, Sample1.size, etc.",
    )
    parser.add_argument(
        "--label",
        type=str,
        required=True,
        help="Sample label or type (e.g., Cancer, Normal) for the final MDS output.",
    )
    parser.add_argument(
        "--extend-size",
        type=int,
        default=10,
        help="The size of the motif to extract from each end. Default: 10.",
    )
    parser.add_argument(
        "--keep-intermediate",
        action="store_true",
        help="If set, keeps intermediate files (*.motif.bed, *.size, *.4mer.motif.freq).",
    )
    args = parser.parse_args()

    # --- Setup paths ---
    output_prefix = args.output_prefix
    sample_name = Path(output_prefix).name
    motif_bed_file = Path(f"{output_prefix}.motif.bed")
    size_file = Path(f"{output_prefix}.size")
    size_ratio_file = Path(f"{output_prefix}.size_ratio")
    freq_file = Path(f"{output_prefix}.4mer.motif.freq")
    mds_file = Path(f"{output_prefix}.4mer.motif.MDS")

    try:
        # --- Run Pipeline ---
        # Load reference genome once
        fasta_data = read_fasta(args.ref_fasta)

        # Step 2: Extract end motifs
        step2_extract_end_motifs(args.fragment_bed, motif_bed_file, fasta_data, args.extend_size)

        # Step 3: Calculate size and 4-mer frequencies
        step3_calculate_size_and_4mer(motif_bed_file, output_prefix)

        # Step 4: Calculate size ratio
        step4_calculate_size_ratio(size_file, size_ratio_file)

        # Step 5: Calculate MDS
        step5_calculate_mds(freq_file, mds_file, args.label, sample_name)

    finally:
        # --- Cleanup ---
        if not args.keep_intermediate:
            logging.info("Cleaning up intermediate files...")
            for f in [motif_bed_file, size_file, freq_file]:
                if f.exists():
                    f.unlink()
                    logging.info(f"Removed {f}")
    
    logging.info("Pipeline completed successfully!")
    print("\n========================================================")
    print(f"Final outputs for sample '{sample_name}':")
    print(f"  - Size Ratio: {size_ratio_file}")
    print(f"  - Motif MDS:  {mds_file}")
    print("========================================================")


if __name__ == "__main__":
    main()