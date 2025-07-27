#!/usr/bin/env python3
import argparse
import sys
from collections import defaultdict
import pysam

def bam_to_fragment_bed(bam_file: str, out_file: str, min_mapq: int, max_mismatch: int):
    """
    Processes a BAM file to generate a fragment BED file.

    This script replicates the logic of two combined Perl scripts:
    1. It filters for high-quality, properly paired reads.
    2. It converts these read pairs into fragments in 4-column BED format.

    The script performs a single pass over the BAM file for efficiency,
    collecting potential read pairs in memory and processing them.

    Args:
        bam_file (str): Path to the input BAM file. Must be sorted by queryname.
        out_file (str): Path to the output BED file. Use '-' for stdout.
        min_mapq (int): Minimum mapping quality required for a read.
        max_mismatch (int): Maximum allowed NM (edit distance) tag value.
    """
    # A dictionary to hold reads waiting for their mates.
    # Key: read name (QNAME), Value: pysam AlignedSegment object
    read_pairs = defaultdict(list)
    
    # Open the input BAM file
    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
    except ValueError as e:
        print(f"Error opening BAM file: {e}", file=sys.stderr)
        print("Please ensure the BAM file is indexed.", file=sys.stderr)
        sys.exit(1)

    # Prepare output stream
    if out_file == '-':
        output = sys.stdout
    else:
        output = open(out_file, 'w')

    print(f"Processing {bam_file}...", file=sys.stderr)
    
    processed_pairs = 0
    total_reads = 0

    for read in bam:
        total_reads += 1
        
        # --- Basic Filtering (applied to each read individually) ---
        # 1. Skip if duplicate, secondary, or supplementary
        if read.is_duplicate or read.is_secondary or read.is_supplementary:
            continue
        # 2. Skip if either read or its mate is unmapped
        if read.is_unmapped or read.mate_is_unmapped:
            continue
        # 3. Check mapping quality
        if read.mapping_quality < min_mapq:
            continue
        # 4. Must be a proper pair (ensures FR orientation and reasonable insert size)
        #    This flag covers many of the original Perl script's checks.
        if not read.is_proper_pair:
            continue
        # 5. Must map to the same chromosome
        if read.reference_name != read.next_reference_name:
            continue
        # 6. Check insert size (TLEN in column 9)
        if abs(read.template_length) > 600:
            continue
        
        # Store the read, waiting for its mate
        read_pairs[read.query_name].append(read)

        # If we have found both mates, process them
        if len(read_pairs[read.query_name]) == 2:
            read1, read2 = read_pairs.pop(read.query_name)

            # --- Paired-read Filtering ---
            # 1. Check CIGAR string for both reads. Must be simple match/mismatch only.
            # The CIGAR tuple for "150M" is (0, 150), where 0 is M_OP
            if len(read1.cigartuples) != 1 or read1.cigartuples[0][0] != 0:
                continue
            if len(read2.cigartuples) != 1 or read2.cigartuples[0][0] != 0:
                continue

            # 2. Check NM tag (edit distance) for both reads
            try:
                nm1 = read1.get_tag("NM")
                nm2 = read2.get_tag("NM")
                if nm1 > max_mismatch or nm2 > max_mismatch:
                    continue
            except KeyError:
                # If NM tag is missing, we cannot check mismatch, so skip the pair
                print(f"Warning: NM tag not found for read {read.query_name}. Skipping pair.", file=sys.stderr)
                continue

            # --- Conversion to BED format ---
            # All filters passed, now calculate fragment coordinates
            
            # Fragment start is the 5' position of the leftmost read
            # pysam positions are 0-based
            fragment_start = min(read1.reference_start, read2.reference_start)
            
            # Fragment end is the 3' position of the rightmost read
            fragment_end = max(read1.reference_end, read2.reference_end)
            
            # Fragment length
            fragment_len = fragment_end - fragment_start
            
            # Write to output file in 4-column BED format
            # Use read1's chromosome, as they are on the same one.
            output.write(f"{read1.reference_name}\t{fragment_start}\t{fragment_end}\t{fragment_len}\n")
            processed_pairs += 1

    # --- Finalization ---
    print(f"Finished processing.", file=sys.stderr)
    print(f"Total reads processed: {total_reads}", file=sys.stderr)
    print(f"Fragments written to BED: {processed_pairs}", file=sys.stderr)
    if read_pairs:
        print(f"Warning: {len(read_pairs)} reads were left without a passing mate.", file=sys.stderr)

    # Close files
    bam.close()
    if out_file != '-':
        output.close()


def main():
    parser = argparse.ArgumentParser(
        description="Filter a BAM file and convert high-quality pairs to a fragment BED file.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "bam_file",
        help="Input BAM file. MUST be sorted by read name (queryname)."
    )
    parser.add_argument(
        "output_bed",
        help="Path for the output 4-column BED file. Use '-' to print to standard output."
    )
    parser.add_argument(
        "--min_mapq",
        type=int,
        default=30,
        help="Minimum mapping quality (MAPQ) for a read to be considered. [Default: 30]"
    )
    parser.add_argument(
        "--max_mismatch",
        type=int,
        default=5,
        help="Maximum number of mismatches (NM tag) allowed in each read. [Default: 5]"
    )
    args = parser.parse_args()
    
    # Important Prerequisite!
    # The Perl script implicitly worked on name-sorted data because it read the whole file.
    # This Python script requires name-sorted data to work efficiently in a single pass.
    # We will check the sort order.
    with pysam.AlignmentFile(args.bam_file, "rb") as bam:
        header = bam.header.to_dict()
        if 'HD' not in header or 'SO' not in header['HD'] or header['HD']['SO'] != 'queryname':
            print("Error: The input BAM file must be sorted by read name (queryname).", file=sys.stderr)
            print("Please sort it using: samtools sort -n -o name_sorted.bam your_input.bam", file=sys.stderr)
            sys.exit(1)

    bam_to_fragment_bed(args.bam_file, args.output_bed, args.min_mapq, args.max_mismatch)


if __name__ == "__main__":
    main()
