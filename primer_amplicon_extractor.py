#!/usr/bin/env python3

import argparse
import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

from Bio import SeqIO
from Bio.Seq import Seq


# ============================================================
# BASIC FUNCTIONS
# ============================================================

def normalize_primer(primer):
    """Normalize and validate a DNA primer sequence."""

    primer = primer.upper().replace(" ", "")

    if not primer:
        raise ValueError("Primer sequence cannot be empty.")

    allowed = set("ACGT")
    invalid = set(primer) - allowed

    if invalid:
        raise ValueError(
            "Primer contains unsupported nucleotide(s): "
            + ", ".join(sorted(invalid))
        )

    return primer


def complement_base(base):
    """Return the complement of a nucleotide."""

    return str(Seq(base).complement())


# ============================================================
# GROUP IDENTIFICATION
# ============================================================

def get_group_name(description, group_mode):
    """
    Determine the group identifier.

    species mode:

        >NC_001 Aspergillus fumigatus chromosome 1

    becomes:

        Aspergillus fumigatus

    record mode:

        NC_001
    """

    parts = description.split()

    if not parts:
        return "unknown"

    if group_mode == "species":

        if len(parts) >= 3:
            return f"{parts[1]} {parts[2]}"

        if len(parts) >= 2:
            return parts[1]

    return parts[0]


def get_accession(description):
    """Return the first token of a FASTA description."""

    parts = description.split()

    if parts:
        return parts[0]

    return "unknown"


# ============================================================
# FASTA HEADER
# ============================================================

def make_header(
    description,
    start,
    end,
    group,
    reverse=False
):
    """
    Create FASTA header.

    Example:

        >PS006470.1:403842-404267 Aspergillus fruticulosus

    Reverse orientation:

        >ABC123.1:c500-200 Aspergillus fumigatus
    """

    accession = get_accession(
        description
    )

    # Convert Python 0-based coordinate to 1-based.
    start_1based = start + 1

    if reverse:

        coordinate = (
            f"c{end}-{start_1based}"
        )

    else:

        coordinate = (
            f"{start_1based}-{end}"
        )

    return (
        f"{accession}:{coordinate} {group}"
    )


# ============================================================
# PRIMER MATCHING
# ============================================================

def find_primer_matches(
    seq,
    primer,
    max_mismatch
):
    """
    Search both orientations of a primer.

    The original 3'-terminal mismatch rule is preserved.
    """

    matches = []

    variants = [
        (
            primer,
            "forward"
        ),
        (
            str(
                Seq(primer).reverse_complement()
            ),
            "reverse"
        )
    ]

    primer_length = len(
        primer
    )

    bad_3prime_pairs = {
        ("G", "A"),
        ("A", "G"),
        ("G", "G"),
        ("C", "C")
    }

    for variant, direction in variants:

        for i in range(
            len(seq) - primer_length + 1
        ):

            window = seq[
                i:i + primer_length
            ]

            mismatches = 0
            valid = True

            for j, (
                subject_base,
                primer_base
            ) in enumerate(
                zip(
                    window,
                    variant
                )
            ):

                if subject_base == primer_base:
                    continue

                is_real_3prime = (
                    (
                        direction == "forward"
                        and
                        j == primer_length - 1
                    )
                    or
                    (
                        direction == "reverse"
                        and
                        j == 0
                    )
                )

                if is_real_3prime:

                    real_primer_3_base = (
                        primer[-1]
                    )

                    if direction == "forward":

                        real_target_base = (
                            complement_base(
                                subject_base
                            )
                        )

                    else:

                        real_target_base = (
                            subject_base
                        )

                    pair = (
                        real_primer_3_base,
                        real_target_base
                    )

                    if pair in bad_3prime_pairs:

                        valid = False
                        break

                mismatches += 1

                if mismatches > max_mismatch:

                    valid = False
                    break

            if valid:

                matches.append(
                    (
                        mismatches,
                        i,
                        variant,
                        direction
                    )
                )

    return matches


# ============================================================
# PROCESS ONE FASTA RECORD
# ============================================================

def process_record(args):
    """
    Process one FASTA record.

    Correctly oriented primer-defined amplicons are returned
    BEFORE length filtering.

    This allows NO_AMPLICON and LENGTH_EXCLUDED to be
    distinguished correctly.
    """

    (
        seq,
        description,
        group,
        primer_1,
        primer_2,
        max_mismatch
    ) = args

    matches_1 = find_primer_matches(
        seq,
        primer_1,
        max_mismatch
    )

    matches_2 = find_primer_matches(
        seq,
        primer_2,
        max_mismatch
    )

    candidates = []

    if not matches_1 or not matches_2:

        return {
            "group": group,
            "candidates": []
        }

    for (
        mm1,
        start1,
        var1,
        dir1
    ) in matches_1:

        for (
            mm2,
            start2,
            var2,
            dir2
        ) in matches_2:

            end1 = (
                start1 + len(var1)
            )

            end2 = (
                start2 + len(var2)
            )

            # --------------------------------------------
            # Primer orientation
            # --------------------------------------------

            if (
                dir1 == "forward"
                and
                dir2 == "reverse"
            ):

                reverse = False

            elif (
                dir1 == "reverse"
                and
                dir2 == "forward"
            ):

                reverse = True

            else:

                # Same-direction primer pairs are ignored.
                continue

            # --------------------------------------------
            # Region between primers
            # --------------------------------------------

            inner_start = min(
                end1,
                end2
            )

            inner_end = max(
                start1,
                start2
            )

            # Primers overlap or are incorrectly positioned.
            if inner_start >= inner_end:
                continue

            inner_seq = seq[
                inner_start:inner_end
            ]

            # Normalize the extracted sequence to the
            # biological forward orientation.
            if reverse:

                inner_seq = str(
                    Seq(
                        inner_seq
                    ).reverse_complement()
                )

            inner_seq = (
                inner_seq
                .replace(" ", "")
                .replace("\t", "")
                .replace("\n", "")
                .upper()
            )

            candidates.append(
                {
                    "group": group,

                    "description": (
                        description
                    ),

                    "total_mismatch": (
                        mm1 + mm2
                    ),

                    "primer1_mismatch": mm1,
                    "primer2_mismatch": mm2,

                    "primer1_direction": dir1,
                    "primer2_direction": dir2,

                    "primer1_start": start1,
                    "primer1_end": end1,

                    "primer2_start": start2,
                    "primer2_end": end2,

                    "inner_start": (
                        inner_start
                    ),

                    "inner_end": (
                        inner_end
                    ),

                    "sequence": (
                        inner_seq
                    ),

                    "length": len(
                        inner_seq
                    ),

                    "reverse": (
                        reverse
                    )
                }
            )

    # --------------------------------------------------------
    # Remove duplicate detections of exactly the same
    # amplicon from the SAME FASTA record and coordinates.
    #
    # This does NOT collapse identical sequences occurring
    # in different chromosomes/contigs/records.
    # --------------------------------------------------------

    unique_hits = {}

    for candidate in candidates:

        key = (
            candidate["description"],
            candidate["inner_start"],
            candidate["inner_end"],
            candidate["reverse"],
            candidate["sequence"]
        )

        if key not in unique_hits:

            unique_hits[
                key
            ] = candidate

        else:

            previous = unique_hits[
                key
            ]

            if (
                candidate["total_mismatch"]
                <
                previous["total_mismatch"]
            ):

                unique_hits[
                    key
                ] = candidate

    candidates = sorted(
        unique_hits.values(),
        key=lambda x: (
            x["inner_start"],
            x["inner_end"],
            x["total_mismatch"]
        )
    )

    return {
        "group": group,
        "candidates": candidates
    }


# ============================================================
# COMMAND-LINE ARGUMENTS
# ============================================================

def parse_arguments():

    parser = argparse.ArgumentParser(
        description=(
            "Primer-based amplicon extraction "
            "with group-level analysis."
        )
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input FASTA file."
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output directory."
    )

    parser.add_argument(
        "--forward",
        required=True,
        help="Forward primer sequence (5'-3')."
    )

    parser.add_argument(
        "--reverse",
        required=True,
        help="Reverse primer sequence (5'-3')."
    )

    parser.add_argument(
        "--min-length",
        required=True,
        type=int,
        help=(
            "Minimum accepted length of the sequence "
            "between primer-binding sites."
        )
    )

    parser.add_argument(
        "--max-length",
        required=True,
        type=int,
        help=(
            "Maximum accepted length of the sequence "
            "between primer-binding sites."
        )
    )

    parser.add_argument(
        "--max-mismatch",
        required=True,
        type=int,
        help=(
            "Maximum number of mismatches allowed "
            "per primer."
        )
    )

    parser.add_argument(
        "--group-mode",
        choices=[
            "species",
            "record"
        ],
        default="species",
        help=(
            "Grouping method. "
            "'species' uses the first two words "
            "after the FASTA accession. "
            "'record' treats each FASTA record "
            "independently. "
            "Default: species."
        )
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help=(
            "Number of worker processes. "
            "Default: 1."
        )
    )

    return parser.parse_args()


# ============================================================
# MAIN
# ============================================================

def main():

    args = parse_arguments()

    # --------------------------------------------------------
    # Validate arguments
    # --------------------------------------------------------

    if not os.path.isfile(
        args.input
    ):

        raise FileNotFoundError(
            f"Input FASTA file not found: "
            f"{args.input}"
        )

    if args.min_length < 0:

        raise ValueError(
            "--min-length cannot be negative."
        )

    if (
        args.max_length
        <
        args.min_length
    ):

        raise ValueError(
            "--max-length must be greater than "
            "or equal to --min-length."
        )

    if args.max_mismatch < 0:

        raise ValueError(
            "--max-mismatch cannot be negative."
        )

    if args.threads < 1:

        raise ValueError(
            "--threads must be at least 1."
        )

    primer_1 = normalize_primer(
        args.forward
    )

    primer_2 = normalize_primer(
        args.reverse
    )

    os.makedirs(
        args.output,
        exist_ok=True
    )

    # --------------------------------------------------------
    # Output files
    # --------------------------------------------------------

    fasta_output = os.path.join(
        args.output,
        "extracted_amplicons.fasta"
    )

    summary_output = os.path.join(
        args.output,
        "group_summary.txt"
    )

    attrition_output = os.path.join(
        args.output,
        "attrition_report.txt"
    )

    # --------------------------------------------------------
    # Read input FASTA
    # --------------------------------------------------------

    tasks = []

    groups = set()

    fasta_record_count = 0

    for record in SeqIO.parse(
        args.input,
        "fasta"
    ):

        fasta_record_count += 1

        sequence = str(
            record.seq
        ).upper()

        group = get_group_name(
            record.description,
            args.group_mode
        )

        groups.add(
            group
        )

        tasks.append(
            (
                sequence,
                record.description,
                group,
                primer_1,
                primer_2,
                args.max_mismatch
            )
        )

    if not tasks:

        raise ValueError(
            "No FASTA records were found "
            "in the input file."
        )

    # --------------------------------------------------------
    # Process records
    # --------------------------------------------------------

    group_candidates = defaultdict(
        list
    )

    with ProcessPoolExecutor(
        max_workers=args.threads
    ) as executor:

        for result in executor.map(
            process_record,
            tasks
        ):

            group = result[
                "group"
            ]

            group_candidates[
                group
            ].extend(
                result[
                    "candidates"
                ]
            )

    # --------------------------------------------------------
    # Group-level evaluation
    # --------------------------------------------------------

    group_results = {}

    for group in sorted(
        groups
    ):

        all_candidates = (
            group_candidates[
                group
            ]
        )

        # Correctly oriented primer-defined amplicons
        # before the length filter.
        raw_amplicons = len(
            all_candidates
        )

        # --------------------------------------------
        # Length filter
        # --------------------------------------------

        valid_candidates = [
            candidate
            for candidate in all_candidates
            if (
                args.min_length
                <=
                candidate["length"]
                <=
                args.max_length
            )
        ]

        # --------------------------------------------
        # Count unique biological sequences
        #
        # Identical sequences from different FASTA
        # records count as ONE unique amplicon here.
        #
        # They are NOT removed from the FASTA output.
        # --------------------------------------------

        unique_sequences = set()

        for candidate in valid_candidates:

            unique_sequences.add(
                candidate[
                    "sequence"
                ]
            )

        unique_amplicons = len(
            unique_sequences
        )

        # --------------------------------------------
        # Group classification
        # --------------------------------------------

        if raw_amplicons == 0:

            status = (
                "NO_AMPLICON"
            )

        elif len(
            valid_candidates
        ) == 0:

            status = (
                "LENGTH_EXCLUDED"
            )

        elif unique_amplicons > 1:

            status = (
                "CONFLICTING"
            )

        else:

            status = (
                "INCLUDED"
            )

        group_results[
            group
        ] = {
            "raw_amplicons": (
                raw_amplicons
            ),

            "valid_amplicons": len(
                valid_candidates
            ),

            "unique_amplicons": (
                unique_amplicons
            ),

            "status": (
                status
            ),

            "valid_candidates": (
                valid_candidates
            )
        }

    # ========================================================
    # WRITE FASTA
    # ========================================================

    fasta_written = 0

    with open(
        fasta_output,
        "w"
    ) as handle:

        for group in sorted(
            group_results
        ):

            result = group_results[
                group
            ]

            # Groups without any valid-length amplicon
            # have nothing to write.
            if result["status"] in {
                "NO_AMPLICON",
                "LENGTH_EXCLUDED"
            }:

                continue

            # IMPORTANT:
            #
            # Every valid amplicon is written.
            #
            # Example:
            #
            # chromosome 1 -> AAAA
            # chromosome 2 -> AAAA
            # chromosome 3 -> AAAT
            #
            # All three entries are retained in FASTA.
            #
            # Unique amplicons = 2
            # Status = CONFLICTING

            for candidate in result[
                "valid_candidates"
            ]:

                header = make_header(
                    candidate[
                        "description"
                    ],

                    candidate[
                        "inner_start"
                    ],

                    candidate[
                        "inner_end"
                    ],

                    group,

                    candidate[
                        "reverse"
                    ]
                )

                sequence = candidate[
                    "sequence"
                ]

                handle.write(
                    f">{header}\n"
                )

                handle.write(
                    f"{sequence}\n"
                )

                fasta_written += 1

    # ========================================================
    # WRITE GROUP SUMMARY
    # ========================================================

    with open(
        summary_output,
        "w"
    ) as handle:

        header = (
            f"{'Group':<40}"
            f"{'Raw':>10}"
            f"{'Valid':>10}"
            f"{'Unique':>10}"
            f"{'Status':>20}\n"
        )

        handle.write(
            header
        )

        handle.write(
            "-" * 90
            + "\n"
        )

        for group in sorted(
            group_results
        ):

            result = group_results[
                group
            ]

            line = (
                f"{group:<40}"
                f"{result['raw_amplicons']:>10}"
                f"{result['valid_amplicons']:>10}"
                f"{result['unique_amplicons']:>10}"
                f"{result['status']:>20}\n"
            )

            handle.write(
                line
            )

    # ========================================================
    # ATTRITION COUNTS
    # ========================================================

    status_counts = {
        "NO_AMPLICON": 0,
        "LENGTH_EXCLUDED": 0,
        "CONFLICTING": 0,
        "INCLUDED": 0
    }

    for result in (
        group_results.values()
    ):

        status_counts[
            result["status"]
        ] += 1

    initial_groups = len(
        groups
    )

    # ========================================================
    # WRITE ATTRITION REPORT
    # ========================================================

    report_lines = [
        "Primer-based Amplicon Extraction",
        "Group-level Attrition Report",
        "=" * 55,
        "",
        (
            f"Input FASTA: "
            f"{args.input}"
        ),
        (
            f"Grouping mode: "
            f"{args.group_mode}"
        ),
        (
            f"Forward primer: "
            f"{primer_1}"
        ),
        (
            f"Reverse primer: "
            f"{primer_2}"
        ),
        (
            "Maximum mismatches per primer: "
            f"{args.max_mismatch}"
        ),
        (
            "Accepted internal amplicon length: "
            f"{args.min_length}-"
            f"{args.max_length} bp"
        ),
        "",
        (
            f"Input FASTA records: "
            f"{fasta_record_count}"
        ),
        (
            f"Initial groups: "
            f"{initial_groups}"
        ),
        "",
        (
            "Primer mismatch / no amplicon: "
            f"{status_counts['NO_AMPLICON']}"
        ),
        (
            "Length outside accepted range: "
            f"{status_counts['LENGTH_EXCLUDED']}"
        ),
        (
            "Multiple/conflicting amplicons: "
            f"{status_counts['CONFLICTING']}"
        ),
        (
            "Final retained: "
            f"{status_counts['INCLUDED']}"
        ),
        "",
        (
            "Groups check: "
            f"{status_counts['NO_AMPLICON']} + "
            f"{status_counts['LENGTH_EXCLUDED']} + "
            f"{status_counts['CONFLICTING']} + "
            f"{status_counts['INCLUDED']} = "
            f"{initial_groups}"
        ),
        "",
        (
            "Valid FASTA entries written: "
            f"{fasta_written}"
        )
    ]

    report_text = "\n".join(
        report_lines
    )

    with open(
        attrition_output,
        "w"
    ) as handle:

        handle.write(
            report_text
            + "\n"
        )

    # ========================================================
    # TERMINAL OUTPUT
    # ========================================================

    print(
        report_text
    )

    print(
        "\nOutput files:"
    )

    print(
        f"  {fasta_output}"
    )

    print(
        f"  {summary_output}"
    )

    print(
        f"  {attrition_output}"
    )


if __name__ == "__main__":
    main()
