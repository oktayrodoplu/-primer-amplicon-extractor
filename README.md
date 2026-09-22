# Primer Amplicon Extractor

A Python command-line tool for extracting primer-defined amplicons from nucleotide FASTA files and evaluating amplification results at the species or individual-record level.

The program searches both orientations of user-defined forward and reverse primers, allows a configurable number of mismatches, evaluates primer orientation, extracts the internal sequence between primer-binding sites, applies a user-defined length filter, and reports group-level amplicon attrition.

## Features

- Searches both orientations of forward and reverse primers
- Allows a user-defined maximum number of mismatches per primer
- Applies a specific mismatch rule at the true 3′ end of each primer
- Accepts only oppositely oriented primer pairs
- Extracts the sequence located between the primer-binding sites
- Normalizes reverse-orientation amplicons to the biological forward orientation
- Filters amplicons according to a user-defined length range
- Supports species-level or individual-record-level grouping
- Retains all valid amplicon occurrences in the FASTA output
- Detects multiple distinct amplicon sequences within the same group
- Generates a group-level summary and an attrition report
- Supports multiprocessing for large FASTA datasets

## Requirements

- Python 3
- Biopython

Install Biopython with:

```bash
pip install biopython
```

## Usage

```bash
python primer_amplicon_extractor.py \
    --input fungi.fna \
    --forward GGTAACCAAATCGGTGCTGCTTTC \
    --reverse GACCCTCAGTGTAGTGACCCTTGGC \
    --min-length 300 \
    --max-length 700 \
    --max-mismatch 3 \
    --group-mode species \
    --threads 8 \
    --output results
```

## Command-line arguments

| Argument | Description |
|---|---|
| `-i`, `--input` | Input nucleotide FASTA file |
| `-o`, `--output` | Output directory |
| `--forward` | Forward primer sequence in 5′→3′ orientation |
| `--reverse` | Reverse primer sequence in 5′→3′ orientation |
| `--min-length` | Minimum accepted length of the internal sequence between primer-binding sites |
| `--max-length` | Maximum accepted length of the internal sequence between primer-binding sites |
| `--max-mismatch` | Maximum number of mismatches allowed per primer |
| `--group-mode species` | Groups FASTA records using the first two words following the accession |
| `--group-mode record` | Treats each FASTA record independently |
| `--threads` | Number of worker processes |

## Primer matching

Each primer is searched in both its supplied orientation and its reverse-complement orientation.

Primer matches containing up to the user-defined maximum number of mismatches are accepted, subject to the implemented 3′-terminal mismatch rule.

Primer pairs are considered candidate amplicons only when the two primers occur in opposite orientations. Same-direction primer pairs are ignored.

## 3′-terminal mismatch rule

The program evaluates the true 3′ terminal nucleotide of each primer separately from the general mismatch count.

The following primer–target base combinations at the true 3′ end are treated as unacceptable:

- G–A
- A–G
- G–G
- C–C

If one of these combinations occurs at the true 3′ end, that primer match is rejected even when the total number of mismatches remains within the specified mismatch threshold.

## Amplicon extraction

For each correctly oriented primer pair, the sequence located between the two primer-binding sites is extracted.

The primer sequences themselves are not included in the extracted internal sequence.

For reverse-orientation amplicons, the extracted sequence is reverse-complemented so that all reported sequences are represented in a consistent biological forward orientation.

## Length filtering

The length filter is applied to the internal sequence between the primer-binding sites.

For example:

```text
--min-length 300
--max-length 700
```

retains internal sequences between 300 and 700 bp, inclusive.

## Group-level classification

Each species or record is assigned one of four statuses:

| Status | Meaning |
|---|---|
| `NO_AMPLICON` | No correctly oriented primer-defined amplicon was detected |
| `LENGTH_EXCLUDED` | Primer-defined amplicon(s) were detected, but none satisfied the specified length range |
| `CONFLICTING` | More than one distinct valid amplicon sequence was detected within the group |
| `INCLUDED` | Valid amplicon(s) were detected and all valid occurrences corresponded to one unique sequence |

Identical valid sequences detected in different FASTA records are retained as separate entries in the FASTA output but are counted as a single unique amplicon for group-level classification.

## Output files

The program creates three output files in the specified output directory.

### `extracted_amplicons.fasta`

Contains all amplicon occurrences that pass the length filter.

Identical sequences originating from different FASTA records are retained rather than globally deduplicated.

Example header:

```text
>PS006470.1:403842-404267 Aspergillus fruticulosus
```

Reverse-orientation hits are indicated using `c` in the coordinate field.

### `group_summary.txt`

Provides group-level counts of:

- raw primer-defined amplicons
- valid-length amplicons
- unique valid amplicon sequences
- final classification status

### `attrition_report.txt`

Provides an overall summary of the extraction process, including:

- number of input FASTA records
- number of initial groups
- groups without detectable amplicons
- groups excluded by the length criterion
- groups containing conflicting amplicons
- groups finally retained
- number of FASTA entries written

This report can be used to document dataset attrition during primer-based in silico amplicon extraction.

## Grouping

In `species` mode, the program assumes FASTA descriptions follow a structure similar to:

```text
>NC_001 Aspergillus fumigatus chromosome 1
```

and assigns the record to:

```text
Aspergillus fumigatus
```

In `record` mode, the FASTA accession is used as the group identifier.

## Notes

Species-level grouping depends on the structure of the FASTA headers. Users should verify that their input headers are compatible with the expected format before using `--group-mode species`.

The mismatch and amplicon-length parameters are user-defined and should be selected according to the primer pair, target locus, and intended analysis.

## License

This project is distributed under the MIT License. See the `LICENSE` file for details.
