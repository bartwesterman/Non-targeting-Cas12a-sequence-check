# Non-targeting-Cas12a-sequence-check
_(Checking Cas12a sgRNA as non-targeting sequences using the human GRCh38 reference genome)_

Classify Cas12a guide RNA spacer sequences as exact human-genome targets, olfactory-receptor targets, or non-targeting controls by running a local short-sequence BLAST search against the human GRCh38 reference genome.

The script downloads the NCBI RefSeq **GRCh38.p14** assembly and gene annotation, builds a local BLAST database, searches both DNA strands, evaluates exact and near-exact genomic matches, annotates overlapping genes, identifies olfactory-receptor loci, and reports Cas12a PAM compatibility.

## Key features

- Searches the human **GRCh38.p14** reference assembly (`GCF_000001405.40`)
- Uses NCBI `blastn-short` for short guide sequences
- Searches **both genomic strands**
- Detects exact full-length, ungapped genomic matches
- Reports full-length matches with one or two mismatches
- Reports exact near-full-length partial alignments
- Annotates hits with genes from the corresponding NCBI GFF3 file
- Normalizes BLAST and GFF3 sequence identifiers before interval matching
- Identifies olfactory-receptor genes from gene symbols and descriptions
- Evaluates adjacent PAM sequences for wild-type AsCas12a or enAsCas12a
- Supports headered and headerless CSV/TSV input files
- Automatically installs portable NCBI command-line tools on Windows when needed
- Reuses downloaded reference files and the local BLAST database on later runs

## Classification approach

The main classification is based on **exact genomic sequence mapping**, independently of PAM compatibility.

| Classification | Meaning |
|---|---|
| `non_targeting_GRCh38_exact` | No exact full-length, ungapped match was found in GRCh38 |
| `olfactory_receptor_targeting_by_sequence` | At least one exact full-length match overlaps an annotated olfactory-receptor gene |
| `other_human_targeting_by_sequence` | Exact full-length match(es) exist, but none overlap an olfactory-receptor gene |
| `ambiguous_OR_and_other_exact_targets` | Exact matches occur in both an olfactory-receptor gene and one or more non-OR loci |
| `invalid_spacer_sequence` | The input is not an unambiguous 18–30 nucleotide A/C/G/T spacer |

PAM-compatible targetability is reported separately in the `pam_compatible_classification` column.

## Important interpretation

`non_targeting_GRCh38_exact` means that no exact full-length, ungapped match was detected in the selected GRCh38 assembly.

It does **not** prove that a guide:

- has no partial similarity to the human genome;
- has no genomic match containing mismatches;
- cannot produce mismatch-tolerant Cas12a off-target activity;
- is absent from every human haplotype, alternate assembly, transcript database, or future genome release.

The script reports one- and two-mismatch full-length BLAST hits separately, but these do not change the exact-mapping classification.

## Requirements

### Python

- Python 3.9 or newer
- `pandas`

Install the Python dependency with:

```bash
python -m pip install pandas
```

### NCBI tools

The workflow requires:

- `blastn`
- `makeblastdb`
- `blastdbcmd`
- `datasets`

On Windows, the script can automatically download portable copies of NCBI BLAST+ and NCBI Datasets into an `NCBI_portable_tools` folder.

On Linux or macOS, install the NCBI tools separately and make them available on `PATH`, or configure their locations in the script.

### Disk space

The first run downloads and extracts GRCh38 and builds a local BLAST database. At least approximately **12 GB of free disk space** is required by the current script settings. More free space is preferable.

## Repository structure

A minimal repository can use the following layout:

```text
.
├── classify_cas12a_guides_GRCh38_v6_annotation_fix.py
├── README.md
└── example_guides.csv
```

Downloaded tools, reference files, BLAST databases, and result files should normally be excluded from version control.

Suggested `.gitignore` entries:

```gitignore
NCBI_portable_tools/
GRCh38_NCBI_reference/
*.nhr
*.nin
*.nsq
*.ndb
*.not
*.ntf
*.nto
*.njs
*_queries.fasta
*_raw_blast.tsv
*_classified_v*.csv
*_exact_hits_annotated_v*.csv
*_candidate_hits_annotated_v*.csv
*_all_BLAST_hits_v*.csv
__pycache__/
```

## Input format

The input file may be comma-separated, semicolon-separated, or tab-separated.

A headered input file is recommended:

```csv
guide_id,sequence
NCT1,ACTGACTGACTGACTGACTGACT
NCT2,CTGAAGGTGTCTGGCAGAGCTTA
```

The sequence column may use names such as:

- `sequence`
- `guide_sequence`
- `spacer`
- `guide`
- `crrna`

Guide identifiers may use names such as:

- `guide_id`
- `id`
- `name`
- `guide_name`
- `label`

Headerless files are also supported. For example:

```csv
NCT1,ACTGACTGACTGACTGACTGACT
NCT2,CTGAAGGTGTCTGGCAGAGCTTA
```

The script attempts to recognize the sequence-containing column automatically. Spacer sequences must contain only unambiguous `A`, `C`, `G`, and `T` bases and must be 18–30 nucleotides long. RNA `U` characters are converted to `T` during normalization, but sequences containing ambiguous `N` bases are not submitted as valid spacers.

## Configuration

Open the script and edit the settings near the top.

### Input file

```python
INPUT_FILE = Path(
    r"C:\path\to\your\guide_sequences.csv"
)
```

All result files are written to the same folder as the input file.

### Reference assembly

The default reference is:

```python
ASSEMBLY_ACCESSION = "GCF_000001405.40"
```

This corresponds to NCBI RefSeq GRCh38.p14.

### Cas12a variant

The default configuration is enAsCas12a:

```python
CAS12A_VARIANT = "enAsCas12a"
```

Available configurations are:

```python
CAS12A_VARIANT = "enAsCas12a"
CAS12A_VARIANT = "wild_type_AsCas12a"
```

The configured PAM regular expressions are:

```python
PAM_REGEX_SETS = {
    "wild_type_AsCas12a": [
        r"TTT[ACG]",          # TTTV
    ],
    "enAsCas12a": [
        r"TT[CT][ACGT]",      # TTYN
        r"[ACG]TT[ACG]",      # VTTV
        r"T[AG]T[ACG]",       # TRTV
    ],
}
```

PAM sequences are evaluated in the 5′→3′ orientation of the protospacer strand.

### Optional executable locations

Leave these values as `None` when the executables are available on `PATH` or when using automatic Windows installation:

```python
BLAST_BIN_DIR = None
DATASETS_BIN_DIR = None
```

Alternatively, specify the executable folders:

```python
BLAST_BIN_DIR = Path(r"C:\Tools\ncbi-blast\bin")
DATASETS_BIN_DIR = Path(r"C:\Tools\ncbi-datasets")
```

### Gene-flanking region

By default, a guide must overlap the annotated gene interval:

```python
GENE_FLANK_BP = 0
```

Increase this value only when promoter or flanking-region targeting is intentionally included in the definition of a gene target.

## Usage

Run the script from a terminal or Anaconda Prompt:

```bash
python classify_cas12a_guides_GRCh38_v6_annotation_fix.py
```

Example on Windows:

```bat
conda activate oncoproai
python "C:\path\to\classify_cas12a_guides_GRCh38_v6_annotation_fix.py"
```

During the first run, the script will:

1. inspect and normalize the input file;
2. locate or install the NCBI command-line tools;
3. download the GRCh38 FASTA and GFF3 annotation;
4. build a local nucleotide BLAST database;
5. run `blastn-short` against both strands;
6. annotate exact and near-exact full-length hits;
7. evaluate PAM compatibility;
8. create summary and detailed output files.

Later runs reuse the downloaded tools, reference files, and BLAST database when available.

## BLAST settings

The script uses a sensitive short-sequence search:

```text
-task blastn-short
-strand both
-word_size 7
-dust no
-soft_masking false
-evalue 1000
-max_target_seqs 100000
```

The main exact-target classification requires:

- alignment beginning at guide position 1;
- alignment ending at the final guide nucleotide;
- alignment length equal to the complete spacer length;
- zero mismatches;
- zero gaps;
- 100% sequence identity.

Full-length ungapped hits with one or two mismatches are retained as diagnostic candidate hits.

## Output files

The script writes the following files beside the input CSV.

### Summary classification

```text
Control guide Cas12 RNAs_GRCh38_classified_v6.csv
```

Contains one row per input guide, including:

- normalized spacer sequence;
- spacer length and validity;
- main sequence-based classification;
- classification rationale;
- PAM-compatible classification;
- exact hit counts;
- plus- and minus-strand hit counts;
- olfactory-receptor hit counts;
- one- and two-mismatch hit counts;
- overlapping gene symbols;
- exact genomic loci;
- best BLAST hit summary;
- warnings for near matches.

### Exact annotated hits

```text
Control guide Cas12 RNAs_GRCh38_exact_hits_annotated_v6.csv
```

Contains exact full-length, ungapped, 100%-identity matches, including:

- genomic accession and coordinates;
- strand;
- adjacent PAM;
- PAM acceptance;
- overlapping genes;
- olfactory-receptor annotation.

### Candidate annotated hits

```text
Control guide Cas12 RNAs_GRCh38_candidate_hits_annotated_v6.csv
```

Contains full-length ungapped matches with zero, one, or two mismatches and their gene and PAM annotations.

### All BLAST hits

```text
Control guide Cas12 RNAs_GRCh38_all_BLAST_hits_v6.csv
```

Contains all BLAST alignments retained by the search, including partial hits and calculated alignment diagnostics.

### Raw BLAST output

```text
Control_guide_Cas12_RNAs_GRCh38_raw_blast.tsv
```

Contains the unmodified tabular output produced by `blastn`.

## Olfactory-receptor annotation

A genomic hit is considered an olfactory-receptor hit when it overlaps an annotated gene for which either:

- the gene symbol begins with `OR` followed by a digit, such as `OR1A2`; or
- the GFF3 gene description contains `olfactory receptor`.

BLAST accession identifiers and GFF3 sequence identifiers are normalized before overlap testing. This accommodates identifier forms such as:

```text
ref|NC_000017.11|
```

and:

```text
NC_000017.11
```

The script prints an annotation identifier audit showing how many BLAST sequence identifiers are shared with the GFF3 annotation.

## Troubleshooting

### The first guide is missing

Use the current script version, which detects headerless files and rereads them with the first row retained as data. A headered file is nevertheless recommended.

### `blastn.exe` or `makeblastdb.exe` was not found

On Windows, keep:

```python
AUTO_INSTALL_NCBI_TOOLS = True
```

The script should download portable NCBI tools automatically. Otherwise, install BLAST+ manually and set `BLAST_BIN_DIR`.

### NCBI Datasets reports a stream error

The script retries the Datasets download and can fall back to direct resumable downloads from the corresponding NCBI RefSeq assembly folder.

### BLAST fails because the Windows path contains spaces

The script stages BLAST input and database files in a space-free workspace, normally:

```text
C:\Users\Public\Cas12a_GRCh38_work
```

### A manually identified genomic hit has no gene annotation

Inspect:

```text
Control guide Cas12 RNAs_GRCh38_candidate_hits_annotated_v6.csv
```

Relevant columns include:

- `sseqid`
- `saccver`
- `annotation_seqid`
- `annotation_seqid_in_GFF`
- `target_start`
- `target_end`
- `overlapping_gene_symbols`
- `olfactory_receptor_symbols`

A genomic alignment can be real while failing gene annotation if the alignment uses an alternate contig not represented by the expected GFF identifier, lies outside the annotated gene interval, or corresponds to a partial or mismatched alignment rather than an exact full-length genomic match.

### No exact target is found, but near matches exist

Check:

- `n_full_length_1_mismatch_BLAST_hits`
- `n_full_length_2_mismatch_BLAST_hits`
- `n_near_full_length_exact_partial_hits`
- `near_match_warning`
- `best_BLAST_hit_summary`

These fields distinguish exact non-targeting classification from the absence of all genomic similarity.

## Reproducibility

The script records the following settings in its output:

- reference assembly accession;
- configured Cas12a variant;
- accepted PAM regular expressions;
- exact and mismatch hit counts;
- genomic coordinates and strands;
- gene annotations.

For reproducible analyses, report the script version, GRCh38 assembly accession, Cas12a variant, PAM definition, and date of analysis.

Suggested Methods wording:

> Cas12a spacer sequences were searched against the NCBI RefSeq GRCh38.p14 assembly (`GCF_000001405.40`) using local `blastn-short` searches of both strands. Guides lacking an exact full-length ungapped genomic match were classified as non-targeting relative to GRCh38. Exact matches were intersected with the corresponding NCBI GFF3 gene annotation to identify olfactory-receptor targets. Full-length alignments with one or two mismatches and PAM compatibility were reported separately.

## Limitations

- The analysis is reference-genome based and does not represent all human genetic variation.
- The main classification uses exact full-length sequence matching.
- One- and two-mismatch hits are reported but are not scored as functional off-targets.
- The script does not model guide activity, chromatin accessibility, nuclease concentration, mismatch position, bulges, or cell-specific genome variation.
- PAM recognition depends on the selected regular-expression definitions.
- Gene assignment is based on overlap with the downloaded NCBI GFF3 annotation.
- Alternate contigs, transcript-only records, or other genome builds may yield different results.
- This workflow should not replace a dedicated experimental or computational off-target validation strategy.

## Citation

When using this workflow, cite the relevant NCBI BLAST+, NCBI Datasets, human reference assembly, and Cas12a/enAsCas12a source publications appropriate to your analysis.

Add a project-specific citation here after the repository is archived or published.

## License

No license is implied by publication of the source code. Add an explicit license file, such as `LICENSE`, before distributing or reusing the software under defined terms.
