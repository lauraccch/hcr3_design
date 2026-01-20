# Overview

This package provides tools for designing HCR probe sets, filtering them via BLAST, and generating final probe sequences in FASTA and TSV formats.


# Installation

1. **Clone this GitHub repository:**

```bash
git clone https://github.com/lauraccch/hcr3_design.git
cd hcr3_design
```

2. **Create and activate a conda environment:**

```bash
conda create -n HCR_design python=3.13
conda activate HCR_design
```

3. **Install the package:**

```bash
pip install -e .
```

4. **Install BLAST:**

This package relies on NCBI BLAST for probe filtering. Install BLAST+ separately by following the installtion instructions here: https://www.ncbi.nlm.nih.gov/books/NBK52640/
After installation, if necessary, add the path for the blastn binary to your path.

# Usage

1. Import the package modules in your notebooks or scripts:

```python
from hcr3_design import maker
```
2. Inputs

```python
name = "gene_name"``` -> specify the name of the target gene
- fullseq = "ATG[...]" -> enter the sequence of your target gene
- amplifier = "B2" -> select the HCR amplifier you want to use
- pause = 10 -> the number of nucleotides that should be skipped in the beginning and end of the gene (no probes designed for the first and last x nucleotides)
- polyAT = 4 -> the maximum allowed number of consecutive As or Ts
- polyCG = 4 -> the maximum allowed number of consecutive Cs or Gs
- numbr = 30 -> the maximum number of probe pairs that should be designed
- BlastProbes = "y" -> if you want the probes to be blasted ('Y') or not ('N')
- target_organism_db = "/home/user/.../complete_genome.fasta" -> path for the complete genome sequence for the organism that contains the gene that will be targeted
- background_organism_db = "/home/user /.../complete_genome.fasta" -> path for the complete genome of the second organism that will be present in the sample, but does not contain the sequence that is targeted with FISH
- dropout = "Y" -> if bad probes should be removed ('Y') or not ('N')
- report = "Y" -> if in the end a summary of used inputs etc should be given ('Y') or not ('N')

# Output Files

## FASTA Files

* **prelim_probes.fa**: All probe sequences initially designed, **before** BLAST filtering and probe number limitation.
* **probes.fa**: Final selected probes after BLAST filtering and number limitation.

> Both FASTA files contain only the sequences that will hybridize to the target (with `NN` in the middle), not the amplifiers
> When manually performing BLAST, use `blastn`. The sequences will show in the **opposite direction** compared to the gene.

## TSV Files

* Contain the detailed results of BLAST searches.
* If there are multiple TSV files, this indicates probes were blasted against multiple organisms (e.g., the target organism and another present in the sample).


## To DO
- create ready-to-order output file (oPool)
- remove show parameter
-
