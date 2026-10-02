# Overview

This package provides tools for designing HCR probe sets, filtering them via BLAST and melting temperature, selecting evenly spaced probe pairs, and generating final probe sequences in FASTA and IDT-compatible XLSX formats.


# Installation

1. **Clone this GitHub repository:**

```bash
git clone https://github.com/lauraccch/hcr3_design.git
cd hcr3_design
```

2. **Create and activate a conda environment:**

```bash
conda env create -n HCR_design -f environment.yml
conda activate HCR_design
```

3. **Install the package:**

```bash
pip install -e .
```

4. **Install BLAST:**

This package relies on NCBIs BLAST+ for probe filtering. Install BLAST+ separately by following the installtion instructions here: https://www.ncbi.nlm.nih.gov/books/NBK52640/.
After installation, if necessary, add the path for the blastn binary to your path.

# Usage

1. **Import the package modules in your notebooks or scripts:**

```python
from hcr3_design import maker
```
2. **Inputs:**

```python
name = "gene_name"
```
-> specify the name of the target gene
```python
fullseq = "ATG[...]"
```
-> enter the sequence of your target gene
```python
amplifier = "B2"
```
-> select the HCR amplifier you want to use
```python
pause = 10
```
-> the number of nucleotides that should be skipped in the beginning and end of the gene (no probes designed for the first and last x nucleotides)
```python
polyAT = 4
```
-> the maximum allowed number of consecutive As or Ts
```python
polyCG = 4
```
-> the maximum allowed number of consecutive Cs or Gs
```python
numbr = 30
```
-> the maximum number of probe pairs that should be designed
```python
BlastProbes = "y"
```
-> if you want the probes to be blasted ('y') or not ('n')
```python
target_organism_db = "/home/user/.../complete_genome.fasta"
```
-> path for the complete genome sequence for the organism that contains the gene that will be targeted
```python
background_organism_db = "/home/user/.../complete_genome.fasta"
```
-> path for the complete genome of the second organism that will be present in the sample, but does not contain the sequence that is targeted with FISH
```python
dropout = "Y"
```
-> if bad probes should be removed ('Y') or not ('N')
```python
report = "Y"
```
-> if in the end a summary of used inputs etc should be given ('Y') or not ('N')
```python
min_arm_tm=30
```
-> temperature cutoff to filter out probes (°C)
```python
Na=975
```
-> Na+ concentration in mM to calculate the melting temperature
```python
formamide=30
```
-> formamide concentration in % to calculate the melting temperature
```python
dnac1=4
```
-> probe concentration in nM to calculate the melting temperature
```python
dnac2=0
```
-> target concentration in nM to calculate the melting temperature (negligible compared to the probe -> defaults to 0)

# Output Files

## FASTA Files
* **<name>_final_hyb_seqs.fa**: RNA-hybridizing regions of the BLASTed and temperature-filtered probes. The regions of the 2 splits are separated by "NN", so if you BLAST these manually again, choose the blastn algorithm for somewhat similar sequences. When looking at the graphical result, the sequences will show in the **opposite direction** compared to the gene.
* **<name>_<amplifier>_final_probes.fa**: final probes with amplifier sequences attached. If you BLAST these manually again to double check orientation etc, choose the discontiguous megablast algorithm for more dissimilar sequences

## TSV Files
Contain the detailed results of BLAST searches. If there are multiple TSV files, this indicates probes were blasted against multiple organisms (e.g., the target organism and another present in the sample).
* **_blast_target_<target_organism>.tsv**: BLAST results of the probes against the target organism (organism containing the gene)
* **_blast_background_<bg_organism>.tsv**: BLAST results of the probes against the background organism, if specified

## Excel Files
* * **<name>_<amplifier>_IDT.xlsx**: final probes in IDT ordering format


## To DO
- remove show parameter
