# HCR3 Probe Design Package

## Overview

This package provides tools for designing HCR probe sets, filtering them via BLAST, and generating final probe sequences in FASTA and TSV formats.


## Installation

1. **Clone this GitHub repository**:

```bash
git clone https://github.com/lauraccch/hcr3_design.git
cd hcr3_design
```

2. **Create and activate a conda environment** with Python 3.13:

```bash
conda create -n HCR_design python=3.13
conda activate HCR_design
```

3. **Install the package**:

```bash
pip install -e .
```

## Usage

Import the package modules in your notebooks or scripts:

```python
from hcr3_design import maker
maker()
```


## Output Files

### FASTA Files

* **prelim_probes.fa**: All probe sequences initially designed, **before** BLAST filtering and probe number limitation.
* **probes.fa**: Probes after BLAST filtering and final number limitation.

> Both FASTA files contain only the sequences that hybridize to the target (with `NN` in the middle).
> When manually performing BLAST, use `blastn`. The sequences will show in the **opposite direction** compared to the gene.

### TSV Files

* Contain the results of BLAST searches.
* If there are multiple TSV files, this indicates probes were blasted against multiple organisms (e.g., the target organism and another present in the sample).
