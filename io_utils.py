from .utils import amp, BOLD, END

def write_probe_fasta(seqs, outfile, name=None):
    """
    Write probe sequences to a FASTA file.

    Parameters
    ----------
    seqs : dict or list
        - dict: values are probe entries [start, seq, stop, id, idx]
        - list: elements are probe entries [start, seq, stop, id, idx]
    outfile : str
        Output FASTA filename
    name : str or None
        Optional gene name prefix for FASTA headers
    """
    if not seqs:
        print("No probes to write. FASTA not created.")
        return

    # Normalize to list
    if isinstance(seqs, dict):
        probes = list(seqs.values())
    else:
        probes = seqs

    with open(outfile, "w") as f:
        for i, p in enumerate(probes):
            header = f"{name}_{i+1}" if name else str(i+1)
            f.write(f">{header}\n{p[1]}\n")

    print(f"FASTA written: {outfile}")

def print_table(columns, rows, pad=2):
    widths = [
        max(len(str(col)), max(len(str(row[i])) for row in rows))
        for i, col in enumerate(columns)]
    fmt = "".join(f"{{:<{w + pad}}}" for w in widths)
    print("\n", fmt.format(*columns))
    print("-" * sum(w + pad for w in widths))
    for row in rows:
        print(fmt.format(*row))

# Output formatting
def output(cdna, g, fullseq, count, amplifier, name, seqs):
    amplifier = amplifier.upper()
    uspc, dspc, upinit, dninit = amp(amplifier)

    if len(seqs) == 0:
        print("No probes to display.")
        return

    # Figure Layout
    print(f"{BOLD}{amplifier}_{name}{END}")
    headers = ["Pair#", "Initiator", "Spacer", "Probe", "Probe", "Spacer", "Initiator"]
    rows = []
    for i, s in enumerate(seqs, start=1):
        rows.append([i, upinit, uspc, s[1][27:52], s[1][0:25], dspc, dninit])
    print_table(headers, rows)

    # Hybridizing sequences
    headers2 = ["Pair#", "cDNAcoord", "Probe", "cDNAcoord", "cDNAcoord", "Probe", "cDNAcoord"]
    rows = []
    for i in reversed(range(len(seqs))):
        pair = i + 1
        coord1 = cdna - int(seqs[i][0])
        coord2 = coord1 - 25
        coord3 = coord2 - 2
        coord4 = cdna - int(seqs[i][2])
        rows.append([pair, coord1, seqs[i][1][0:25], coord2, coord3, seqs[i][1][27:52], coord4])
    print_table(headers2, rows)

    # Sense / anti-sense sequences
    print(f"\n{BOLD}In-place localization of probe pairs along full-length sense cDNA:{END}\n")
    print(f">{name} Sense Strand\n{g}")
    print(f"\n{BOLD}Anti-sense sequence used to create probes:{END}\n")
    print(f">{name} Anti-Sense Strand\n{fullseq}\n")