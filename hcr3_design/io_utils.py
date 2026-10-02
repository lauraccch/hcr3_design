import os
from .utils import amp, BOLD, END
from openpyxl import Workbook
from openpyxl.styles import Font

def write_idt_fasta(seqs, name, amplifier, upinit, uspc, dspc, dninit):
    pool_name = f"{name}_{amplifier}"
    records = []

    for i, probe in enumerate(seqs, start=1):
        raw_arm1, raw_arm2 = probe[1].split("NN")

        # Target-facing order is raw_arm2 -> raw_arm1
        left_arm = raw_arm2
        right_arm = raw_arm1

        records.append(
            f">{pool_name}_{i}_1\n"
            f"{upinit}{uspc}{left_arm}"
        )
        records.append(
            f">{pool_name}_{i}_2\n"
            f"{right_arm}{dspc}{dninit}"
        )

    output_path = f"{name}_{amplifier}_IDT.fa"

    with open(output_path, "w") as f:
        f.write("\n".join(records) + "\n")

    print(f"IDT FASTA written: {output_path}")

def write_idt_xlsx(seqs, name, amplifier, upinit, uspc, dspc, dninit):
    from openpyxl import Workbook
    from openpyxl.styles import Font

    pool_name = f"{name}_{amplifier}"

    wb = Workbook()
    ws = wb.active
    ws.title = "Sheet1"
    ws.append(["Pool name", "Sequence"])

    for probe in seqs:
        raw_arm1, raw_arm2 = probe[1].split("NN")
        left_arm = raw_arm2
        right_arm = raw_arm1

        ws.append([pool_name, f"{upinit}{uspc}{left_arm}"])
        ws.append([pool_name, f"{right_arm}{dspc}{dninit}"])

    for row in ws.iter_rows():
        for cell in row:
            cell.font = Font(name="Arial", size=10)

    ws.column_dimensions["A"].width = 23.08
    ws.column_dimensions["B"].width = 64.05

    output_path = f"{name}_{amplifier}_IDT.xlsx"
    wb.save(output_path)
    print(f"IDT Excel file written: {output_path}")

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
        return None

    probes = list(seqs.values()) if isinstance(seqs, dict) else seqs

    with open(outfile, "w") as f:
        for p in probes:
            probe_id = p[3]  # <-- critical
            f.write(f">{probe_id}\n{p[1]}\n")

    abs_path = os.path.abspath(outfile)
    print(f"FASTA written: {abs_path}")
    return abs_path

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

def print_idt_order(seqs, target_name, upinit, uspc, dspc, dninit, amplifier):
    """
    Print sequences formatted for direct IDT ordering.
    Amplifier overhangs are placed toward the inner gap of the probe pair.
    """

    full_name = f"{target_name}_{amplifier}"

    print(f"{BOLD}\nIDT ordering format:\n{END}")

    for row in seqs:
        full_probe = row[1]

        if "NN" not in full_probe:
            continue

        raw_arm1, raw_arm2 = full_probe.split("NN")

        # Target-facing order is raw_arm2 -> raw_arm1
        left_probe = raw_arm2
        right_probe = raw_arm1

        left_seq = upinit + uspc + left_probe
        right_seq = right_probe + dspc + dninit

        print(f"{full_name},{left_seq}")
        print(f"{full_name},{right_seq}")