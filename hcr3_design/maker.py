from .utils import amp
from .design import limit_probes_evenly, blast_filter
from .io_utils import write_probe_fasta, output
from Bio.Seq import Seq
import numpy as np
import pandas as pd
from datetime import date

def maker(name, fullseq, amplifier, pause, polyAT, polyCG, BlastProbes, target_organism_db, background_organism_db, dropout, show, report, numbr):
    pd.set_option('display.max_columns', 500)
    pd.set_option('display.max_rows',5000)
    pd.set_option('display.width', 80)

    name = str(name)
    fullseq = str(Seq(fullseq).reverse_complement())
    cdna = len(fullseq)
    pause = int(pause)
    amplifier = str(amplifier.upper())
    uspc, dspc, upinit, dninit = amp(amplifier)

    hpA = "A" * (polyAT + 1)
    hpT = "T" * (polyAT + 1)
    hpC = "C" * (polyCG + 1)
    hpG = "G" * (polyCG + 1)

    start_limit = pause
    end_limit = cdna - pause - 52  # 52 = probe length
    table = np.vstack([np.arange(start_limit, end_limit + 1),
                       np.arange(start_limit + 52, end_limit + 1 + 52)])
    seqs = {}
    pos = []
    for a in range(len(table[0])):
        s = fullseq[table[0][a]:table[1][a]]
        if s.find(hpA) + s.find(hpT) + s.find(hpC) + s.find(hpG) <= -4:
            pos.append([table[0][a], table[1][a]])

    # Collapse consecutive positions
    newlist = []
    strt, stp = pos[0]
    newlist.append([strt, stp])
    for p in pos[1:]:
        if p[0] > (stp + 2):
            strt, stp = p
            newlist.append([strt, stp])

    newlist = np.array(newlist)
    count = str(len(newlist))

    # Prepare seqs
    for i, nl in enumerate(newlist):
        seqs[i] = [nl[0], fullseq[nl[0]:nl[0]+25] + "NN" + fullseq[nl[0]+27:nl[1]], nl[1], str(i), i]

    # BLAST Filtering
    if BlastProbes in ['y', 'Y']:
        # Automatically construct BLAST output file name from gene name and database names
        # Target BLAST output
        target_db_basename = target_organism_db.split('/')[-1].split('.')[0]  # e.g., "phikz"
        out_prefix_target = f"{name}_blast_target_{target_db_basename}"
        
        # Background BLAST output (if specified)
        if background_organism_db:
            background_db_basename = background_organism_db.split('/')[-1].replace(".fasta", "").replace(".fa","")  # e.g., "P_aeruginosa_PAO1"
            out_prefix_bg = f"{name}_blast_background_{background_db_basename}"
        else:
            out_prefix_bg = None
        
        seqs = blast_filter(seqs, name, target_organism_db, out_prefix_target, out_prefix_bg, background_organism_db=background_organism_db, target_min_cov=95.0, target_min_id=95.0, target_max_e=1e-13, background_max_cov=60.0, background_max_e=1e-12, save_outputs=True)

    # limit number of probes evenly and write to fasta
    numbr = min(numbr, len(seqs))
    seqs = limit_probes_evenly(seqs, numbr)
    count = str(len(seqs))
    write_probe_fasta(seqs, f"{name}_final_probes.fa", name=name)

    # --- Build in-place localization ---
    graphic = ['n'] * cdna
    for s in seqs:
        start, stop = s[0], s[2]
        graphic[start:stop] = str(s[1])
    g = Seq(''.join(graphic)).reverse_complement()

    output(cdna, g, fullseq, count, amplifier, name, seqs)

    # Report
    if report == 'y':
        print("Run", str(date.today()), "\nSettings:")
        print(f"\t5'Pause: {pause}")
        print(f"\tLength of acceptable polyA/polyT runs: {polyAT}")
        print(f"\tLength of acceptable polyC/polyG runs: {polyCG}")
        print(f"\tBLASTn of probes: {BlastProbes}")
        print(f"\tRemoval of probes with low quality BLAST hits: {dropout}")
        print(f"\tNumber of final probes: {count}")

    
