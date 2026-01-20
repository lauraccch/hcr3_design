from Bio.Seq import Seq
from Bio.Blast.Applications import NcbiblastnCommandline as bn
import io
import numpy as np
from .io_utils import write_probe_fasta

def limit_probes_evenly(seqs, numbr):
    n_total = len(seqs)
    if isinstance(seqs, dict):
        seqs_list = list(seqs.values())
    else:
        seqs_list = list(seqs)
    if numbr is None or numbr <= 0:
        return seqs_list
    if numbr > n_total:
        print(f"{BOLD}Requested {numbr} probes, but only {n_total} could be designed. Using all available probes.{END}")
        return seqs_list
    if numbr == n_total:
        return seqs_list
    indices = np.linspace(0, n_total - 1, numbr, dtype=int)
    return [seqs_list[i] for i in indices]


def blast_filter(seqs, name, target_organism_db, out_prefix_target, out_prefix_bg, background_organism_db=None, target_min_cov=85.0, target_min_id=85.0, target_max_e=1e-13, background_max_cov=45.0, background_max_e=1e-8, save_outputs=True):
    if target_organism_db is None:
        print("No Database for Target Organism BLAST was given. No sequences were filtered.")
        return seqs
    if background_organism_db is None:
        print("No Database for Background Organism BLAST was given. BLASTing only against Target Organism")

    # Write probe fasta
    fasta_file = write_probe_fasta(seqs, f"{name}_prelim_probes.fa", name=name)
    
    
    # --- BLAST against target ---
    filtered_seqs = {}
    cline_target = bn(query=fasta_file, subject=target_organism_db, outfmt='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs', task="blastn-short")
    stdout_target, stderr_target = cline_target()
    if save_outputs:
        header = "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqcovs\n"
        with open(f"{out_prefix_target}.tsv", "w") as f:
            f.write(header)
            f.write(stdout_target)
    
    dt = [
        (np.str_, 8),   # qseqid
        (np.str_, 40),  # sseqid
        (np.float64, 1), # pident (percent identity)
        (np.int32, 1),  # length
        (np.int32, 1),  # mismatch
        (np.int32, 1),  # gapopen
        (np.int32, 1),  # qstart
        (np.int32, 1),  # qend
        (np.int32, 1),  # sstart
        (np.int32, 1),  # send
        (np.float64, 1),# evalue
        (np.float64, 1),# bitscore
        (np.float64, 1) # qcovs
    ]
    blast_target = np.genfromtxt(io.StringIO(stdout_target), delimiter="\t", dtype=dt)

    good_target = set()
    for row in blast_target:
        qcovs = row[12]
        pident = row[2]
        evalue = row[10]
        if qcovs >= target_min_cov and pident >= target_min_id and evalue <= target_max_e:
            good_target.add(row[0])

    # --- BLAST against background ---
    bad_bg = set()
    if background_organism_db:
        print("Background Organism Specified, BLASTing in Progress.")
        cline_bg = bn(query=fasta_file, subject=background_organism_db, outfmt='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs', task="blastn-short")
        stdout_bg, stderr_bg = cline_bg()
        if save_outputs:
            header = "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqcovs\n"
            with open(f"{out_prefix_bg}.tsv", "w") as f:
                f.write(header)
                f.write(stdout_bg)
        blast_bg = np.genfromtxt(io.StringIO(stdout_bg), delimiter="\t", dtype=dt)
        for row in blast_bg:
            qcovs = row[11]
            evalue = row[10]
            if qcovs >= background_max_cov and evalue <= background_max_e:
                bad_bg.add(row[0])

    # --- Filter probes ---
    final_ids = good_target - bad_bg
    for k in seqs:
        if seqs[k][3] in final_ids:
            filtered_seqs[k] = seqs[k]

    print(f"Total probes before filtering: {len(seqs)}")
    print(f"Probes passing target & background criteria: {len(filtered_seqs)}")

    return filtered_seqs
